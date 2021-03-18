import pandas as pd
import viola
from viola.core.db import Vcf
from viola.core.cohort import MultiVcf
from collections import OrderedDict
from typing import (
    List,
    Optional
)
class TmpVcfForMerge(MultiVcf):
    _internal_attrs = [
        "_df_id",
        "_df_svpos",
        "_odict_df_info",
        "_ls_infokeys",
        "_odict_alltables",
        "_repr_config",
        "_sig_criteria"
    ]
    _internal_attrs_set = set(_internal_attrs)
    _repr_column_names = [
        "id",
        "bp1",
        "bp2",
        "strand",
        "qual",
        "svtype",
    ]
    _repr_column_names_set = set(_repr_column_names)
    def __init__(
        self,
        ls_vcf: List[Vcf] = None, 
        ls_patient_names: List[str] = None, 
        direct_tables: Optional[List[pd.DataFrame]] = None
        ):
        if direct_tables is None:
            df_id, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers  = self.__init__from_ls_vcf(ls_vcf, ls_patient_names)
            self.__init__common(df_id, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers)
        else:
            self.__init__common(*direct_tables)

    def __init__from_ls_vcf(self, ls_vcf, ls_patient_names):
        ls_df_id = []
        ls_df_svpos = []
        ls_df_filters = []
        odict_ls_df_info = OrderedDict() 
        ls_df_formats = []
        odict_ls_df_headers = OrderedDict()

        # Header Integration
        def _integrate_infos_meta(value):
            for idx, df in enumerate(value):
                if idx == 0:
                    df_merged = df
                    continue
                on = ['id', 'number', 'type']
                df_merged = df_merged.merge(df, how='outer', on=on)
            
            # construct description
            ser_description = ''
            for col_name in df_merged.columns:
                if col_name.startswith('description_'):
                    ser_description = ser_description + df_merged[col_name].astype(str) + ', '
                    df_merged.drop(col_name, axis=1, inplace=True)
            df_merged['description'] = ser_description
            df_merged['source'] = None
            df_merged['version'] = None
            return df_merged

        for vcf, patient_name in zip(ls_vcf, ls_patient_names):
            for key, value in vcf._odict_df_headers.items():
                value = value.copy()
                if key == 'infos_meta':
                    description_name = 'description_' + patient_name
                    value[description_name] = patient_name + ': ' + value['description'].astype(str)
                    value.drop(['description', 'source', 'version'], axis=1, inplace=True)
                if odict_ls_df_headers.get(key) is None:
                    odict_ls_df_headers[key] = [value]
                else:
                    odict_ls_df_headers[key].append(value)

        odict_df_headers = OrderedDict()
        for key, value in odict_ls_df_headers.items():
            if key == 'infos_meta':
                df_merged = _integrate_infos_meta(value)
            else:
                for idx, df in enumerate(value):
                    if idx == 0:
                        df_merged = df
                        continue
                    on = list(df_merged.columns)
                    df_merged = df_merged.merge(df, how='outer', on=on)
            
            odict_df_headers[key] = df_merged
        

        # /Header Integration

        for vcf, patient_name in zip(ls_vcf, ls_patient_names):
            df_svpos = vcf.get_table('positions')
            df_filters = vcf.get_table('filters')
            df_formats = vcf.get_table('formats')
            df_id = df_svpos[['id']].copy()
            df_id['patients'] = patient_name
            df_id['global_id'] = str(patient_name) + '_' + df_id['id'].astype(str)
            df_id = df_id[['global_id', 'patients', 'id']]
            ls_df_id.append(df_id)

            df_svpos['id'] = str(patient_name) + '_' + df_svpos['id'].astype(str)
            ls_df_svpos.append(df_svpos)

            df_filters['id'] = str(patient_name) + '_' + df_filters['id'].astype(str)
            ls_df_filters.append(df_filters)

            df_formats['id'] = str(patient_name) + '_' + df_formats['id'].astype(str)
            ls_df_formats.append(df_formats)

            for info in odict_df_headers['infos_meta'].id:
                df_info_ = vcf._odict_df_info.get(info, None)
                if df_info_ is None:
                    df_info = pd.DataFrame(columns=('id', 'value_idx', info.lower()))
                else:
                    df_info = df_info_.copy()
                    df_info['id'] = str(patient_name) + '_' + df_info['id'].astype(str)
                if odict_ls_df_info.get(info) is None:
                    odict_ls_df_info[info] = [df_info]
                else:
                    odict_ls_df_info[info].append(df_info)

            
        df_concat_id = pd.concat(ls_df_id, ignore_index=True)
        df_concat_svpos = pd.concat(ls_df_svpos, ignore_index=True)
        df_concat_filters = pd.concat(ls_df_filters, ignore_index=True)
        df_concat_formats = pd.concat(ls_df_formats, ignore_index=True)
        odict_df_info = OrderedDict()

        for key, value in odict_ls_df_info.items():
            odict_df_info[key] = pd.concat(value)
        
        
        return (df_concat_id, df_concat_svpos, df_concat_filters, odict_df_info, df_concat_formats, odict_df_headers)

    def __init__common(self, df_id, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers = {}):
        self._df_id = df_id
        self._df_svpos = df_svpos
        self._df_filters = df_filters
        self._odict_df_info = odict_df_info
        self._df_formats = df_formats
        self._odict_df_headers = odict_df_headers
        self._ls_infokeys = list(odict_df_info.keys())
        ls_keys = ['global_id', 'positions', 'filters'] + self._ls_infokeys + ['formats'] + \
        list(odict_df_headers.keys())
        ls_values = [df_id, df_svpos, df_filters] + list(odict_df_info.values()) + [df_formats] + list(odict_df_headers.values())
        self._odict_alltables = OrderedDict([(k, v) for k, v in zip(ls_keys, ls_values)])
        self._repr_config = {
            'info': None,
        }
