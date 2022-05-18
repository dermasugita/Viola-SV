import pandas as pd
import viola
from intervaltree import IntervalTree, Interval
from viola.core.vcf import Vcf
from viola.core.bedpe import Bedpe
from viola.core.cohort import MultiVcf
from collections import OrderedDict
from typing import (
    List,
    Optional
)
def merge(ls_inputs, ls_caller_names=None, mode='distance', threshold=100, integration=True):
    """
    merge(ls_inputs:list, ls_caller_names:list, threshold=100, integration=False)
    merge Vcf objects or Bedpe objects
    Return a merged Vcf or Bedpe

    Parameters
    ----------
    ls_inputs:list
        A list of Vcf or Bedpe objects to be merged
    ls_caller_names:list
        A list of names of callers(str).
        Only needed for Bedpes.
    mode: {'distance', 'confidence_intervals'}, default 'distance'
        The mode of the merging strategy. 

        * ``'distance'``: Merge SV records by representative SV positions, that is, coordinates of POS field or that of END in the INFO field. If multiple SV positions are within the distance specified in ``threshold`` each other, they will be merged.
        * ``'confidence_intervals'``: Merge SV records according to the confidence intervals reported by SV callers. If confidence intervals of multiple SV records share their genomic coordinates at least 1bp, the will be merged.

    threshold:int, default 100
        Two SVs of mutual distance is under 
        this threshold are cosidered to be identical.
    integration:bool, default True
        If true, only the results from one caller will be retained  for each SV event that is predicted to be identical, according to the caller priority. 
        Otherwise, all SV records in the input VCF file are retained and the "mergedid" INFO is added.
        The caller priority of integration is the same as 
        the order of callers in ls_inputs.
        For now it only works for Vcf.
        Detail explanation is described in the User Guide :ref:`"VCF Merging"<merge>`.
        
    Returns
    ----------
    Vcf or Bedpe

    See Also
    ----------
    :ref:`VCF Merging<merge>`
    """
    first_object = ls_inputs[0]
    
    if isinstance(first_object, Vcf):
        return first_object.merge(ls_vcf = ls_inputs, mode=mode, integration = integration, threshold = threshold)

    if isinstance(first_object, Bedpe):
        return first_object.merge(ls_bedpe = ls_inputs, ls_caller_names = ls_caller_names, mode=mode, threshold = threshold)        

class TmpVcfForMerge(MultiVcf):
    _internal_attrs = [
        "_df_id",
        "_df_patients",
        "_df_svpos",
        "_odict_df_info",
        "_ls_infokeys",
        "_ls_patients",
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
            df_id, df_patients, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers  = self.__init__from_ls_vcf(ls_vcf, ls_patient_names)
            self.__init__common(df_id, df_patients, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers)
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
                
                # resolving conflict of 'number' field in the headers of different callers.
                df.loc[df['id'] == 'HOMLEN', 'number'] = 1
                df.loc[df['id'] == 'PE', 'number'] = 1
                df.loc[df['id'] == 'SR', 'number'] = 1
                # /resolving conflict of 'number' field in the headers of different callers.

                if idx == 0:
                    df_merged = df
                    continue
                on = ['id', 'number', 'type']
                df_merged = df_merged.merge(df, how='outer', on=on)
            
            mask = df_merged['id'].duplicated(keep=False)
            
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
                if key == 'contigs_meta':
                    if len(ls_patient_names) > 1:
                        lumpy_index = [i for i, x in enumerate(ls_patient_names) if x == 'lumpy']
                        if len(lumpy_index) == 1:
                            lumpy_index = lumpy_index[0]
                            value = [v for i, v in enumerate(value) if i != lumpy_index]
                for idx, df in enumerate(value):
                    if 'number' in df.columns:
                        df['number'] = df['number'].astype(object)
                    if idx == 0:
                        df_merged = df
                        continue
                    on = list(df_merged.columns)
                    df_merged = df_merged.merge(df, how='outer', on=on)
                    if 'number' in on:
                        df_merged['number'] = df_merged['number'].astype(object)
            
            odict_df_headers[key] = df_merged
        

        # /Header Integration

        ls_patient_id = [i for i in range(len(ls_patient_names))]
        df_patients = pd.DataFrame({'id': ls_patient_id, 'patients': ls_patient_names})
        for vcf, patient_id, patient_name in zip(ls_vcf, ls_patient_id, ls_patient_names):
            df_svpos = vcf.get_table('positions')
            df_filters = vcf.get_table('filters')
            df_formats = vcf.get_table('formats')
            df_id = df_svpos[['id']].copy()
            df_id['patient_id'] = patient_id
            df_id['global_id'] = str(patient_name) + '_' + df_id['id'].astype(str)
            df_id = df_id[['global_id', 'patient_id', 'id']]
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
                if odict_ls_df_info.get(info, None) is None:
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
        
        
        return (df_concat_id, df_patients, df_concat_svpos, df_concat_filters, odict_df_info, df_concat_formats, odict_df_headers)

    def __init__common(self, df_id, df_patients, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers = {}):
        self._df_id = df_id
        self._df_patients = df_patients
        self._df_svpos = df_svpos
        self._df_filters = df_filters
        self._odict_df_info = odict_df_info
        self._df_formats = df_formats
        self._odict_df_headers = odict_df_headers
        self._ls_patients = df_patients['patients'].to_list()
        self._ls_infokeys = [x.lower() for x in list(odict_df_info.keys())]
        ls_keys = ['global_id', 'patients', 'positions', 'filters'] + self._ls_infokeys + ['formats'] + \
        list(odict_df_headers.keys())
        ls_values = [df_id, df_patients, df_svpos, df_filters] + list(odict_df_info.values()) + [df_formats] + list(odict_df_headers.values())
        self._odict_alltables = OrderedDict([(k, v) for k, v in zip(ls_keys, ls_values)])
        self._repr_config = {
            'info': None,
        }

class IntervalTreeForMerge():
    def __init__(self, df, header):
        self._df = df
        self._header = header
        self._tree = self._intervaltree_init(df)
    
    def _intervaltree_init(self, df):
        dict_tree = dict()
        for row in df.itertuples():
            chrom = getattr(row, 'chrom')
            st = getattr(row, 'chromStart')
            en = getattr(row, 'chromEnd') + 1
            if dict_tree.get(chrom) is None:
                dict_tree[chrom] = IntervalTree()
            dict_tree[chrom][st:en] = row
        
        return dict_tree
    
    def query(self, chrom, st, en=None) -> pd.DataFrame:
        if self._tree.get(chrom) is None:
            return pd.DataFrame(columns = self._df.columns)
        tree = self._tree[chrom]
        if en is None:
            set_out = tree[st]
        else:
            set_out = tree[st:en]
        if len(set_out) == 0:
            return pd.DataFrame(columns = self._df.columns)
        ls_out = list(set_out)
        ls_ser_out = [x.data for x in ls_out]
        df_out = pd.DataFrame(ls_ser_out)
        df_out.set_index('Index', inplace=True)
        return df_out
