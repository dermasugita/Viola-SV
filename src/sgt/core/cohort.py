import pandas as pd
from sgt.core.db import Vcf, Bedpe
from typing import (
    List,
    Optional,
)

class MultiBedpe(Bedpe):
    """
    A database-like object that contains information of multiple BEDPE files.
    In this class, main keys in most tables are "global id" instead of using
    "SV id" from SV callers. "global id" is unique ID of all the SV record
    across all the samples.
    """
    _internal_attrs = [
        "_df_id",
        "_df_svpos",
        "_dict_df_info",
        "_ls_infokeys",
        "_dict_alltables",
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
        ls_bedpe: List[Bedpe] = None, 
        ls_patient_names: List[str] = None, 
        direct_tables: Optional[List[pd.DataFrame]] = None
        ):
        if direct_tables is None:
            df_id, df_svpos, dict_df_info = self.__init__from_ls_bedpe(ls_bedpe, ls_patient_names)
            self.__init__common(df_id, df_svpos, dict_df_info)
        else:
            self.__init__common(*direct_tables)
    
    def __init__from_ls_bedpe(self, ls_bedpe, ls_patient_names):
        ls_df_id = []
        ls_df_svpos = []
        dict_ls_df_info = dict() 
        for bedpe, patient_name in zip(ls_bedpe, ls_patient_names):
            df_svpos = bedpe.get_table('positions')
            df_id = df_svpos[['id']].copy()
            df_id['patients'] = patient_name
            df_id['global_id'] = str(patient_name) + '_' + df_id['id'].astype(str)
            df_id = df_id[['global_id', 'patients', 'id']]
            ls_df_id.append(df_id)

            df_svpos['id'] = str(patient_name) + '_' + df_svpos['id'].astype(str)
            ls_df_svpos.append(df_svpos)

            for key, value in bedpe._dict_df_info.items():
                value = value.copy()
                value['id'] = str(patient_name) + '_' + value.astype(str)
                if dict_ls_df_info.get(key) is None:
                    dict_ls_df_info[key] = [value]
                else:
                    dict_ls_df_info[key].append(value)
        df_concat_id = pd.concat(ls_df_id, ignore_index=True)
        df_concat_svpos = pd.concat(ls_df_svpos, ignore_index=True)
        dict_df_info = dict()
        for key, value in dict_ls_df_info.items():
            dict_df_info[key] = pd.concat(value)
        
        return (df_concat_id, df_concat_svpos, dict_df_info)
    
    def __init__common(self, df_id, df_svpos, dict_df_info):
        self._df_id = df_id
        self._df_svpos = df_svpos 
        self._dict_df_info = dict_df_info
        self._ls_infokeys = [x.lower() for x in dict_df_info.keys()]
        ls_keys = ['global_id', 'positions'] + self._ls_infokeys
        ls_values = [df_id, df_svpos] + list(dict_df_info.values())
        self._dict_alltables = {k: v for k, v in zip(ls_keys, ls_values)}
        self._repr_config = {
            'info': None,
        }
        

    def filter_by_id(self, arrlike_id):
        """
        filter_by_id(arrlike_id)
        Filter MultiBedpe object according to the list of SV ids.
        Return object is also an instance of the MultiBedpe object

        Parameters
        ---------------
        arrlike_id: list-like
            Global ids which you would like to keep.
        
        Returns
        ---------------
        MultiBedpe
            A MultiBedpe object with the SV id specified in the arrlike_id argument.
            All records associated with SV ids that are not in the arrlike_id will be discarded.
        
        """
        df_global_id = self.get_table('global_id')
        out_global_id = df_global_id.loc[df_global_id['global_id'].isin(arrlike_id)].reset_index(drop=True)
        out_svpos = self._filter_by_id('positions', arrlike_id)
        out_dict_df_info = {k: self._filter_by_id(k, arrlike_id) for k in self._ls_infokeys}
        return MultiBedpe(direct_tables=[out_global_id, out_svpos, out_dict_df_info])
    

    def classify_manual_svtype(self, ls_conditions, ls_names, ls_order=None, return_data_frame=True):
        """
        classify_manual_svtype(ls_conditions, ls_names, ls_order=None)
        Classify SV records by user-defined criteria. A new INFO table named
        'manual_sv_type' will be created.
        """
        set_ids_current = set(self.ids)
        obj = self
        ls_ids = []
        ls_result_names = []
        for func, name in zip(ls_conditions, ls_names):
            obj = obj.filter_by_id(set_ids_current)
            set_ids = func(obj)
            set_ids_intersection = set_ids_current & set_ids
            ls_ids += list(set_ids_intersection)
            ls_result_names += [name for i in range(len(set_ids_intersection))]
            set_ids_current = set_ids_current - set_ids_intersection
        ls_ids += list(set_ids_current)
        ls_result_names += ['others' for i in range(len(set_ids_current))]
        ls_zeros = [0 for i in range(len(self.ids))]
        df_result = pd.DataFrame({'id': ls_ids, 'value_idx': ls_zeros, 'manual_sv_type': ls_result_names})
        self.create_info_table('manual_sv_type', df_result)
        if return_data_frame:
            if ls_order is None:
                pd_ind_reindex = pd.Index(ls_names + ['others'])
            else:
                pd_ind_reindex = pd.Index(ls_order)
            df_feature_counts = self.get_feature_count_as_data_frame(ls_order=pd_ind_reindex)
            return df_feature_counts
    
    def get_feature_count_as_data_frame(self, feature='manual_sv_type', ls_order=None):
        df_feature = self.get_table(feature)
        df_id = self.get_table('global_id')
        df_merged = pd.merge(df_feature, df_id, left_on='id', right_on='global_id')
        df_feature_counts = df_merged.pivot_table('global_id', index='patients', columns=feature, aggfunc='count', fill_value=0)
        if ls_order is not None:
            pd_ind_reindex = pd.Index(ls_order, name=feature)
            df_feature_counts = df_feature_counts.reindex(columns=pd_ind_reindex, fill_value=0)
        return df_feature_counts