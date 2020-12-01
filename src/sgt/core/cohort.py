import pandas as pd
from sgt.core.db import Vcf, Bedpe
from typing import (
    List,
)

class MultiBedpe(Bedpe):
    def __init__(self, ls_bedpe: List[Bedpe], ls_patient_names: List[str]):
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
        self._df_id = pd.concat(ls_df_id)
        self._df_svpos = pd.concat(ls_df_svpos)
        dict_df_info = dict()
        for key, value in dict_ls_df_info.items():
            dict_df_info[key] = pd.concat(value)
        self._dict_df_info = dict_df_info
        self._ls_infokeys = [x.lower() for x in dict_df_info.keys()]
        ls_keys = ['global_id', 'positions'] + self._ls_infokeys
        ls_values = [self._df_id, df_svpos] + list(dict_df_info.values())
        self._dict_alltables = {k: v for k, v in zip(ls_keys, ls_values)}
        self._repr_config = {
            'info': None,
        }