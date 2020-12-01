import pandas as pd
import numpy as np
from sgt.core.db import Bedpe, Vcf
from typing import (
    Set,
    Union,
)
from sgt._typing import (
    IntOrStr,
)

def get_id_by_slicing_info(
    bedpe_or_vcf: Union[Bedpe, Vcf],
    info: str, 
    value_idx: int = 0,
    st = None, 
    en = None,
    svtype: str = 'any') -> Set[IntOrStr]:
    """
    get_id_by_slicing_info(bedpe_or_vcf, info, value_idx, st, en, svtype)
    Return id of bedpe_or_vcf within the range [st, en).
    """
    if svtype != 'any':
        bedpe_or_vcf = bedpe_or_vcf.filter('svtype == {}'.format(svtype))
    df_info = bedpe_or_vcf.get_table(info)
    if df_info.empty:
        return set()
    df_info = df_info.loc[df_info['value_idx'] == value_idx]
    if st is None:
        if en is None:
            return set(df_info['id'].tolist())
        condition = df_info[info] < en
    elif en is None:
        condition = df_info[info] >= st
    else:
        condition = (df_info[info] >= st) & (df_info[info] < en)
    df_info = df_info.loc[condition]
    if df_info.empty:
        return set()
    return set(df_info['id'].tolist())

def get_id_by_boolean_info(
    bedpe_or_vcf: Union[Bedpe, Vcf],
    info: str,
    true_or_false: bool = True,
    svtype: str = 'any') -> Set[str]:
    if svtype != 'any':
        bedpe_or_vcf = bedpe_or_vcf.filter('svtype == {}'.format(svtype))
    df_info = bedpe_or_vcf.get_table(info)
    true_ids = set(df_info['id'].tolist())
    if true_or_false:
        return true_ids
    else:
        all_ids = bedpe_or_vcf.ids
        return all_ids - true_ids
