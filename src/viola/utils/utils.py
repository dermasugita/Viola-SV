import pandas as pd
import numpy as np
import re
from typing import (
    Set,
    Union,
)
from viola._typing import (
    IntOrStr,
)

def is_url(x):
    regex = re.compile(
            r'^(?:http|ftp)s?://' # http:// or https://
            r'(?:(?:[A-Z0-9](?:[A-Z0-9-]{0,61}[A-Z0-9])?\.)+(?:[A-Z]{2,6}\.?|[A-Z0-9-]{2,}\.?)|' #domain...
            r'localhost|' #localhost...
            r'\d{1,3}\.\d{1,3}\.\d{1,3}\.\d{1,3})' # ...or ip
            r'(?::\d+)?' # optional port
            r'(?:/?|[/?]\S+)$', re.IGNORECASE)
    return re.match(regex, x) is not None

def get_inslen_and_insseq_from_alt(alt):
    """
    parse_alt(str)
    
    This function is the modification of Reader._parse_alt of PyVCF package.
    """
    alt_pattern = re.compile('[\[\]]')
    if alt_pattern.search(alt) is not None:
        # Paired breakend
        items = alt_pattern.split(alt)
        mate_be = items[1].split(':')
        chrom2 = mate_be[0]
        if chrom2[0] == '<':
            chrom2 = chrom2[1:-1]
            within_main_assembly = False
        else:
            within_main_assembly = True
        pos2 = mate_be[1]
        orientation = (alt[0] == '[' or alt[0] == ']') # If True, strand of be1 is '-'.
        remote_orientation = (re.search('\[', alt) is not None)
        if orientation:
            insertion_sequence = items[2][:-1]
        else:
            insertion_sequence = items[0][1:]
        #print(connecting_sequence)
        return (len(insertion_sequence), insertion_sequence)
    else:
        return (0, '')


def get_id_by_slicing_info(
    bedpe_or_vcf: "Union[Bedpe, Vcf]",
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
    bedpe_or_vcf: "Union[Bedpe, MultiBedpe, Vcf]",
    info: str,
    true_or_false: bool = True,
    svtype: str = 'any') -> Set[str]:
    """
    get_id_by_boolean_info(bedpe_or_vcf, info, true_or_false, svtype)
    Return SV id that specified
    """
    if svtype != 'any':
        bedpe_or_vcf = bedpe_or_vcf.filter('svtype == {}'.format(svtype))
    df_info = bedpe_or_vcf.get_table(info)
    true_ids = set(df_info['id'].tolist())
    if true_or_false:
        return true_ids
    else:
        all_ids = bedpe_or_vcf.ids
        return all_ids - true_ids
