import viola
from viola.core.bedpe import Bedpe
from viola.core.vcf import Vcf
import pandas as pd


def assert_bedpe_equal(left, right, check_like=True):
    ls_left_table = left.table_list
    ls_right_table = right.table_list
    assert len(ls_left_table) == len(ls_right_table)
    assert set(ls_left_table) == set(ls_right_table)
    for tablename in ls_left_table:
        df_left = left.get_table(tablename)
        df_right = right.get_table(tablename)
        id_left = df_left.columns[0]
        id_right = df_right.columns[0]
        df_left.sort_values(id_left, inplace=True, ignore_index=True)
        df_right.sort_values(id_right, inplace=True, ignore_index=True)
        if df_left.empty & df_right.empty:
            continue
        try:
            pd.testing.assert_frame_equal(df_left, df_right, check_like=check_like, check_exact=True)
        except AssertionError:
            print('\nwhen asserting {} table, following error occured!'.format(tablename))
            raise 

    # assert patient name
    if left._patient_name is None:
        assert right._patient_name is None
    elif right._patient_name is None:
        assert False
    else:
        assert left._patient_name == right._patient_name
    
        

def assert_vcf_equal(left, right, check_like=True):
    ls_left_table = left.table_list
    ls_right_table = right.table_list
    assert len(ls_left_table) == len(ls_right_table)
    assert set(ls_left_table) == set(ls_right_table)
    for tablename in ls_left_table:
        df_left = left.get_table(tablename)
        df_right = right.get_table(tablename)
        id_left = df_left.columns[0]
        id_right = df_right.columns[0]
        df_left.sort_values(id_left, inplace=True, ignore_index=True)
        df_right.sort_values(id_right, inplace=True, ignore_index=True)
        if df_left.empty & df_right.empty:
            continue
        try:
            pd.testing.assert_frame_equal(df_left, df_right, check_like=check_like, check_exact=True)
        except AssertionError:
            print('\nwhen asserting {} table, following error occured!'.format(tablename))
            raise 
    
    # assert metadata
    assert left._metadata == right._metadata

    # assert patient name
    if left._patient_name is None:
        assert right._patient_name is None
    elif right._patient_name is None:
        assert False
    else:
        assert left._patient_name == right._patient_name
