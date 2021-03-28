import viola
import pandas as pd
from io import StringIO
import sys, os
HERE = os.path.abspath(os.path.dirname(__file__))
data_expected = """vcf1_test1	0	small_del
vcf2_test1	0	small_del
vcf1_test2	0	small_del
vcf2_test2	0	small_del
vcf1_test3	0	large_del
vcf2_test3	0	large_del
vcf1_test4	0	large_del
vcf2_test4	0	large_del
vcf1_test5	0	large_del
vcf2_test5	0	large_del
vcf1_test6	0	small_dup
vcf2_test6	0	small_dup
vcf1_test7	0	small_inv
vcf2_test7	0	small_inv
vcf1_test8	0	others
vcf2_test8	0	others
vcf1_test9	0	small_inv
vcf2_test9	0	small_inv
vcf1_viola_breakpoint:0	0	tra
vcf2_viola_breakpoint:0	0	tra
vcf1_viola_breakpoint:1	0	tra
vcf2_viola_breakpoint:1	0	tra
"""

DEFINITIONS = """name 'small_del'
0 SVLEN > -100
1 SVTYPE == DEL
logic 0 & 1

name 'large_del'
0 SVTYPE == DEL
logic 0

name 'small_dup'
0 SVLEN < 100
1 SVTYPE == DUP
logic 0 & 1

name 'large_dup'
0 SVTYPE == DUP
logic 0

name 'small_inv'
0 SVLEN < 100
1 SVTYPE == INV
logic 0 & 1

name 'tra'
0 SVTYPE == TRA
logic 0
"""

def small_del(x):
    return x.filter(['svlen > -100', 'svtype == DEL']).ids
def large_del(x):
    return x.filter(['svtype == DEL']).ids
def small_dup(x):
    return x.filter(['svlen < 100', 'svtype == DUP']).ids
def large_dup(x):
    return x.filter(['svtype == DUP']).ids
def small_inv(x):
    return x.filter(['svlen < 100', 'svtype == INV']).ids
def tra(x):
    return x.filter('svtype == TRA').ids

def test_classify_manual_svtype():
    vcf = viola.read_vcf(os.path.join(HERE, 'data/manta1.vcf'))
    vcf2 =  vcf.copy()
    vcf = vcf.breakend2breakpoint()
    vcf2 = vcf2.breakend2breakpoint()
    multi_vcf = viola.MultiVcf([vcf, vcf2], ['vcf1', 'vcf2'])
    ls_conditions = [small_del, large_del, small_dup, large_dup, small_inv, tra]
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    result = multi_vcf.classify_manual_svtype(ls_conditions=ls_conditions, ls_names=ls_names)

    manual_sv_type = multi_vcf.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.DataFrame([[2, 3, 1, 0, 2, 2, 1], [2, 3, 1, 0, 2, 2, 1]])
    result_expected.columns = ls_names + ['others']
    result_expected.columns.name = 'manual_sv_type'
    result_expected.index = ['vcf1', 'vcf2']
    result_expected.index.name = 'patients'
    pd.testing.assert_frame_equal(result, result_expected)

"""
def test_classify_manual_svtype_empty():
    vcf = viola.read_vcf(os.path.join(HERE, 'data/manta1.vcf'))
    vcf_empty = viola.read_vcf(os.path.join(HERE, 'data/manta1_empty.vcf'))
    vcf2 =  vcf.copy()
    vcf_empty2 = vcf_empty.copy()
    vcf = vcf.breakend2breakpoint()
    vcf2 = vcf2.breakend2breakpoint()
    multi_vcf = viola.MultiVcf([vcf, vcf_empty, vcf2, vcf_empty2], ['vcf1', 'empty1', 'vcf2', 'empty2'])
    ls_conditions = [small_del, large_del, small_dup, large_dup, small_inv, tra]
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    result = multi_vcf.classify_manual_svtype(ls_conditions=ls_conditions, ls_names=ls_names)

    manual_sv_type = multi_vcf.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.DataFrame([[2, 3, 1, 0, 2, 2, 1], [0, 0, 0, 0, 0, 0, 0], [2, 3, 1, 0, 2, 2, 1], [0, 0, 0, 0, 0, 0, 0]])
    result_expected.columns = ls_names + ['others']
    result_expected.columns.name = 'manual_sv_type'
    result_expected.index = ['vcf1', 'empty1', 'vcf2', 'empty2']
    result_expected.index.name = 'patients'
    pd.testing.assert_frame_equal(result, result_expected)
"""

def test_classify_manual_svtype_from_definitions():
    vcf = viola.read_vcf(os.path.join(HERE, 'data/manta1.vcf'))
    vcf2 =  vcf.copy()
    vcf = vcf.breakend2breakpoint()
    vcf2 = vcf2.breakend2breakpoint()
    multi_vcf = viola.MultiVcf([vcf, vcf2], ['vcf1', 'vcf2'])
    ls_conditions = [small_del, large_del, small_dup, large_dup, small_inv, tra]
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    result = multi_vcf.classify_manual_svtype(definitions=StringIO(DEFINITIONS))

    manual_sv_type = multi_vcf.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.DataFrame([[2, 3, 1, 0, 2, 2, 1], [2, 3, 1, 0, 2, 2, 1]])
    result_expected.columns = ls_names + ['others']
    result_expected.columns.name = 'manual_sv_type'
    result_expected.index = ['vcf1', 'vcf2']
    result_expected.index.name = 'patients'
    pd.testing.assert_frame_equal(result, result_expected)

def test_classify_manual_svtype_from_file():
    vcf = viola.read_vcf(os.path.join(HERE, 'data/manta1.vcf'))
    vcf2 =  vcf.copy()
    vcf = vcf.breakend2breakpoint()
    vcf2 = vcf2.breakend2breakpoint()
    multi_vcf = viola.MultiVcf([vcf, vcf2], ['vcf1', 'vcf2'])
    ls_conditions = [small_del, large_del, small_dup, large_dup, small_inv, tra]
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    path = os.path.join(HERE, '../vcf/classify_manual_svtype/data/example_definition.txt')
    result = multi_vcf.classify_manual_svtype(definitions=path)

    manual_sv_type = multi_vcf.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.DataFrame([[2, 3, 1, 0, 2, 2, 1], [2, 3, 1, 0, 2, 2, 1]])
    result_expected.columns = ls_names + ['others']
    result_expected.columns.name = 'manual_sv_type'
    result_expected.index = ['vcf1', 'vcf2']
    result_expected.index.name = 'patients'
    pd.testing.assert_frame_equal(result, result_expected)