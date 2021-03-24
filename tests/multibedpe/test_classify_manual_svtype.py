import os
import viola
import pandas as pd
from io import StringIO
HERE = os.path.abspath(os.path.dirname(__file__))
data = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2
chr1	10	11	chr1	20	21	test1	60	+	-
chr1	10	11	chr1	25	26	test2	60	+	-
chr1	100	101	chr1	250	251	test3	60	+	-
chr1	105	106	chr1	290	291	test4	60	+	-
chr1	150	151	chr1	300	301	test5	60	+	-
chr2	10	11	chr2	20	21	test6	60	-	+
chr2	10	11	chr2	20	21	test7	60	-	-
chr2	100	101	chr2	280	281	test8	60	-	-
chr2	180	181	chr2	2000	2001	test9	60	-	-
chr2	10	11	chr2	40	41	test10	60	+	+
chr2	10	11	chr5	20	21	test11	60	+	-
chr3	10	11	chr4	20	21	test12	60	-	-
"""
data_expected = """bedpe1_test1	0	small_del
bedpe2_test1	0	small_del
bedpe1_test2	0	small_del
bedpe2_test2	0	small_del
bedpe1_test3	0	large_del
bedpe2_test3	0	large_del
bedpe1_test4	0	large_del
bedpe2_test4	0	large_del
bedpe1_test5	0	large_del
bedpe2_test5	0	large_del
bedpe1_test6	0	small_dup
bedpe2_test6	0	small_dup
bedpe1_test7	0	small_inv
bedpe2_test7	0	small_inv
bedpe1_test8	0	others
bedpe2_test8	0	others
bedpe1_test9	0	others
bedpe2_test9	0	others
bedpe1_test10	0	small_inv
bedpe2_test10	0	small_inv
bedpe1_test11	0	tra
bedpe2_test11	0	tra
bedpe1_test12	0	tra
bedpe2_test12	0	tra
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
0 SVTYPE == BND
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
    return x.filter('svtype == BND').ids

def test_classify_manual_svtype():
    bedpe1 = viola.read_bedpe(StringIO(data))
    bedpe2 = viola.read_bedpe(StringIO(data))
    ls_conditions = [small_del, large_del, small_dup, large_dup, small_inv, tra]
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    multibedpe = viola.MultiBedpe([bedpe1, bedpe2], ['bedpe1', 'bedpe2'])
    result = multibedpe.classify_manual_svtype(ls_conditions=ls_conditions, ls_names=ls_names)
    manual_sv_type = multibedpe.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.DataFrame([[2, 3, 1, 0, 2, 2, 2],[2, 3, 1, 0, 2, 2, 2]])
    result_expected.columns = ls_names + ['others']
    result_expected.columns.name = 'manual_sv_type'
    result_expected.index = ['bedpe1', 'bedpe2']
    result_expected.index.name = 'patients'
    pd.testing.assert_frame_equal(result, result_expected)
    
def test_classify_manual_svtype_from_definitions():
    bedpe1 = viola.read_bedpe(StringIO(data))
    bedpe2 = viola.read_bedpe(StringIO(data))
    ls_conditions = [small_del, large_del, small_dup, large_dup, small_inv, tra]
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    multibedpe = viola.MultiBedpe([bedpe1, bedpe2], ['bedpe1', 'bedpe2'])
    result = multibedpe.classify_manual_svtype(definitions=StringIO(DEFINITIONS))
    manual_sv_type = multibedpe.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.DataFrame([[2, 3, 1, 0, 2, 2, 2],[2, 3, 1, 0, 2, 2, 2]])
    result_expected.columns = ls_names + ['others']
    result_expected.columns.name = 'manual_sv_type'
    result_expected.index = ['bedpe1', 'bedpe2']
    result_expected.index.name = 'patients'
    pd.testing.assert_frame_equal(result, result_expected)

def test_classify_manual_svtype_from_file():
    bedpe1 = viola.read_bedpe(StringIO(data))
    bedpe2 = viola.read_bedpe(StringIO(data))
    ls_conditions = [small_del, large_del, small_dup, large_dup, small_inv, tra]
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    multibedpe = viola.MultiBedpe([bedpe1, bedpe2], ['bedpe1', 'bedpe2'])
    path = os.path.join(HERE, '../bedpe/data/example_definition.txt')
    result = multibedpe.classify_manual_svtype(definitions=path)
    manual_sv_type = multibedpe.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.DataFrame([[2, 3, 1, 0, 2, 2, 2],[2, 3, 1, 0, 2, 2, 2]])
    result_expected.columns = ls_names + ['others']
    result_expected.columns.name = 'manual_sv_type'
    result_expected.index = ['bedpe1', 'bedpe2']
    result_expected.index.name = 'patients'
    pd.testing.assert_frame_equal(result, result_expected)