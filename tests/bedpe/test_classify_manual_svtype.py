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
data_expected = """test1	0	small_del
test2	0	small_del
test3	0	large_del
test4	0	large_del
test5	0	large_del
test6	0	small_dup
test7	0	small_inv
test8	0	others
test9	0	others
test10	0	small_inv
test11	0	tra
test12	0	tra
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
    bedpe = viola.read_bedpe(StringIO(data), patient_name="patient1")
    ls_conditions = [small_del, large_del, small_dup, large_dup, small_inv, tra]
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    result = bedpe.classify_manual_svtype(ls_conditions=ls_conditions, ls_names=ls_names)
    manual_sv_type = bedpe.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.Series([2, 3, 1, 0, 2, 2, 2])
    result_expected.index = ls_names + ['others']
    result_expected.name = 'manual_sv_type'
    pd.testing.assert_series_equal(result, result_expected)
    
def test_classify_manual_svtype_from_definitions():
    bedpe = viola.read_bedpe(StringIO(data), patient_name="patient1")
    result = bedpe.classify_manual_svtype(definitions=StringIO(DEFINITIONS))
    manual_sv_type = bedpe.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.Series([2, 3, 1, 0, 2, 2, 2])
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    result_expected.index = ls_names + ['others']
    result_expected.name = 'manual_sv_type'
    pd.testing.assert_series_equal(result, result_expected)

def test_classify_manual_svtype_from_default():
    bedpe = viola.read_bedpe(StringIO(data), patient_name="patient1")
    try:
        result = bedpe.classify_manual_svtype(definitions="default")
    except TypeError:
        pass

def test_classify_manual_svtype_from_article():
    bedpe = viola.read_bedpe(StringIO(data), patient_name="patient1")
    try:
        result = bedpe.classify_manual_svtype(definitions="article")
    except TypeError:
        pass

def test_classify_manual_svtype_from_file():
    bedpe = viola.read_bedpe(StringIO(data), patient_name="patient1")
    path = os.path.join(HERE, 'data/example_definition.txt')
    result = bedpe.classify_manual_svtype(definitions=path)
    manual_sv_type = bedpe.manual_sv_type
    manual_sv_type.set_index('id', inplace=True)
    manual_sv_type_expected = pd.read_csv(StringIO(data_expected), sep='\t', names=('id', 'value_idx', 'manual_sv_type'))
    manual_sv_type_expected.set_index('id', inplace=True)
    pd.testing.assert_frame_equal(manual_sv_type, manual_sv_type_expected, check_like=True)

    result_expected = pd.Series([2, 3, 1, 0, 2, 2, 2])
    ls_names = ['small_del', 'large_del', 'small_dup', 'large_dup', 'small_inv', 'tra']
    result_expected.index = ls_names + ['others']
    result_expected.name = 'manual_sv_type'
    pd.testing.assert_series_equal(result, result_expected)
