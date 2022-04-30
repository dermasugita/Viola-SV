import viola
import pandas as pd
from io import StringIO
DATA = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test1
chr1	10	11	chr1	20	21	test1	60	+	-	True
chr1	10	11	chr1	25	26	test2	60	+	-	False
chr1	100	101	chr1	250	251	test3	60	+	-	True
chr1	105	106	chr1	290	291	test4	60	+	-	False
chr1	150	151	chr1	300	301	test5	60	+	-	True
chr2	10	11	chr2	20	21	test6	60	-	+	False
chr2	10	11	chr2	20	21	test7	60	-	-	True
chr2	100	101	chr2	280	281	test8	60	-	-	False
chr2	180	181	chr2	2000	2001	test9	60	-	-	False
chr2	10	11	chr2	40	41	test10	60	+	+	True
chr2	10	11	chr5	20	21	test11	60	+	-	False
chr3	10	11	chr4	20	21	test12	60	-	-	True
"""
data_expected1 = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test1
chr1	10	11	chr1	20	21	test1	60	+	-	True
chr1	10	11	chr1	25	26	test2	60	+	-	False
chr2	10	11	chr2	20	21	test6	60	-	+	False
"""
data_expected2 = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test1
chr1	10	11	chr1	20	21	test1	60	+	-	True
chr1	10	11	chr1	25	26	test2	60	+	-	False
"""
data_expected3 = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test1
chr2	10	11	chr2	20	21	test6	60	-	+	False
chr2	10	11	chr2	20	21	test7	60	-	-	True
chr2	100	101	chr2	280	281	test8	60	-	-	False
chr2	180	181	chr2	2000	2001	test9	60	-	-	False
chr2	10	11	chr2	40	41	test10	60	+	+	True
"""
data_expected4 = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test1
chr2	10	11	chr2	20	21	test6	60	-	+	False
chr2	10	11	chr2	20	21	test7	60	-	-	True
chr2	100	101	chr2	280	281	test8	60	-	-	False
chr2	180	181	chr2	2000	2001	test9	60	-	-	False
chr2	10	11	chr2	40	41	test10	60	+	+	True
chr2	10	11	chr5	20	21	test11	60	+	-	False
"""
data_expected5 = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test1
chr1	10	11	chr1	20	21	test1	60	+	-	True
chr1	10	11	chr1	25	26	test2	60	+	-	False
chr1	100	101	chr1	250	251	test3	60	+	-	True
chr1	105	106	chr1	290	291	test4	60	+	-	False
chr1	150	151	chr1	300	301	test5	60	+	-	True
chr2	10	11	chr5	20	21	test11	60	+	-	False
chr3	10	11	chr4	20	21	test12	60	-	-	True
"""

def test_filter():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    bedpe_expected = viola.read_bedpe(StringIO(data_expected1), patient_name='patient1')
    f_bedpe = bedpe.filter(['svtype == DEL', 'svtype == DUP', 'svlen > -100'], query_logic='(0 | 1) & 2')
    viola.testing.assert_bedpe_equal(f_bedpe, bedpe_expected)

def test_filter_and():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    bedpe_expected = viola.read_bedpe(StringIO(data_expected2), patient_name='patient1')
    f_bedpe = bedpe.filter(['svtype == DEL', 'svlen > -100'], query_logic='and')
    viola.testing.assert_bedpe_equal(f_bedpe, bedpe_expected)

def test_filter_or():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    bedpe_expected = viola.read_bedpe(StringIO(data_expected3), patient_name='patient1')
    f_bedpe = bedpe.filter(['svtype == DUP', 'svtype == INV'], query_logic='or')
    viola.testing.assert_bedpe_equal(f_bedpe, bedpe_expected)

def test_filter_pos():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    bedpe_expected = viola.read_bedpe(StringIO(data_expected4), patient_name='patient1')
    f_bedpe = bedpe.filter(['be1 chr2'])
    viola.testing.assert_bedpe_equal(f_bedpe, bedpe_expected)

def test_filter_pos_exclude():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    bedpe_expected = viola.read_bedpe(StringIO(data_expected5), patient_name='patient1')
    f_bedpe = bedpe.filter(['be2 !chr2'])
    viola.testing.assert_bedpe_equal(f_bedpe, bedpe_expected)