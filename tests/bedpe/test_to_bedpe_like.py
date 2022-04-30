
from distutils.sysconfig import customize_compiler
import viola
import pandas as pd
from io import StringIO
data = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2
chr1	10	11	chr1	20	21	test1	60	+	-
chr1	10	11	chr1	25	26	test2	60	+	-
chr1	100	101	chr1	250	251	test3	60	-	+
chr1	105	106	chr1	290	291	test4	60	+	+
chr1	150	151	chr2	300	301	test5	60	+	-
"""
data2 = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	svtype_0
chr1	10	11	chr1	20	21	test1	60	+	-	DEL
chr1	10	11	chr1	25	26	test2	60	+	-	DEL
chr1	100	101	chr1	250	251	test3	60	-	+	DUP
chr1	105	106	chr1	290	291	test4	60	+	+	INV
chr1	150	151	chr2	300	301	test5	60	+	-	BND
"""
data3 = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2
chr1	10	12	chr1	20	23	test1	60	+	-
chr1	10	14	chr1	25	22	test2	60	+	-
chr1	100	103	chr1	250	255	test3	60	-	+
chr1	105	103	chr1	290	293	test4	60	+	+
chr1	150	156	chr2	300	302	test5	60	+	-
"""

def test_to_bedpe_like():
    bedpe = viola.read_bedpe(StringIO(data), patient_name='test')
    df_bedpe = bedpe.to_bedpe_like()
    df_expected = pd.read_table(StringIO(data), index_col=None)
    pd.testing.assert_frame_equal(df_bedpe, df_expected)

def test_to_bedpe_like_custom_infonames():
    bedpe = viola.read_bedpe(StringIO(data), patient_name='test')
    df_bedpe = bedpe.to_bedpe_like(custom_infonames=['svtype'])
    df_expected = pd.read_table(StringIO(data2), index_col=None)
    pd.testing.assert_frame_equal(df_bedpe, df_expected)

def test_to_bedpe_like_confidence_intervals():
    bedpe = viola.read_bedpe(StringIO(data3), patient_name='test')
    df_bedpe = bedpe.to_bedpe_like(confidence_intervals=True)
    df_expected = pd.read_table(StringIO(data3), index_col=None)
    pd.testing.assert_frame_equal(df_bedpe, df_expected)