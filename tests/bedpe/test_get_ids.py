import viola
from io import StringIO
data = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2
chr1	10	11	chr1	20	21	test1	60	+	-
chr1	10	11	chr1	25	26	test2	60	+	-
chr1	100	101	chr1	250	251	test3	60	+	-
chr1	105	106	chr1	290	291	test4	60	+	-
chr1	150	151	chr1	300	301	test5	60	+	-
"""
buf = StringIO(data)

def test_ged_ids():
    bedpe = viola.read_bedpe(buf, patient_name='test')
    assert bedpe.get_ids() == {'test1', 'test2', 'test3', 'test4', 'test5'}