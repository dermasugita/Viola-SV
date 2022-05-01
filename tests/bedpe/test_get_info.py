import viola
import pandas as pd
import pytest
from io import StringIO
from viola._exceptions import InfoNotFoundError
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


def test_get_info():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    print(bedpe)
    info = bedpe.get_info('svtype')
    pd.testing.assert_frame_equal(info, bedpe.get_table('svtype'))


def test_get_info_not_found():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    with pytest.raises(InfoNotFoundError):
        bedpe.get_info('qwerty')
