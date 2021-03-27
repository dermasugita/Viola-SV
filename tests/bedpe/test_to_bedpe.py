import viola
import pandas as pd
from io import StringIO
import os
import filecmp
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

def test_to_bedpe():
    bedpe = viola.read_bedpe(StringIO(data))
    bedpe.to_bedpe(os.path.join(HERE, 'data/output.bedpe'))
    bedpe.to_bedpe(os.path.join(HERE, 'data/output.svlen.bedpe'), custom_infonames=['svlen'])
    assert filecmp.cmp(os.path.join(HERE, 'data/output.bedpe'), os.path.join(HERE, 'data/bedpe_expected.bedpe'))
    assert filecmp.cmp(os.path.join(HERE, 'data/output.svlen.bedpe'), os.path.join(HERE, 'data/output.svlen.expected.bedpe'))