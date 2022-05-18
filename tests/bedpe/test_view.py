import viola
import pandas as pd
import pytest
from io import StringIO
from viola._exceptions import InfoNotFoundError
DATA = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test1
chr1	10	11	chr1	20	21	test1	60	+	-	True
chr2	10	11	chr2	20	21	test2	60	-	+	False
chr2	10	11	chr2	40	41	test3	60	+	+	True
chr3	10	11	chr4	20	21	test4	60	-	-	True
"""

def test_view():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    view = bedpe.view()
    view_expected = """INFO=svlen,svtype,cipos,ciend,test1
Documentation of Bedpe object ==> https://dermasugita.github.io/ViolaDocs/docs/html/reference/bedpe.html
"""
    txt = """id	be1	be2	strand	qual	svtype
test1	chr1:11	chr1:21	+-	60	DEL
test2	chr2:11	chr2:21	-+	60	DUP
test3	chr2:11	chr2:41	++	60	INV
test4	chr3:11	chr4:21	--	60	BND
"""
    df_expected = pd.read_table(StringIO(txt))
    view_expected += str(df_expected)
    assert view == view_expected