import viola
import pandas as pd
from io import StringIO
import pytest
from viola._exceptions import (
    InfoNotFoundError,
    DestructiveTableValueError,
    TableValueConfliction,
)

DATA = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test1	test2
chr1	10	11	chr1	20	21	test1	60	+	-	True	1
chr2	10	11	chr2	20	21	test2	60	-	+	False	23
chr2	10	11	chr2	40	41	test3	60	+	+	True	32
chr3	10	11	chr4	20	21	test4	60	-	-	True	13
"""
DATA_expected = """chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	test3	test2
chr1	10	11	chr1	20	21	test1	60	+	-	True	1
chr2	10	11	chr2	20	21	test2	60	-	+	False	23
chr2	10	11	chr2	40	41	test3	60	+	+	True	32
chr3	10	11	chr4	20	21	test4	60	-	-	True	13
"""


def test_rename_info():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    bedpe_expected = viola.read_bedpe(StringIO(DATA_expected), patient_name="patient1")
    bedpe.rename_info('test1', 'test3')
    viola.testing.assert_bedpe_equal(bedpe, bedpe_expected)

def test_rename_info_info_not_found():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    with pytest.raises(InfoNotFoundError):
        bedpe.rename_info('test3', 'test4')
    
def test_rename_info_destructive_value():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    with pytest.raises(DestructiveTableValueError):
        bedpe.rename_info('svlen', 'test4')
    with pytest.raises(DestructiveTableValueError):
        bedpe.rename_info('svtype', 'test4')
    with pytest.raises(DestructiveTableValueError):
        bedpe.rename_info('cipos', 'test4')
    with pytest.raises(DestructiveTableValueError):
        bedpe.rename_info('ciend', 'test4')

def test_rename_info_escaping_destructive_value():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    bedpe.rename_info('svlen', 'test4', False)

def test_rename_info_name_conflict():
    bedpe = viola.read_bedpe(StringIO(DATA), patient_name="patient1")
    with pytest.raises(TableValueConfliction):
        bedpe.rename_info('test2', 'test1')