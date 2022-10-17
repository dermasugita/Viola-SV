import viola
import pandas as pd
from io import StringIO
import pytest
from viola._exceptions import (
    InfoNotFoundError,
    DestructiveTableValueError,
    TableValueConfliction,
)

DATA_manta = """##fileformat=VCFv4.1
##contig=<ID=chr1,length=195471971>
##INFO=<ID=TEST,Number=1,Type=String,Description="Test INFO">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##ALT=<ID=DEL,Description="Deletion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mouse1_N	mouse1_T
chr1	82550461	test1	G	<DEL>	.	.	TEST=test;SVTYPE=DEL;SVLEN=-3764;END=82554225;CIPOS=-51,52;CIEND=-51,52	PR	21,0	43,4
"""
DATA_manta_expected = """##fileformat=VCFv4.1
##contig=<ID=chr1,length=195471971>
##INFO=<ID=RENAMED,Number=1,Type=String,Description="Test INFO">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##ALT=<ID=DEL,Description="Deletion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mouse1_N	mouse1_T
chr1	82550461	test1	G	<DEL>	.	.	RENAMED=test;SVTYPE=DEL;SVLEN=-3764;END=82554225;CIPOS=-51,52;CIEND=-51,52	PR	21,0	43,4
"""

def test_rename_info():
    vcf = viola.read_vcf(StringIO(DATA_manta), patient_name="patient1")
    vcf_expected = viola.read_vcf(StringIO(DATA_manta_expected), patient_name="patient1")
    vcf.rename_info('test', 'renamed')
    viola.testing.assert_vcf_equal(vcf, vcf_expected)

def test_rename_info_protected():
    vcf = viola.read_vcf(StringIO(DATA_manta), patient_name="patient1")
    with pytest.raises(DestructiveTableValueError):
        vcf.rename_info('svtype', 'renamed')
    with pytest.raises(DestructiveTableValueError):
        vcf.rename_info('svlen', 'renamed')
    with pytest.raises(DestructiveTableValueError):
        vcf.rename_info('cipos', 'renamed')
    with pytest.raises(DestructiveTableValueError):
        vcf.rename_info('ciend', 'renamed')

def test_rename_info_not_found():
    vcf = viola.read_vcf(StringIO(DATA_manta), patient_name="patient1")
    with pytest.raises(InfoNotFoundError):
        vcf.rename_info('svtypes', 'renamed')

def test_rename_info_protected_excape():
    vcf = viola.read_vcf(StringIO(DATA_manta), patient_name="patient1")
    vcf.rename_info('svtype', 'renamed', False)

def test_rename_info_name_conflict():
    vcf = viola.read_vcf(StringIO(DATA_manta), patient_name="patient1")
    with pytest.raises(TableValueConfliction):
        vcf.rename_info('test', 'cipos')
    with pytest.raises(TableValueConfliction):
        vcf.rename_info('test', 'positions')