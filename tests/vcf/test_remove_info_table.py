import viola
from viola.testing import assert_vcf_equal
import os
HERE = os.path.abspath(os.path.dirname(__file__))


def test_remove_info_table():
    vcf = viola.read_vcf(
        os.path.join(HERE, 'data/test.info.added.manta.vcf'),
        variant_caller='manta',
        patient_name='test')
    vcf_expected = viola.read_vcf(
        os.path.join(HERE, 'data/test.manta.vcf'),
        variant_caller='manta',
        patient_name='test')
    vcf.remove_info_table('test')
    assert_vcf_equal(vcf, vcf_expected)

def test_remove_info_table2():
    vcf = viola.read_vcf2(
        os.path.join(HERE, 'data/test.info.added.manta.vcf'),
        variant_caller='manta',
        patient_name='test')
    vcf_expected = viola.read_vcf2(
        os.path.join(HERE, 'data/test.manta.vcf'),
        variant_caller='manta',
        patient_name='test')
    vcf.remove_info_table('test')
    assert_vcf_equal(vcf, vcf_expected)
