import viola
from viola.testing import assert_vcf_equal
import os
HERE = os.path.abspath(os.path.dirname(__file__))


def test_add_info_table():
    vcf = viola.read_vcf(os.path.join(HERE, 'data/test.manta.vcf'))
    vcf_expected = viola.read_vcf(os.path.join(HERE, 'data/test.info.added.manta.vcf'))
    table = vcf.svlen
    table.columns = ['id', 'value_idx', 'test']
    vcf.add_info_table('test', table, 1, 'Integer', 'test info')
    assert_vcf_equal(vcf, vcf_expected)
