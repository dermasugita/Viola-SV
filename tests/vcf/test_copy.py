import viola
from viola.testing import assert_vcf_equal
import os
HERE = os.path.abspath(os.path.dirname(__file__))


def test_copy():
    manta_path = os.path.join(HERE, '../io/data/test.manta.vcf')
    vcf = viola.read_vcf(manta_path)
    vcf_copy = vcf.copy()
    assert_vcf_equal(vcf, vcf_copy)