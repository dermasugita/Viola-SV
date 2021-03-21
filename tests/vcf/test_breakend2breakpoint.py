import viola
from viola.testing import assert_vcf_equal
import os
HERE = os.path.abspath(os.path.dirname(__file__))


def test_breakend2breakpoint():
    manta_path = os.path.join(HERE, '../io/data/test.manta.vcf')
    vcf = viola.read_vcf(manta_path)
    vcf_result = vcf.breakend2breakpoint()