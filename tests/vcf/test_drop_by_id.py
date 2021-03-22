import viola
from viola.testing import assert_vcf_equal
import os
HERE = os.path.abspath(os.path.dirname(__file__))

def test_drop_by_id():
    vcf = viola.read_vcf(os.path.join(HERE, 'data/test.manta.vcf'))
    vcf_dropped_expected = viola.read_vcf(os.path.join(HERE, 'data/test.dropped.manta.vcf'))
    vcf_dropped2_expected = viola.read_vcf(os.path.join(HERE, 'data/test.dropped2.manta.vcf'))
    vcf_dropped = vcf.drop_by_id('test1')
    vcf_dropped2 = vcf.drop_by_id(['test2', 'test3'])
    assert_vcf_equal(vcf_dropped, vcf_dropped_expected)
    assert_vcf_equal(vcf_dropped2, vcf_dropped2_expected)
