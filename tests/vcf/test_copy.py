import viola
from viola.testing import assert_vcf_equal
import os
HERE = os.path.abspath(os.path.dirname(__file__))


def test_copy():
    manta_path = os.path.join(HERE, '../io/data/test.manta.vcf')
    delly_path = os.path.join(HERE, '../io/data/test.delly.vcf')
    lumpy_path = os.path.join(HERE, '../io/data/test.lumpy.vcf')
    gridss_path = os.path.join(HERE, '../io/data/test.gridss.vcf')
    manta_vcf = viola.read_vcf(manta_path)
    delly_vcf = viola.read_vcf(delly_path, variant_caller='delly')
    lumpy_vcf = viola.read_vcf(lumpy_path, variant_caller='lumpy')
    gridss_vcf = viola.read_vcf(gridss_path, variant_caller='gridss')
    assert_vcf_equal(manta_vcf, manta_vcf.copy())
    assert_vcf_equal(delly_vcf, delly_vcf.copy())
    assert_vcf_equal(lumpy_vcf, lumpy_vcf.copy())
    assert_vcf_equal(gridss_vcf, gridss_vcf.copy())
