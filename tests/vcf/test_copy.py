import viola
from viola.testing import assert_vcf_equal
import os
HERE = os.path.abspath(os.path.dirname(__file__))


def test_copy():
    manta_path = os.path.join(HERE, '../io/data/test.manta.vcf')
    delly_path = os.path.join(HERE, '../io/data/test.delly.vcf')
    lumpy_path = os.path.join(HERE, '../io/data/test.lumpy.vcf')
    gridss_path = os.path.join(HERE, '../io/data/test.gridss.vcf')
    manta_vcf = viola.read_vcf(
        manta_path,
        variant_caller='manta',
        patient_name='test')
    delly_vcf = viola.read_vcf(
        delly_path, 
        variant_caller='delly',
        patient_name='test')
    lumpy_vcf = viola.read_vcf(
        lumpy_path, 
        variant_caller='lumpy',
        patient_name='test')
    gridss_vcf = viola.read_vcf(
        gridss_path,
        variant_caller='gridss',
        patient_name='test')
    assert_vcf_equal(manta_vcf, manta_vcf.copy())
    assert_vcf_equal(delly_vcf, delly_vcf.copy())
    assert_vcf_equal(lumpy_vcf, lumpy_vcf.copy())
    assert_vcf_equal(gridss_vcf, gridss_vcf.copy())

def test_copy2():
    manta_path = os.path.join(HERE, '../io/data/test.manta.vcf')
    delly_path = os.path.join(HERE, '../io/data/test.delly.vcf')
    lumpy_path = os.path.join(HERE, '../io/data/test.lumpy.vcf')
    gridss_path = os.path.join(HERE, '../io/data/test.gridss.vcf')
    manta_vcf = viola.read_vcf2(
        manta_path,
        variant_caller='manta',
        patient_name='test')
    delly_vcf = viola.read_vcf2(
        delly_path, 
        variant_caller='delly',
        patient_name='test')
    lumpy_vcf = viola.read_vcf2(
        lumpy_path, 
        variant_caller='lumpy',
        patient_name='test')
    gridss_vcf = viola.read_vcf2(
        gridss_path,
        variant_caller='gridss',
        patient_name='test')
    assert_vcf_equal(manta_vcf, manta_vcf.copy())
    assert_vcf_equal(delly_vcf, delly_vcf.copy())
    assert_vcf_equal(lumpy_vcf, lumpy_vcf.copy())
    assert_vcf_equal(gridss_vcf, gridss_vcf.copy())
