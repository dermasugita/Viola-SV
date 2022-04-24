import pytest
import viola
import sys, os
import filecmp
import numpy as np
HERE = os.path.abspath(os.path.dirname(__file__))

class TestWriteVcf:
    # get vcf class
    #manta_path = os.path.join(HERE, 'data/manta1.inv.vcf')
    manta_path = os.path.join(HERE, 'data/test.manta.vcf')
    lumpy_path = os.path.join(HERE, 'data/test.lumpy.vcf')
    delly_path = os.path.join(HERE, 'data/test.delly.vcf')
    gridss_path = os.path.join(HERE, 'data/test.gridss.vcf')
    vcf_manta = viola.read_vcf(manta_path, patient_name='manta1')
    vcf_lumpy = viola.read_vcf(lumpy_path, variant_caller='lumpy', patient_name='lumpy1')
    vcf_delly = viola.read_vcf(delly_path, variant_caller='delly', patient_name='delly1')
    vcf_gridss = viola.read_vcf(gridss_path, variant_caller='gridss', patient_name='gridss1')
    def test_write_vcf_manta(self):
        # self.df.to_csv('path/to/csv')
        self.vcf_manta.to_vcf('tests/io/output/write_vcf_manta.vcf')
        #assert filecmp.cmp('tests/io/output/write_vcf.vcf', 'tests/io/data/manta1.inv.vcf')
        assert filecmp.cmp('tests/io/output/write_vcf_manta.vcf', 'tests/io/data/test.manta.validation.vcf')

    def test_write_vcf_lumpy(self):
        # self.df.to_csv('path/to/csv')
        self.vcf_lumpy.to_vcf('tests/io/output/write_vcf_lumpy.vcf')
        #assert filecmp.cmp('tests/io/output/write_vcf.vcf', 'tests/io/data/manta1.inv.vcf')
        assert filecmp.cmp('tests/io/output/write_vcf_lumpy.vcf', 'tests/io/data/test.lumpy.validation.vcf')

    def test_write_vcf_delly(self):
        self.vcf_delly.to_vcf('tests/io/output/write_vcf_delly.vcf')
        assert filecmp.cmp('tests/io/output/write_vcf_delly.vcf', 'tests/io/data/test.delly.validation.vcf')

    def test_write_vcf_gridss(self):
        self.vcf_gridss.to_vcf('tests/io/output/write_vcf_gridss.vcf')
        assert filecmp.cmp('tests/io/output/write_vcf_gridss.vcf', 'tests/io/data/test.gridss.validation.vcf')

    def test_write_info(self):
        self.vcf_manta.to_vcf('tests/io/output/write_info.vcf', onlyinfo=True)
        #assert filecmp.cmp('tests/io/output/write_info.vcf', 'tests/io/data/manta1.inv_info.vcf')
        assert filecmp.cmp('tests/io/output/write_info.vcf', 'tests/io/data/test_info.vcf')
