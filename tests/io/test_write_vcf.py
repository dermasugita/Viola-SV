import pytest
import sgt
import sys, os
import filecmp
import numpy as np
HERE = os.path.abspath(os.path.dirname(__file__))

class TestWriteVcf:
    # get vcf class
    manta_path = os.path.join(HERE, 'data/manta1.inv.vcf')
    df = sgt.read_vcf(manta_path)
    def test_write_vcf(self):
        # self.df.to_csv('path/to/csv')
        self.df.to_vcf('tests/io/output/write_vcf.vcf')
        assert filecmp.cmp('tests/io/output/write_vcf.vcf', 'tests/io/data/manta1.inv.vcf')

    def test_write_info(self):
        self.df.to_vcf('tests/io/output/write_info.vcf', onlyinfo=True)
        assert filecmp.cmp('tests/io/output/write_info.vcf', 'tests/io/data/manta1.inv_info.vcf')






