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
        self.df.to_vcf('output/write_vcf.vcf')
        assert filecmp('output/write_vcf.vcf', 'data/manta1.inv.vcf')






