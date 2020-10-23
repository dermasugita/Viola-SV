import pytest
import sgt
import sys, os
import numpy as np
HERE = os.path.abspath(os.path.dirname(__file__))

class TestReadVcf:
    manta_path = os.path.join(HERE, 'data/manta1.inv.vcf')
    result = sgt.read_vcf(manta_path)
    def test_read_vcf_info_cipos(self):
        info_cipos = self.result.get_table('cipos')
        assert list(info_cipos.columns) == ['id', 'value_idx', 'cipos']
        cipos_value_idx = np.unique(info_cipos['value_idx'])
        cipos_expected = np.array([0, 1])
        np.testing.assert_array_equal(cipos_value_idx, cipos_expected)