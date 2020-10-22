import pytest
import sgt
import sys, os
HERE = os.path.abspath(os.path.dirname(__file__))

class TestToBedpe:
    manta_path = os.path.join(HERE, 'data/manta1.inv.vcf')
    result = sgt.read_vcf(manta_path)

    def test_to_bedpe(self):
        bedpe = self.result.to_bedpe_like()
        assert list(bedpe.columns) == ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
                                        'name', 'score', 'strand1', 'strand2']