import viola    
import numpy as np
import os
HERE = os.path.abspath(os.path.dirname(__file__))
manta = viola.read_vcf(os.path.join(HERE, 'data/test.merge.manta.vcf'), variant_caller='manta')
delly = viola.read_vcf(os.path.join(HERE, 'data/test.merge.delly.vcf'), variant_caller='delly')
lumpy = viola.read_vcf(os.path.join(HERE, 'data/test.merge.lumpy.vcf'), variant_caller='lumpy')
gridss = viola.read_vcf(os.path.join(HERE, 'data/test.merge.gridss.vcf'), variant_caller='gridss')
def test_merge():
    merged_all = gridss.merge(ls_vcf=[delly, lumpy, manta, gridss], threshold=100)
    viola.testing.assert_vcf_equal(merged_all, merged_all.copy())