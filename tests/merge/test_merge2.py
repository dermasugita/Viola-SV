import viola    
import numpy as np
"""
delly = viola.read_vcf("/shared_data/share/merge/test.merge.delly.vcf", variant_caller="delly")
gridss = viola.read_vcf("/shared_data/share/merge/test.merge.gridss.vcf", variant_caller="gridss")
lumpy = viola.read_vcf("/shared_data/share/merge/test.merge.lumpy.vcf", variant_caller="lumpy")
manta = viola.read_vcf("/shared_data/share/merge/test.merge.manta.vcf", variant_caller="manta")
def test_merge():
    merged_all = gridss.merge(ls_vcf=[delly, lumpy, manta, gridss], threshold=100)
    viola.testing.assert_vcf_equal(merged_all, merged_all.copy())
"""