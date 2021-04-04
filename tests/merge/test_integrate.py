import viola    
import numpy as np
import os
HERE = os.path.abspath(os.path.dirname(__file__))
delly = viola.read_vcf(os.path.join(HERE, 'data', "test.merge.delly.vcf"), variant_caller="delly")
manta = viola.read_vcf(os.path.join(HERE, 'data', "test.merge.manta.vcf"), variant_caller="manta")
gridss = viola.read_vcf(os.path.join(HERE, 'data', "test.merge.gridss.vcf"), variant_caller="gridss")
lumpy = viola.read_vcf(os.path.join(HERE, 'data', "test.merge.lumpy.vcf"), variant_caller="lumpy")

def test_merge_to_vcf_like():
    merged = viola.merge([manta, gridss, delly, lumpy], integration=True)
    print(merged)
    merged_integrated = merged.filter('supportingcallercount > 1')
    print(merged_integrated)
