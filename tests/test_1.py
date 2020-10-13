import pytest
import sgt
def test_read_vcf():
    path = './resoreces/vcf/lumpy1_sample.vcf'
    sgt_ob = sgt.read_vcf(path)
    sgt_bedpe = sgt_ob.to_bedpe_like()
    print(sgt_bedpe.head(20))