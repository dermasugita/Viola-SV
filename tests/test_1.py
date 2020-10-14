import pytest
import sgt
def test_read_vcf():
    path = './resources/vcf/lumpy1_sample.vcf'
    path2 = './resources/bedpe/manta1.bedpe'
    sgt_ob = sgt.read_vcf(path)
    sgt_ob.get_unique_events()
    #sgt_bedpe = sgt_ob.to_bedpe_like(add_filters=True, unique_events=True)
    #print(sgt_bedpe)
    #print(sgt_bedpe.columns)
    