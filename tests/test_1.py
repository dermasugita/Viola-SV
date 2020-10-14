import pytest
import sgt
def test_read_vcf():
    path = './resources/vcf/lumpy1_sample.vcf'
    sgt_ob = sgt.read_vcf(path)
    #sgt_bedpe = sgt_ob.to_bedpe_like()
    print(sgt_ob.get_table('positions'))
    print(sgt_ob.get_table('positions').columns)
    