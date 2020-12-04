import sgt
def test_read_vcf():
    path = './resources/vcf/lumpy1_sample.vcf'
    sgt_ob = sgt.read_vcf(path)
    print(sgt_ob)
    print(sgt_ob.get_table("positions"))