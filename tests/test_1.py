import viola
def test_read_vcf():
    path = './resources/vcf/lumpy1_sample.vcf'
    viola_ob = viola.read_vcf(path)
    print(viola_ob)
    print(viola_ob.get_table("positions"))