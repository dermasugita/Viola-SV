import pandas
import vcf

def read_vcf(
    filepath,
    variant_caller="manta"
):
    vcf_reader = vcf.Reader(open(filepath, 'r'))