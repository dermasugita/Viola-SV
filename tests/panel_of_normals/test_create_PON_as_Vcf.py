import viola    
import os
HERE = os.path.abspath(os.path.dirname(__file__))
normal1 = viola.read_vcf(os.path.join(HERE, 'data/normal1.vcf'), variant_caller='manta', patient_name='sample1')
normal2 = viola.read_vcf(os.path.join(HERE, 'data/normal2.vcf'), variant_caller='manta', patient_name='sample2')
normal3 = viola.read_vcf(os.path.join(HERE, 'data/normal3.vcf'), variant_caller='manta', patient_name='sample3')
normal4 = viola.read_vcf(os.path.join(HERE, 'data/normal4.vcf'), variant_caller='manta', patient_name='sample4')


def test_create_PON_as_VCF():
    viola.core.cohort.create_PON_as_Vcf([normal1, normal2, normal3, normal4])

def test_create_PON_as_VCF_with_priority():
    viola.core.cohort.create_PON_as_Vcf([normal1, normal2, normal3, normal4], priority=['sample4', 'sample2', 'sample3', 'sample1'])