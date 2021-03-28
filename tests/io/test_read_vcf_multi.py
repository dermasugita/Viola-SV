import viola
import pandas as pd
from io import StringIO
import os
HERE = os.path.abspath(os.path.dirname(__file__))

def test_read_vcf_multi_with_empty():
    multi_vcf = viola.read_vcf_multi(os.path.join(HERE, 'data/multivcf'), variant_caller='manta')
    assert set(multi_vcf._ls_patients) == set(['empty.manta1', 'empty.manta2', 'test.manta1', 'test.manta2'])

def test_read_vcf_multi_without_empty():
    multi_vcf = viola.read_vcf_multi(os.path.join(HERE, 'data/multivcf'), variant_caller='manta', exclude_empty_cases=True)
    assert set(multi_vcf._ls_patients) == set(['test.manta1', 'test.manta2'])