import viola
import sys, os
import pandas as pd
HERE = os.path.abspath(os.path.dirname(__file__))

def test_merge():
    manta = viola.read_vcf(os.path.join(HERE, 'data/test.merge.manta.vcf'), variant_caller='manta')
    delly = viola.read_vcf(os.path.join(HERE, 'data/test.merge.delly.vcf'), variant_caller='delly')
    lumpy = viola.read_vcf(os.path.join(HERE, 'data/test.merge.lumpy.vcf'), variant_caller='lumpy')
    gridss = viola.read_vcf(os.path.join(HERE, 'data/test.merge.gridss.vcf'), variant_caller='gridss')

    merged = manta.merge(threshold=100, ls_caller_names=['manta', 'delly', 'lumpy', 'gridss'], ls_vcf=[manta, delly, lumpy, gridss])

