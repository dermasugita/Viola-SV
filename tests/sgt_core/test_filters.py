import sgt
import pandas as pd
import sys, os

HERE = os.path.abspath(os.path.dirname(__file__))

class TestFilters:
    path = os.path.join(HERE, 'data/manta1.inv.vcf')
    obj = sgt.read_vcf(path)
    
    def test_filter_by_info(self):
        pass