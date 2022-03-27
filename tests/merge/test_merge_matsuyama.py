import viola    
import numpy as np
import os
HERE = os.path.abspath(os.path.dirname(__file__))
manta = viola.read_vcf(os.path.join(HERE, 'data/test.merge.manta.vcf'), variant_caller='manta')
delly = viola.read_vcf(os.path.join(HERE, 'data/test.merge.delly.vcf'), variant_caller='delly')
lumpy = viola.read_vcf(os.path.join(HERE, 'data/test.merge.lumpy.vcf'), variant_caller='lumpy')
gridss = viola.read_vcf(os.path.join(HERE, 'data/test.merge.gridss.vcf'), variant_caller='gridss')

def test_merge():
    merged_all = viola.merge(ls_inputs=[delly, gridss, lumpy, manta], integration=True)
    print(merged_all.table_list)
    info_lst = ["positions", "mergedid", 'supportingid', 'supportingcaller', 'supportingidcount', 'supportingcallercount']
    for i in info_lst:
        print(merged_all.get_table(i))
        print(merged_all.get_table(i).dtypes)
        print()