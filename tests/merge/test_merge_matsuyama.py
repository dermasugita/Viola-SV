import viola    
import numpy as np
delly = viola.read_vcf("/shared_data/share/merge/test.merge.delly.vcf", variant_caller="delly")
gridss = viola.read_vcf("/shared_data/share/merge/test.merge.gridss.vcf", variant_caller="gridss")
lumpy = viola.read_vcf("/shared_data/share/merge/test.merge.lumpy.vcf", variant_caller="lumpy")
manta = viola.read_vcf("/shared_data/share/merge/test.merge.manta.vcf", variant_caller="manta")

def test_merge():
    merged_all = viola.merge(ls_inputs=[delly, gridss, lumpy, manta], integration=True)
    print(merged_all.table_list)
    info_lst = ["positions", "mergedid", 'supportingid', 'supportingcaller', 'supportingidcount', 'supportingcallercount']
    for i in info_lst:
        print(merged_all.get_table(i))
        print(merged_all.get_table(i).dtypes)
        print()