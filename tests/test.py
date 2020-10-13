import sv_parser
import pandas as pd

lumpy = sv_parser.read_vcf('./tests/vcf/lumpy1_sample.vcf')
print(lumpy.to_bedpe_like(how='expand', add_formats=True).iloc[:, -8:].head(20))
lumpy.get_table_list()
lumpy_filtered = lumpy.filter(['1N SR == 0', '1N PE == 0', '1T SR > 0'])
print(lumpy_filtered.to_bedpe_like(how='expand', add_formats=True).iloc[:, -8:].head(20))
lumpy_all = sv_parser.read_vcf('./tests/vcf/lumpy1.vcf')
lumpy_somatic = lumpy_all.filter(['1N SR == 0', '1N PE == 0', '1T SR > 0'])
print(lumpy_somatic.to_bedpe_like(how='expand', add_formats=True).iloc[:, -8:])
