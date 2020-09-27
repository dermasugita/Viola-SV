import pandas as pd
from sv_parser.io.parser import read_vcf
class Sgt_core(object):
    def __init__(self, df_svpos, df_filters, dict_df_info, df_formats, ls_df_headers = []):
        self.df_svpos = df_svpos
        self.df_filters = df_filters
        self.dict_df_info = dict_df_info
        self.df_formats = df_formats
        self.ls_df_headers = ls_df_headers
        ls_keys = ['positions', 'filters'] + [x.lower() for x in list(dict_df_info.keys())] + ['formats'] + \
        ['contigs_meta', 'alts_meta', 'infos_meta', 'formats_meta', 'filters_meta']
        ls_values = [df_svpos, df_filters] + list(dict_df_info.values()) + ls_df_headers
        # self.dict_alltables is a {tablename: table} dictionary
        self.dict_alltables = {k: v for k, v in zip(ls_keys, ls_values)}
        #print(ls_df_headers)

    def get_table_list(self):
        return list(self.dict_alltables.keys())
    def get_table(self, table_name):
        try:
            table = self.dict_alltables[table_name]
        except KeyError:
            print('no such table')
            return
        return table
    def to_bedpe_like(self, how='minimum'):
        df_svpos = self.get_table('positions')
        df_svpos.rename(columns={'pos1': 'start1', 'pos2': 'start2'}, inplace=True)
        df_homlen = self.get_table('homlen')
        df_merged = df_svpos.merge(df_homlen, how='left')
        df_merged['homlen1'] = df_merged['homlen1'].fillna(0).astype(int)
        df_merged['end1'] = df_merged['start1'] + df_merged['homlen1'] + 1# potential risk of name conflict!!!
        df_merged['end2'] = df_merged['start2'] + df_merged['homlen1'] + 1
        
        if how == 'minimum':
            df_out = df_merged[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'id', 'strand1', 'strand2']].copy()
        else:
            df_out = df_merged.copy()
        df_out.rename(columns={'id': 'name'}, inplace=True)
        
        print(df_out)
        

test = read_vcf('../../tests/vcf/manta1.inv.vcf')
test2 = Sgt_core(*test)

print(test2.get_table_list())
test3 = test2.get_table('infos_meta')
test2.to_bedpe_like()
print(test2.dict_alltables)
