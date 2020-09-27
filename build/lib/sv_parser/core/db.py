from sv_parser.io.parser import read_vcf
class Sgt_core(object):
    def __init__(self, df_svpos, df_filters, dict_df_info, df_formats, ls_df_headers = []):
        self.df_svpos = df_svpos
        self.df_filters = df_filters
        self.dict_df_info = dict_df_info
        self.df_formats = df_formats
        self.ls_df_headers = df_headers
        ls_keys = ['POSITIONS', 'FILTERS'] + list(dict_df_info.keys()) + ['FORMATS'] + \
        ['CONTIGS_META', 'ALTS_META', 'INFOS_META', 'FORMATS_META', 'FILTERS_META']
        ls_values = [df_svpos, df_filters] + list(dict_df_info.values()) + ls_df_headers
        # self.dict_alltables is a {tablename: table} dictionary
        self.dict_alltables = {k: v for k, v in zip(ls_keys, ls_values)}

    def get_table_list(self):
        return list(self.dict_alltables.keys())
    def get_table(self, table_name):
        try:
            table = self.dict_alltables[table_name]
        except KeyError:
            print('no such table')
            return
        return table
    pass 

test = read_vcf('../../tests/vcf/manta1.inv.vcf')
print(test)
test2 = Sgt_core(*test)

print(test2.get_table_list())
test3 = test2.get_table('POSITIONS')
print(test3)
