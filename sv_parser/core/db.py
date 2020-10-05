import pandas as pd
from sv_parser.io.parser import read_vcf
class Sgt_core(object):
    def __init__(self, df_svpos, df_filters, dict_df_info, df_formats, dict_df_headers = {}):
        self.df_svpos = df_svpos
        self.df_filters = df_filters
        self.dict_df_info = dict_df_info
        self.df_formats = df_formats
        self.dict_df_headers = dict_df_headers
        self.ls_infokeys = [ x.lower() for x in dict_df_info.keys() ]
        ls_keys = ['positions', 'filters'] + self.ls_infokeys + ['formats'] + \
        list(dict_df_headers.keys())
        ls_values = [df_svpos, df_filters] + list(dict_df_info.values()) + [df_formats] + list(dict_df_headers.values())
        # self.dict_alltables is a {tablename: table} dictionary
        self.dict_alltables = {k: v for k, v in zip(ls_keys, ls_values)}
        #print(ls_df_headers)

    def __repr__(self):
        return 'hello'
    def __str__(self):
        return 'hello'


    def get_table_list(self):
        return list(self.dict_alltables.keys())
    def get_table(self, table_name):
        try:
            table = self.dict_alltables[table_name]
        except KeyError:
            print('no such table')
            return
        return table.copy()
    def to_bedpe_like(self, how='basic', custom_infonames=[]):
        df_svpos = self.get_table('positions')
        df_svpos.rename(columns={'pos1': 'start1', 'pos2': 'start2'}, inplace=True)
        if 'homlen' in self.get_table_list():
            df_homlen = self.get_table('homlen')
            df_merged = df_svpos.merge(df_homlen, how='left')
            df_merged['homlen1'] = df_merged['homlen1'].fillna(0).astype(int)
            df_merged['end1'] = df_merged['start1'] + df_merged['homlen1'] + 1 # potential risk of name conflict!!!
            df_merged['end2'] = df_merged['start2'] + df_merged['homlen1'] + 1
        else:
            df_merged = df_svpos.copy()
            df_merged['end1'] = df_merged['start1'] + 1
            df_merged['end2'] = df_merged['start2'] + 1
        
        if how == 'basic':
            df_out = df_merged[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'id', 'qual', 'strand1', 'strand2']].copy()
        elif how == 'minimum':
            df_out = df_merged[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']].copy() 
        elif how == 'expand':
            if len(custom_infonames) == 0:
                print('error: specify columns to expand')
            df_out = df_merged[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'id', 'qual', 'strand1', 'strand2']].copy()
            df_out = self.append_infos(df_out, custom_infonames)
        df_out.rename(columns={'id': 'name', 'qual': 'score'}, inplace=True)
        return df_out
    def append_infos(self, base_df, ls_tablenames, left_on='id', auto_fillna=True):
        df = base_df.copy()
        if ('infos_meta' in self.get_table_list()) & auto_fillna:
            df_infometa = self.get_table('infos_meta')
            for tablename in ls_tablenames:
                df_to_append = self.get_table(tablename)
                df = pd.merge(df, df_to_append, how='left', left_on=left_on, right_on='id')
                info_dtype = df_infometa.loc[df_infometa['id']==tablename.upper(), 'type'].iloc[0]
                len_info = df_to_append.shape[1] - 1
                ls_ind_fancy = [tablename + str(i+1) for i in range(len_info)]
                if info_dtype == 'Integer':
                    df[ls_ind_fancy] = df[ls_ind_fancy].fillna(0).astype(int)
                elif info_dtype == 'Flag':
                    df[ls_ind_fancy] = df[ls_ind_fancy].fillna(False)
                if left_on != 'id':
                    df.drop('id', axis=1, inplace=True)
            return df        
        else:
            for tablename in ls_tablenames:
                df_to_append = self.get_table(tablename)
                df = pd.merge(df, df_to_append, how='left', left_on=left_on, right_on='id')
                if left_on != 'id':
                    df.drop('id', axis=1, inplace=True) 
            return df
    def append_formats(self, base_df, left_on='id'):
        df_format = self.get_table('formats')
        df_format['format_id'] = df_format['sample'] + '_' + df_format['format'] + '_' + df_format['value_idx'].astype(str) 
        df_format.drop(['sample', 'format', 'value_idx'], axis=1, inplace=True)
        df_format = df_format.pivot(index='id', columns='format_id', values='value')
        df_out = pd.merge(base_df, df_format, how='left', left_on=left_on, right_index=True)
        return df_out


    def append_filters(self, base_df, left_on='id'):
        df_filters = self.get_table('filters')
        df_filters_expand = df_filters['filter'].str.get_dummies()
        df_be_appended = pd.concat([ df_filters['id'], df_filters_expand ], axis=1).replace(to_replace={1: True, 0: False})
        df_out = pd.merge(base_df, df_be_appended, left_on=left_on, right_on='id')
        return df_out


    def _filter_by_id(self, tablename, arrlike_id):
        df = self.get_table(tablename)
        return df.loc[df['id'].isin(arrlike_id)].reset_index(drop=True)

    def filter_by_id(self, arrlike_id):
        out_svpos = self._filter_by_id('positions', arrlike_id)
        out_filters = self._filter_by_id('filters', arrlike_id)
        out_dict_df_info = {k: self._filter_by_id(k, arrlike_id) for k in self.ls_infokeys}
        out_formats = self._filter_by_id('formats', arrlike_id)
        out_dict_df_headers = self.dict_df_headers.copy()
        return Sgt_core(out_svpos, out_filters, out_dict_df_info, out_formats, out_dict_df_headers)

    def _filter_pos_table(self, item, operator, threshold):
        df = self.get_table('positions')
        e = "df.loc[df[item] {0} threshold]['id']".format(operator)
        return set(eval(e))

    def _filter_filters(self, _filter):
        df = self.get_table('filters')
        return set(df.loc[df['filter'] ==  _filter]['id'])

    def _filter_infos(self, infoname, value_idx=1, operator=None, threshold=None):## returning result ids
        df = self.get_table(infoname)
        value_idx = str(value_idx)
        e = "df.loc[df[''.join([infoname, value_idx])] {0} threshold]['id']".format(operator)
        return set(eval(e))

    def _filter_formats(self, sample, item, item_idx=1, operator=None, threshold=None):
        df = self.get_table('formats')
        target_q = (df['sample'] == sample) & (df['format'] == item) & (df['value_idx'] == item_idx)
        df_target = df.loc[target_q]
        e = "df_target.loc[df_target['value'] {0} threshold]['id']".format(operator)
        return set(eval(e))

    def _parse_filter_query(self, query):
        pass

    def filter(self, query):
        pass


        

pd.set_option('display.max_columns', 20)
pd.set_option('display.max_colwidth', 20)
pd.set_option('display.width', 1500) 
 


test = read_vcf('../../tests/vcf/manta1.inv.vcf')
test2 = Sgt_core(*test)
test3 = test2._filter_infos('homlen', operator='>', threshold=1)
test4 = test2._filter_infos('svtype', operator='==', threshold='INV') 
test5 = test2._filter_formats('mouse1_T', 'PR', 1, '>', 50)
test6 = test2._filter_filters('PASS')
test7 = test2._filter_pos_table('ref', '==', 'A')

test_filter = test2.filter_by_id(test7)
print(test_filter.get_table('positions'))
#print(test_filter.get_table('filters'))
#print(test_filter.to_bedpe_like('expand', custom_infonames=['homlen']))
#print(test_filter.get_table('formats'))
