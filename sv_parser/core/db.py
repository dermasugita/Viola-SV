import numpy as np
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
    def get_ids(self):
        df = self.get_table('positions')
        return set(df['id'])
    def to_bedpe_like(self, how='basic', custom_infonames=[], add_filters=False, add_formats=False, unique_events=False):
        if unique_events:
            df_svpos = self.get_unique_events().get_table('positions')
        else:
            df_svpos = self.get_table('positions')
        df_svpos.rename(columns={'pos1': 'start1', 'pos2': 'start2'}, inplace=True)
        if 'homlen' in self.get_table_list():
            df_homlen = self.get_table('homlen')
            df_merged = df_svpos.merge(df_homlen, how='left')
            df_merged['homlen0'] = df_merged['homlen0'].fillna(0).astype(int)
            df_merged['end1'] = df_merged['start1'] + df_merged['homlen0'] + 1 # potential risk of name conflict!!!
            df_merged['end2'] = df_merged['start2'] + df_merged['homlen0'] + 1
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
            if add_filters:
                df_out = self.append_filters(df_out)
            if add_formats:
                df_out = self.append_formats(df_out)
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
                ls_ind_fancy = [tablename + str(i) for i in range(len_info)]
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
        df_be_appended = pd.concat([ df_filters['id'], df_filters_expand ], axis=1)
        df_be_appended = df_be_appended.groupby('id').sum().replace(to_replace={1: True, 0: False})
        df_out = pd.merge(base_df, df_be_appended, how='left', left_on=left_on, right_on='id')
        return df_out
    

    def _parse_filter_query(self, q):
        sq = q.split(' ')

        # for flag informations and filters
        if sq[0].startswith('!'):
            sq0 = sq[0][1:]
        else:
            sq0 = sq[0]

        if sq0 in self.ls_infokeys:
            df_infometa = self.get_table('infos_meta')
            row_mask = df_infometa['id'].str.contains(sq0.upper())
            sq_dtype = df_infometa.loc[row_mask, 'type'].iloc[0]
            if sq_dtype == 'Integer':
                sq[-1] = int(sq[-1])
            elif sq_dtype == 'String':
                sq[-1] = str(sq[-1])
            elif sq_dtype =='Flag':
                if len(sq) == 1:
                    if sq[0].startswith('!'):
                        flag = False
                    else:
                        flag = True
                else:
                    flag = True if sq[-1] == 'True' else False
                exclude = not flag
                set_out = self._filter_infos_flag(sq0, exclude=exclude)
                return set_out

            if len(sq) == 3:
                set_out = self._filter_infos(sq[0], 0, sq[1], sq[2])
            else:
                set_out = self._filter_infos(*sq)
            #print(set_out)
            return set_out

        # is_filter?
        arr_filters = self.get_table('filters_meta')['id'].values
        ls_filters = list(arr_filters) + ['PASS']
        if sq0 in ls_filters:
            if len(sq) == 1:
                if sq[0].startswith('!'):
                    flag = False
                else:
                    flag = True
            else:
                flag = True if sq[-1] == 'True' else False
            exclude = not flag
            set_out = self._filter_filters(sq0, exclude=exclude)
            return set_out 

        # is_format?
        if sq0 in self.get_table('samples_meta').values:
            df_formatmeta = self.get_table('formats_meta')
            row_mask = df_formatmeta['id'].str.contains(sq[1])
            sq_dtype = df_formatmeta.loc[row_mask, 'type'].iloc[0]
            if sq_dtype == 'Integer':
                sq[-1] = int(sq[-1])
            elif sq_dtype == 'String':
                sq[-1] = str(sq[-1])
            if len(sq) == 4:
                set_out = self._filter_formats(sq[0], sq[1], 0, sq[2], sq[3])
            else:
                sq[2] = int(sq[2])
                set_out = self._filter_formats(*sq)
            return set_out

    def filter(self, ls_query, query_logic='and'):
        if isinstance(ls_query, str):
            ls_query = [q]
        if query_logic == 'and':
            set_result = self.get_ids()
            for query in ls_query:
                set_query = self._parse_filter_query(query)
                set_result = set_result & set_query
        out = self.filter_by_id(set_result)
        return out


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

    def _filter_filters(self, _filter, exclude=False):
        df = self.get_table('filters')
        set_out = set(df.loc[df['filter'] ==  _filter]['id'])
        if exclude:
            set_out = self.get_ids() - set_out
        return set_out

    def _filter_infos(self, infoname, value_idx=0, operator=None, threshold=None):## returning result ids
        df = self.get_table(infoname)
        value_idx = str(value_idx)
        e = "df.loc[df[''.join([infoname, value_idx])] {0} threshold]['id']".format(operator)
        set_out = set(eval(e))
        return set_out

    def _filter_infos_flag(self, infoname, exclude=False):
        df = self.get_table(infoname)
        set_out = set(df['id'])
        if exclude:
            set_out = self.get_ids() - set_out
        return set_out
        
    def _filter_formats(self, sample, item, item_idx=0, operator=None, threshold=None):
        df = self.get_table('formats')
        target_q = (df['sample'] == sample) & (df['format'] == item) & (df['value_idx'] == item_idx)
        df_target = df.loc[target_q]
        e = "df_target.loc[df_target['value'] {0} threshold]['id']".format(operator)
        return set(eval(e))

    def get_unique_events(self):
        if 'mateid' not in self.get_table_list():
            print("Can't find mateid table")
            return
        df = self.get_table('mateid')
        df2 = df.reset_index().set_index('mateid0').loc[df['id']]
        arr_mask = df2['index'].values > np.arange(df2.shape[0])
        set_to_subtract = set(df2.loc[arr_mask]['id'])
        set_all_ids = self.get_ids()
        set_result_ids = set_all_ids - set_to_subtract
        return self.filter_by_id(set_result_ids)

    def get_detail_svtype(self):
        if 'svtype' not in self.get_table_list():
            print("Can't find svtype table")
            return

        pass


pd.set_option('display.max_columns', 20)
pd.set_option('display.max_colwidth', 35)
pd.set_option('display.width', 1500) 
 
test = read_vcf('../../tests/vcf/manta1.inv.vcf')
test2 = Sgt_core(*test)
q = "mouse1_T SR 1 > 0"
q2 = "PASS"
q3 = "mouse1_N PR 1 == 0"
q4 = "mouse1_N SR 1 == 0"
q5 = "svlen > 10000"
test8 = test2.filter([q, q2, q3, q4, q5])
print(test8.get_table('formats_meta')['description'].values)
#test_filter2 = test_filter.get_unique_events()
#print(test2.get_table('samples_meta'))
#print(test2.get_table_list())
print(test8.to_bedpe_like(how='expand', custom_infonames=['svtype', 'svlen'], add_filters=True, add_formats=True, unique_events=True))
#print(test_filter.get_table('filters'))
#print(test_filter.to_bedpe_like('expand', custom_infonames=['homlen']))
#print(test_filter.get_table('formats'))
