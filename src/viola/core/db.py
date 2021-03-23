import os,sys
import numpy as np
import pandas as pd
import re
from functools import reduce
from typing import (
    List,
    Set,
    Iterable,
)
from collections import OrderedDict
import viola
from viola.core.indexing import Indexer
from viola.core.bed import Bed
from viola.core.fasta import Fasta
from viola.utils.microhomology import get_microhomology_from_positions
from viola.utils.utils import get_inslen_and_insseq_from_alt
from viola._typing import (
    IntOrStr,
    StrOrIterableStr,
)
from viola._exceptions import (
    TableNotFoundError,
    InfoNotFoundError,
    ContigNotFoundError,
)

from sklearn.cluster import AgglomerativeClustering

class Bedpe(Indexer):
    """
    Relational database-like object containing SV position dataframes and INFO dataframes.
    The instances of this class have information equal to the BEDPE files.
    ...

    Attributes
    ----------
    sv_count: int
        Number of SV records
    table_list
        List of names of all tables included in the object 
    ids
        List of all SV id.
    
    Parameters
    ----------
    df_svpos: DataFrame
        DataFrame containing information such as position, strand, svtype, etc.
        Columns should be following:
        ['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype']
        Main key is 'id'.
    odict_df_info: dict[str, DataFrame]
        OrderedDict of DataFrames which contain additional information on SV record (equivalent to INFO field of vcf).
        Each item of the dictionary contains single INFO.
        The dictionary key is the name of each INFO and should be in lowercase.
        Columns of the DataFrame should be following:
        ['id', 'value_idx', 'infoname']
        The 'value_idx' column contains 0-origin indice of INFO values.
        This is important when one SV record has multiple values of an INFO (eg. cipos). 
        Main key is the combination of ('id', 'value_idx') and 'id' is the foreign key coming from df_svpos table.

    Methods
    ----------
    get_table(table_name)
        Return a table specified in the argument as pandas DataFrame object.
    to_bedpe_like(custom_infonames=[], confidence_intervals=False)
        Return a DataFrame in bedpe-like format.
    filter(ls_query, query_logic="and")
        Filter Bedpe object by the list of queries.
        Return object is also an instance of the Bedpe object
    filter_by_id(arrlike_id)
        Filter Bedpe object according to the list of SV ids.
        Return object is also an instance of the Bedpe object

    """
    _internal_attrs = [
        "_df_svpos",
        "_odict_df_info",
        "_ls_infokeys",
        "_odict_alltables",
        "_repr_config",
        "_sig_criteria"
    ]
    _internal_attrs_set = set(_internal_attrs)
    _repr_column_names = [
        "id",
        "be1",
        "be2",
        "strand",
        "qual",
        "svtype",
    ]
    _repr_column_names_set = set(_repr_column_names)

    def __init__(self, df_svpos: pd.DataFrame, odict_df_info: 'OrderedDict[str, pd.DataFrame]'):
        if not isinstance(odict_df_info, OrderedDict):
            raise TypeError('the type of the argument "odict_df_info" should be collections.OrderedDict')
        self._df_svpos = df_svpos
        self._odict_df_info = odict_df_info
        self._ls_infokeys = [x.lower() for x in odict_df_info.keys()]
        ls_keys = ['positions'] + self._ls_infokeys
        ls_values = [df_svpos] + list(odict_df_info.values())
        self._odict_alltables = OrderedDict([(k, v) for k, v in zip(ls_keys, ls_values)])
        self._repr_config = {
            'info': None,
        }
    
    @property
    def sv_count(self) -> int:
        """
        Return number of SV records.
        """
        return self._df_svpos.shape[0]
    
    @property
    def table_list(self) -> List[str]:
        """
        Return a list of names of all tables in the object. 
        """
        return list(self._odict_alltables.keys())
    
    @property
    def ids(self):
        """
        Return all SV ids as list.
        """
        return list(self.get_ids())
    
    @property
    def contigs(self) -> List[str]:
        """
        Return a list of contigs(chromosomes) included in the object.
        """
        df_svpos = self.get_table('positions')
        arr_chrom1 = df_svpos['chrom1'].unique()
        arr_chrom2 = df_svpos['chrom2'].unique()
        arr_chrom = np.unique(np.concatenate((arr_chrom1, arr_chrom2)))
        return list(arr_chrom)
    
    @property
    def repr_config(self):
        """
        Return current configuration of __repr__() function.
        """
        return self._repr_config

    def add_info_table(self, table_name: str, df: pd.DataFrame):
        """
        add_info_table(table_name, df)
        Add a new INFO table to self.

        Parameters
        ------------
        table_name: str
            The name of the table to be added.
        df: DataFrame
            The table (pandas DataFrame) to be added.
            The column names should be following:
            1st column: "id"
            2nd column: "value_idx"
            3rd column: table_name (equivalent to the first argument of this function)
        """
        self._ls_infokeys += [table_name]
        self._odict_alltables[table_name] = df
    
    def change_repr_config(self, key, value):
        self._repr_config[key] = value

    def __repr__(self):
        config_info = self._repr_config['info']
        return self.view(custom_infonames=config_info)

    def view(self, custom_infonames=None, return_as_dataframe=False):
        """
        view(custom_infonames, return_as_dataframe)
        Quick view function of the Vcf object.

        Parameters
        -----------
        custom_infonames: list_like or None, default None
            The names of the INFO to show additionally.
        return_as_dataframe: bool, default False
            If true, return as pandas DataFrame.
        """
        df_svpos = self.get_table('positions')
        ser_id = df_svpos['id']
        ser_be1 = df_svpos['chrom1'].astype(str) + ':' + df_svpos['pos1'].astype(str)
        ser_be2 = df_svpos['chrom2'].astype(str) + ':' + df_svpos['pos2'].astype(str)
        ser_strand = df_svpos['strand1'] + df_svpos['strand2']
        ser_qual = df_svpos['qual']
        ser_svtype = df_svpos['svtype']
        ls_ser = [ser_id, ser_be1, ser_be2, ser_strand, ser_qual, ser_svtype]
        ls_key = ['id', 'be1', 'be2', 'strand', 'qual', 'svtype']
        dict_ = {k: v for k, v in zip(ls_key, ls_ser)}
        df_out = pd.DataFrame(dict_)
        if custom_infonames is not None:
            df_out = self.append_infos(df_out, ls_tablenames=custom_infonames)
        str_df_out = str(df_out)
        str_infokeys = ','.join(list(self._ls_infokeys))
        desc = 'INFO='
        out = desc + str_infokeys + '\n' + str_df_out
        if return_as_dataframe:
            return df_out
        return str(out)
    
    def __getattr__(self, value):
        if value in self._internal_attrs_set:
            return object.__getattribute__(self, value)
        else:
            return self.get_table(value)
    
    def __getitem__(self, value):
        if value in self._repr_column_names:
            return self.view(return_as_dataframe=True)[value]
        return self.get_table(value)



    def get_table(self, table_name: str) -> pd.DataFrame:
        """
        get_table(table_name: str)
        Return a table specified in the argument as pandas DataFrame object.

        Parameters
        ----------
        table_name: str
            The name of the table to return.

        Returns
        ----------
        DataFrame
            A table specified in the table_name argument.
        
        Raises
        ----------
        TableNotFoundError
            If the table_name doesn't exist in the object.
        """
        if table_name not in self.table_list:
            raise TableNotFoundError(table_name)
        table = self._odict_alltables[table_name]
        return table.copy()
    
    def replace_table(self, table_name: str, table: pd.DataFrame):
        """
        replace_table(table_name, table)
        Replace existing table into new table.
        
        Parameters
        ------------
        table_name: str
            The name of the table to be replaced
        table: DataFrame
            A new table to be set
        """
        if table_name not in self.table_list:
            raise TableNotFoundError(table_name)
        self._odict_alltables[table_name] = table


    def get_ids(self) -> Set[IntOrStr]:
        """
        Return all SV ids as the set type.
        """
        df = self.get_table('positions')
        return set(df['id'])

    def to_bedpe_like(
        self,
        custom_infonames: Iterable[str] = [],
        confidence_intervals: bool = False
    ) -> pd.DataFrame:
        """
        to_bedpe_like(custom_infonames=[], confidence_intervals: bool=False)
        Return a DataFrame in bedpe-like format.
        When specified, you can add INFOs as additional columns.

        Parameters
        ----------
        custom_infonames: list-like[str]
            The table names of INFOs to append.
        confidence_intervals: bool, default False
            Whether or not to consider confidence intervals of the breakpoints.  
            If True, confidence intervals for each breakpoint are represented by [start1, end1) and [start2, end2), respectively.
            Otherwise, breakpoints are represented by a single-nucleotide resolution.
        
        Returns
        ----------
        DataFrame
            A Dataframe in bedpe-like format.
            The columns include at least the following:
            ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2']
        """
        df_svpos = self.get_table('positions')
        if confidence_intervals:
            if 'cipos' in self.table_list and 'ciend' in self.table_list:
                df_svpos = self.append_infos(df_svpos, ['cipos', 'ciend'])
                df_svpos['start1'] = df_svpos['pos1'] - df_svpos['cipos_0'] - 1
                df_svpos['end1'] = df_svpos['pos1'] + df_svpos['cipos_1']
                df_svpos['start2'] = df_svpos['pos2'] - df_svpos['ciend_0'] - 1
                df_svpos['end2'] = df_svpos['pos2'] + df_svpos['ciend_1']
            else:
                pass # raise some exception
        else:
            df_svpos.rename(columns={'pos1': 'end1', 'pos2': 'end2'}, inplace=True)
            df_svpos['start1'] = df_svpos['end1'] - 1
            df_svpos['start2'] = df_svpos['end2'] - 1
        
        df_out = df_svpos[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'id', 'qual', 'strand1', 'strand2']].copy()
        if len(custom_infonames) != 0:
            df_out = self.append_infos(df_out, custom_infonames)
        df_out.rename(columns={'id': 'name', 'qual': 'score'}, inplace=True)
        return df_out

    def append_infos(self,
        base_df: pd.DataFrame,
        ls_tablenames: Iterable[str],
        left_on: str = 'id',
        auto_fillna: bool = False) -> pd.DataFrame:
        """
        append_infos(base_df, ls_tablenames, left_on='id')
        Append INFO tables to the right of the base_df, based on the SV id columns.
        If the name of the SV id column in base_df is not 'id', specify column name into left_on argument. 

        Parameters
        ---------------
        base_df: DataFrame
            The DataFrame to which the INFO tables are appended.
        ls_tablenames: list-like
            The list of INFO table names to be appended.
        left_on: str
            The name of SV id column of base_df
        
        Returns
        ---------------
        DataFrame
            A DataFrame which the INFO tables are added.
        
        Note
        ----
        This function seems to have a bug when a value is 
        passed to the ls_tablenames argument.
        Can someone please fix it?
        """
        df = base_df.copy()
        for tablename in ls_tablenames:
            df_to_append_pre = self.get_table(tablename)
            df_to_append_pre['new_column_names'] = tablename + '_' + df_to_append_pre['value_idx'].astype(str)
            df_to_append = df_to_append_pre.pivot(index='id', columns='new_column_names', values=tablename)
            df = pd.merge(df, df_to_append, how='left', left_on=left_on, right_on='id')
            if left_on != 'id':
                df.drop('id', axis=1, inplace=True) 
            if pd.api.types.is_bool_dtype(df_to_append_pre.iloc[:, 2]):
                column_names = np.unique(df_to_append_pre['new_column_names'])
                df[column_names] = df[column_names].fillna(False)
        return df

    def _parse_filter_query(self, q):
        # sq: split query
        sq = q.split(' ')

        # for flag informations and filters
        if sq[0].startswith('!'):
            sq0 = sq[0][1:]
        else:
            sq0 = sq[0]

        if sq0 in self._ls_infokeys:
            sqtail = sq[-1]
            def _isfloat(value):
                try:
                    float(value)
                    return True
                except ValueError:
                    return False
            def _isflag(value):
                if value == 'True':
                    return True
                elif value == 'False':
                    return True
                else:
                    return False
            if sqtail.isdigit():
                sq_dtype = 'Integer'
            elif _isfloat(sqtail):
                sq_dtype = 'Float'
            elif _isflag(sqtail):
                sq_dtype = 'Flag'
            else:
                sq_dtype = 'String'

            if sq_dtype == 'Integer':
                sq[-1] = int(sq[-1])
            elif sq_dtype == 'Float':
                sq[-1] = float(sq[-1])
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
        
        # is_locus?
        if sq0 in ['be1', 'be2', 'pos1', 'pos2']:
            split_locus = sq[1].split(':')
            chrom = split_locus[0]
            if chrom.startswith('!'):
                exclude_flag = True
                chrom = chrom[1:]
            else:
                exclude_flag = False
            
            if chrom not in self.contigs:
                raise ContigNotFoundError(chrom)

            if len(split_locus) == 1:
                st = None
                en = None
            elif len(split_locus) == 2:
                split_locus_coord = split_locus[1].split('-')
                if len(split_locus_coord) == 1:
                    st = int(split_locus_coord[0])
                    en = int(split_locus_coord[0]) + 1
                elif len(split_locus_coord) == 2:
                    st = split_locus_coord[0]
                    en = split_locus_coord[1]
                    if st == '':
                        st = None
                    else:
                        st = int(st)
                    if en == '':
                        en = None
                    else:
                        en = int(en)
            
            if sq0 in ['be1', 'pos1']:
                pos_num = 1
            elif sq0 in ['be2', 'pos2']:
                pos_num = 2
            
            args = [pos_num, chrom, st, en]

            if exclude_flag:
                return self._filter_by_positions_exclude(*args)
            
            return self._filter_by_positions(*args)
            
            

    
    def _filter(self, ls_query, query_logic):
        if isinstance(ls_query, str):
            ls_query = [ls_query]
        if query_logic == 'and':
            set_result = self.get_ids()
            for query in ls_query:
                set_query = self._parse_filter_query(query)
                set_result = set_result & set_query
        elif query_logic == 'or':
            set_result = set()
            for query in ls_query:
                set_query = self._parse_filter_query(query)
                set_result = set_result | set_query
        
        return set_result
    
    def get_info(self, info_name: str) -> pd.DataFrame:
        """
        get_info(info_name: str)
        Return a info specified in the argument as pandas DataFrame object.

        Parameters
        ----------
        info_name: str
            The name of the info to return.

        Returns
        ----------
        DataFrame
            A info specified in the info_name argument.
        
        Raises
        ----------
        InfoNotFoundError
            If the info_name doesn't exist in the object.
        """

        if info_name not in self._ls_infokeys:
            raise InfoNotFoundError(info_name)
        return self.get_table(info_name)
        

    def filter(self,
        ls_query: StrOrIterableStr,
        query_logic: str = 'and'):
        """
        filter(ls_query, query_logic)
        Filter Bedpe object by the list of queries.
        Return object is also an instance of the Bedpe object

        Parameters
        ----------
        ls_query: str or List[str]
            A query or a list of query.
        query_logic: str, default 'and'
            Two options are allowed; ('and', 'or').
            If 'and' is specified, SV records that meet all the queries will be returned.
            If 'or' is specified, SV records that meet at lease one query will be returned.
        
        Returnes
        ----------
        Bedpe
            A Bedpe object that includes SV records filtered by queries.
        """
        ### != operation is dangerous
        set_result = self._filter(ls_query, query_logic)
        out = self.filter_by_id(set_result)
        return out


    def _filter_by_id(self, tablename, arrlike_id):
        """
        _filter_by_id(tablename, arrlike_id)
        Filter pandas DataFrame by SV ids.
        Input DataFrame should have a columns named 'id'.

        Parameters
        ----------
        tablename: str
            The name of the table to be filtered.
        arrlike_id: list-like
            SV ids which you want to keep.
        
        Returns
        --------
        A filtered DataFrame.
        """
        df = self.get_table(tablename)
        return df.loc[df['id'].isin(arrlike_id)].reset_index(drop=True)


    def filter_by_id(self, arrlike_id):
        """
        filter_by_id(arrlike_id)
        Filter Bedpe object according to the list of SV ids.
        Return object is also an instance of the Bedpe object

        Parameters
        ---------------
        arrlike_id: list-like
            SV ids which you want to keep.
        
        Returns
        ---------------
        Bedpe
            A Bedpe object with the SV id specified in the arrlike_id argument.
            All records associated with SV ids that are not in the arrlike_id will be discarded.
        """
        out_svpos = self._filter_by_id('positions', arrlike_id)
        out_odict_df_info = OrderedDict([(k, self._filter_by_id(k, arrlike_id)) for k in self._ls_infokeys])
        return Bedpe(out_svpos, out_odict_df_info)

    def _filter_pos_table(self, item, operator, threshold):
        df = self.get_table('positions')
        e = "df.loc[df[item] {0} threshold]['id']".format(operator)
        return set(eval(e))

    def _filter_infos(self, infoname, value_idx=0, operator=None, threshold=None):## returning result ids
        df = self.get_table(infoname)
        value_idx = int(value_idx)
        df = df.loc[df['value_idx'] == value_idx]
        e = "df.loc[df[infoname] {0} threshold]['id']".format(operator)
        set_out = set(eval(e))
        return set_out

    def _filter_infos_flag(self, infoname, exclude=False):
        df = self.get_table(infoname)
        set_out = set(df['id'])
        if exclude:
            set_out = self.get_ids() - set_out
        return set_out
    
    def annotate_bed(self,
        bed: Bed,
        annotation: str,
        how: str = 'flag',
        suffix=['left', 'right']):
        """
        annotate_bed(bed, annotation, suffix=['left', 'right'])
        Annotate SV breakpoints using Bed class object.
        Annotation is stored as INFO table.
        For each SV record, two annotations will be made.
        
        Parameters
        -----------
        bed: Bed
            A Bed class object for annotation.
        annotation: str
            The label of annotation.
            The suffixes will be attached then added as INFO table.
        how: str ['flag', 'value'], default 'flag'
            If 'flag', Annotate True when a breakend is in Bed, otherwise False.
            If 'value', Annotate values in the Bed.
        suffix: List[str], default ['left', 'right']
            The suffix that attached after annotation label specified above.
        """
        df_svpos = self.get_table('positions')
        ls_left = []
        ls_right = []
        for idx, row in df_svpos.iterrows():
            svid = row['id']
            chrom1 = row['chrom1']
            pos1 = row['pos1']
            chrom2 = row['chrom2']
            pos2 = row['pos2']

            df_bp1 = bed.query(chrom1, pos1)
            df_bp2 = bed.query(chrom2, pos2)

            if how == 'flag':
                if not df_bp1.empty:
                    ls_left.append([svid, 0, True])
                if not df_bp2.empty:
                    ls_right.append([svid, 0, True])
            elif how == 'value':
                if not df_bp1.empty:
                    j = 0
                    for idx_inner, query_result in df_bp1.iterrows():
                        ls_left.append([svid, j, query_result['name']])
                        j += 1
                if not df_bp2.empty:
                    j = 0
                    for idx_inner, query_result in df_bp2.iterrows():
                        ls_right.append([svid, j, query_result['name']])
                        j += 1
        left_name = annotation + suffix[0]
        right_name = annotation + suffix[1]
        df_left = pd.DataFrame(ls_left, columns=('id', 'value_idx', left_name))
        df_right = pd.DataFrame(ls_right, columns=('id', 'value_idx', right_name))
        self.add_info_table(left_name, df_left)
        self.add_info_table(right_name, df_right)

    def get_microhomology(self, fasta, max_homlen=200):
        """
        get_microhomology(fasta, max_homlen=200)
        Infer microhomology length and sequence in each breakpoint.
        The results will be appended as 'HOMLEN' and 'HOMSEQ' INFO, respectively.

        Parameters
        ----------
        fasta: Fasta
            A Fasta object.
        max_homlen: int
            Maximum length of microhomology to be considered.
        """
        df_svpos = self.get_table('positions')
        ls_homlen = []
        ls_homseq = []
        for idx, row in df_svpos.iterrows():
            svid = row['id']
            args = [
                row['chrom1'],
                row['pos1'],
                row['chrom2'],
                row['pos2'],
                row['strand1'],
                row['strand2'],
                fasta,
                max_homlen
            ]
            homlen, homseq = get_microhomology_from_positions(*args)
            ls_homlen.append([svid, 0, homlen])
            ls_homseq.append([svid, 0, homseq.upper()])
        df_homlen = pd.DataFrame(ls_homlen, columns=('id', 'value_idx', 'homlen'))
        df_homseq = pd.DataFrame(ls_homseq, columns=('id', 'value_idx', 'homseq'))
        self.add_info_table('homlen', df_homlen)
        self.add_info_table('homseq', df_homseq)
    
    def calculate_info(self, operation, name):
        """
        calculate_info(operation, name)
        Calculate values of INFO tables according to the 'operation' argument and add a new INFO table as the result.
        """
        ls_matched = re.findall(r'\${[^}]*}', operation)
        ls_infonames = [s[2:-1] for s in ls_matched]
        ls_df_info = []
        operation_replaced = operation
        for idx, infoname in enumerate(ls_infonames):
            # check the info names are valid
            if infoname not in self.table_list:
                raise TableNotFoundError(infoname)
            ls_df_info.append(self.get_table(infoname))
            operation_replaced = operation_replaced.replace(ls_matched[idx], 'df_merged["'+infoname+'"]')
        df_merged = reduce(lambda left, right: pd.merge(left, right, on=['id', 'value_idx']), ls_df_info)
        ser_result = eval(operation_replaced)
        df_merged[name] = ser_result
        df_to_add = df_merged[['id', 'value_idx', name]]
        self.add_info_table(name, df_to_add)
        



    def classify_manual_svtype(self, ls_conditions, ls_names, ls_order=None, return_series=True):
        """
        classify_manual_svtype(ls_conditions, ls_names, ls_order=None)
        Classify SV records by user-defined criteria. A new INFO table named
        'manual_sv_type' will be created.
        """
        set_ids_current = set(self.ids)
        obj = self
        ls_ids = []
        ls_result_names = []
        for func, name in zip(ls_conditions, ls_names):
            obj = obj.filter_by_id(set_ids_current)
            ids = func(obj)
            set_ids = set(ids)
            set_ids_intersection = set_ids_current & set_ids
            ls_ids += list(set_ids_intersection)
            ls_result_names += [name for i in range(len(set_ids_intersection))]
            set_ids_current = set_ids_current - set_ids_intersection
        ls_ids += list(set_ids_current)
        ls_result_names += ['others' for i in range(len(set_ids_current))]
        ls_zeros = [0 for i in range(len(self.ids))]
        df_result = pd.DataFrame({'id': ls_ids, 'value_idx': ls_zeros, 'manual_sv_type': ls_result_names})
        self.add_info_table('manual_sv_type', df_result)
        if return_series:
            if ls_order is None:
                pd_ind_reindex = pd.Index(ls_names + ['others'])
            else:
                pd_ind_reindex = pd.Index(ls_order)
            ser_feature_counts = self.get_feature_count_as_series(ls_order=pd_ind_reindex)
            return ser_feature_counts
    
    def get_feature_count_as_series(self, feature='manual_sv_type', ls_order=None):
        ser_feature_counts = self.get_table(feature)[feature].value_counts()
        if ls_order is not None:
            pd_ind_reindex = pd.Index(ls_order)
            ser_feature_counts = ser_feature_counts.reindex(index=pd_ind_reindex, fill_value=0)
        return ser_feature_counts



    def is_reciprocal(self):
        pass

    def _filter_by_positions(self, position_num, chrom, pos_min=None, pos_sup=None):
        """
        _filter_by_positions(position_num:int, chrom:str, pos_min:int, pos_sup:int)
        Return ids specified in the argument as a set

        Parameters
        ---------------
        positions_num:int
            1 for the first breakend and 2 for the other 
        chrom:str
            "chr1","chr2",...,"chrX","chrY"
        pos_min:int
        pos_sup:int

        Returns
        ---------------
        set
            A set of ids which satisfies the argument
        """
        positions_df = self.get_table("positions")
        positions_df = positions_df[positions_df["chrom{}".format(position_num)]==chrom]
        pos1or2 = "pos{}".format(position_num)
        if (pos_min is not None) and (pos_sup is not None):
            positions_df = positions_df[(pos_min <= positions_df[pos1or2])
            & (positions_df[pos1or2] < pos_sup)]
        elif pos_min is not None:
            positions_df = positions_df[pos_min <= positions_df[pos1or2]]
        elif pos_sup is not None:
            positions_df = positions_df[positions_df[pos1or2] < pos_sup]
        id_list = positions_df["id"].values
        id_set = set(id_list)
        return id_set
    
    def _filter_by_positions_exclude(self, ex_position_num, ex_chrom, ex_pos_min=None, ex_pos_max=None):
        """
        _filter_by_positions_exclude(position_num, chrom_exclude, ex_pos_min, ex_pos_max)
        Return a set of ids except which are specified by the argument
        
        Parameters
        ---------------
        ex_positions_num:int
            1 for the first breakend and 2 for the other
        ex_chrom:str
            "chr1","chr2",...,"chrX","chrY"
        ex_pos_min:int
        ex_pos_max:int

        Returns
        ---------------
        set
            A set of ids except which satisfies the argument
        """
        positions_df = self.get_table("positions")
        whole_id = positions_df["id"].values
        whole_id_set = set(whole_id)
        ex_positions_df = positions_df[positions_df["chrom{}".format(ex_position_num)]==ex_chrom]
        pos1or2 = "pos{}".format(ex_position_num)
        if (ex_pos_min is not None) and (ex_pos_max is not None):
            ex_positions_df = ex_positions_df[(ex_pos_min <= ex_positions_df[pos1or2])
            & (ex_positions_df[pos1or2] <= ex_pos_max)]
        elif ex_pos_min is not None:
            ex_positions_df = ex_positions_df[ex_pos_min <= ex_positions_df[pos1or2]]
        elif ex_pos_max is not None:
            ex_positions_df = ex_positions_df[ex_positions_df[pos1or2] <= ex_pos_max]
        ex_id = ex_positions_df["id"].values
        ex_id_set = set(ex_id)
        id_set = whole_id_set - ex_id_set
        return id_set


    def _nonoverlap(self, param):
        chr = [param["chr1h"], param["chr2h"]]
        pos1 = [param["pos1h"], param["pos1w"]]
        pos2 = [param["pos2h"], param["pos2w"]]
        proposition_chr = chr[0] == chr[1]
        proposition_pos = (max(pos1[0], pos2[0]) < min(pos1[1], pos2[1])) or (max(pos1[1], pos2[1]) < min(pos1[0], pos2[0]))
        proposition = proposition_chr and proposition_pos
        return proposition

    def _necessary_condition4merge(self, param, mode, str_missing=True):
        if mode == "normal":
            chr1 = [param["chr1h"], param["chr1w"]]
            chr2 = [param["chr2h"], param["chr2w"]]
            pos1 = [param["pos1h"], param["pos1w"]]
            pos2 = [param["pos2h"], param["pos2w"]]
            str1 = [param["str1h"], param["str1w"]]
            str2 = [param["str2h"], param["str2w"]]
        elif mode == "reverse":
            chr1 = [param["chr1h"], param["chr2w"]]
            chr2 = [param["chr2h"], param["chr1w"]]
            pos1 = [param["pos1h"], param["pos2w"]]
            pos2 = [param["pos2h"], param["pos1w"]]
            str1 = [param["str1h"], param["str2w"]]
            str2 = [param["str2h"], param["str1w"]]

        proposition_chr = (chr1[0] == chr1[1]) and (chr2[0] == chr2[1])
        proposition_pos = (pos1[0] == pos1[1]) and (pos2[0] == pos2[1])
        if str_missing: 
            proposition_str = ((str1[0] == str1[1]) or (str1[0]==".") or (str1[1]==".")) and ((str2[0] == str2[1]) or (str2[0]==".") or (str2[1]=="."))
        else:
            proposition_str = (str1[0] == str1[1]) and (str2[0] == str2[1])

        proposition = proposition_chr and proposition_pos and proposition_str
        return proposition


    def merge(self, ls_caller_names, threshold, ls_bedpe=[], linkage = "complete", str_missing=True):
        """
        merge(ls_caller_names:list, threshold:float, ls_bedpe=[], linkage = "complete", str_missing=True)
        Return a merged bedpe object from several caller's bedpe objects

        Parameters
        ----------
        ls_caller_names:list
            a list of names of bedpe objects to be merged, which should have self's name as the first element
        threshold:float
            Two SVs whose diference of positions is under this threshold are cosidered to be the same.
        ls_bedpe:list
            a list of bedpe objects to be merged
        linkage:{‘ward’, ‘complete’, ‘average’, ‘single’}, default=’complete’
            the linkage of hierachial clustering
        str_missing:boolean, default="True"
            If True, all the missing strands are considered to be identical to the others. 

        Returns
        ----------
        A Bedpe object
            A set of ids except which satisfies the argument
        """
        if self in ls_bedpe:
            pass
        else:
            ls_bedpe = [self] + ls_bedpe#ls_bedpeにはself入っていなくても良い

        multibedpe = viola.MultiBedpe(ls_bedpe, ls_caller_names)
        positions_table = multibedpe.get_table("positions")
        N = len(positions_table)#the number of samples
        penalty_length = 3e9#whole genome length
        distance_matrix = np.full((N,N), penalty_length)
        columns = ["chrom1", "chrom2", "pos1", "pos2", "strand1", "strand2"]
        for h in range(N):
            for w in range(N):
                param = {}
                for col in columns:#want to create these variables dynamically
                    key_h = col[:3] + col[-1] + "h"
                    value_h = positions_table.at[positions_table.index[h], col]
                    key_w = col[:3] + col[-1] + "w"
                    value_w = positions_table.at[positions_table.index[w], col]
                    param[key_h] = value_h
                    param[key_w] = value_w

                if self._necessary_condition4merge(param = param, mode = "normal", str_missing = str_missing): 
                    if self._nonoverlap(param=param):
                        distance_matrix[h, w] = penalty_length
                    else:
                        distance_matrix[h, w] = max(np.abs(param["pos1h"] - param["pos1w"]), np.abs(param["pos2h"] - param["pos2w"]))                          

                elif self._necessary_condition4merge(param = param, mode = "reverse", str_missing = str_missing):
                    if self._nonoverlap(param=param):
                        distance_matrix[h, w] = penalty_length
                    else:
                        distance_matrix[h, w] = max(np.abs(param["pos1h"] - param["pos2w"]), np.abs(param["pos2h"] - param["pos1w"]))
        
        hcl_clustering_model = AgglomerativeClustering(n_clusters=None, affinity="precomputed", linkage=linkage, distance_threshold=threshold)
        labels = hcl_clustering_model.fit_predict(X = distance_matrix)
        
        bpid_dict = {labels[0]:0}
        ls_bpid = []
        idx_head = 0
        for label in labels:
            if label in bpid_dict:
                ls_bpid.append(bpid_dict[label])
            else:
                idx_head += 1
                bpid_dict[label] = idx_head
                ls_bpid.append(bpid_dict[label])

        value_idx = pd.Series(np.zeros(N, dtype=int))
        df_bpid = pd.DataFrame({"id":positions_table["id"],"value_idx":value_idx, "bpid":pd.Series(ls_bpid)})
        
        originalid = multibedpe.get_table("global_id")["id"]
        df_originalid = pd.DataFrame({"id":positions_table["id"], "value_idx":value_idx, "originalid":originalid})
        
        caller = multibedpe.get_table("global_id")["patients"]
        df_caller = pd.DataFrame({"id":positions_table["id"], "value_idx":value_idx, "caller":caller})

        ls_infokeys = ['svlen', 'svtype', 'cipos', 'ciend']
        ls_df_infos = []
        for i in ls_infokeys:
            ls_df_infos.append(multibedpe.get_table(i))
        ls_infokeys = ls_infokeys + ["bpid", "originalid", "caller"]
        ls_df_infos = ls_df_infos + [df_bpid, df_originalid, df_caller]
        odict_df_infos = OrderedDict([(k, v) for k, v in zip(ls_infokeys, ls_df_infos)])
        args = [positions_table, odict_df_infos]
        merged_bedpe = viola.Bedpe(*args)
        
        return merged_bedpe

class Vcf(Bedpe):
    """
    Relational database-like object containing SV position dataframes,
    FILTER dataframe, INFO dataframes, FORMAT dataframe, and HEADER dataframes.
    The instances of this class have information equal to the VCF files.

    Attributes
    ---------------
    sv_count: int
        Number of SV records
    table_list
        List of names of all tables included in the object 
    ids
        List of all SV id.

    Parameters
    ----------
    df_svpos: DataFrame
        DataFrame containing information such as position, strand, svtype, etc.
        Columns should be following:
        ['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype']
        Main key is 'id'. The 'chrom1' and 'chrom2' are the foreign key from contigs_meta table.
    df_filter: DataFrame
        DataFrame containing FILTER information which locate on the 7th column of the vcf file.
        Columns of the input DataFrame should be following:
        ['id', 'filter']
        Main Key is the combination of ('id', 'filter'). Each column is the foreign key from 
        df_svpos, and filters_meta table, respectively.
    odict_df_info: dict[str, DataFrame]
        OrderedDict of DataFrames which contain additional information on SV record (equivalent to INFO field of vcf).
        Each item of the dictionary contains single INFO.
        The dictionary key is the name of each INFO and should be in lowercase.
        Columns of the DataFrame should be following:
        ['id', 'value_idx', 'infoname']
        The 'value_idx' column contains 0-origin indice of INFO values.
        This is important when one SV record has multiple values of an INFO (eg. cipos). 
        Main key is the combination of ('id', 'value_idx'), and 'id' is the foreign key coming from df_svpos table.
    df_formats: DataFrame
        DataFrame containing FORMAT information of the vcf file.
        Columns of the DataFrame should be following:
        ['id', 'sample', 'format', 'value_idx', 'value']
        Main key is the combination of ('id', 'sample', 'format').
        The ('id', 'sample', 'format') are the foreign key coming from 
        (df_svpos, samples_meta, format_meta) table, respectively.
    """
    _internal_attrs = [
        "_df_svpos",
        "_df_filters",
        "_odict_df_info",
        "_df_formats",
        "_ls_infokeys",
        "_odict_df_headers",
        "_metadata",
        "_odict_alltables",
        "_repr_config",
        "_sig_criteria"
    ]
    _internal_attrs_set = set(_internal_attrs)
    _repr_column_names = [
        "id",
        "be1",
        "be2",
        "strand",
        "qual",
        "svtype",
    ]
    _repr_column_names_set = set(_repr_column_names)
    def __init__(self, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers = {}, metadata = None):
        if not isinstance(odict_df_info, OrderedDict):
            raise TypeError('the type of the argument "odict_df_info" should be collections.OrderedDict')
        if not isinstance(odict_df_headers, OrderedDict):
            raise TypeError('the type of the argument "odict_df_headers" should be collections.OrderedDict')
        df_svpos['alt'] = df_svpos['alt'].astype(str)
        self._df_svpos = df_svpos
        self._df_filters = df_filters
        self._odict_df_info = odict_df_info
        self._df_formats = df_formats
        self._odict_df_headers = odict_df_headers
        self._metadata = metadata
        self._ls_infokeys = [ x.lower() for x in odict_df_headers['infos_meta']['id'].tolist()]
        ls_keys = ['positions', 'filters'] + self._ls_infokeys + ['formats'] + \
        list(odict_df_headers.keys())
        ls_values = [df_svpos, df_filters] + list(odict_df_info.values()) + [df_formats] + list(odict_df_headers.values())
        # self._odict_alltables is a {tablename: table} dictionary
        self._odict_alltables = OrderedDict([(k, v) for k, v in zip(ls_keys, ls_values)])
        self._repr_config = {
            'info': None,
        }
    
    @property
    def contigs(self) -> List[str]:
        """
        Return a list of contigs (chromosomes) listed in the header of the VCF file.
        """
        df_contigs_meta = self.get_table('contigs_meta')
        arr_contigs = df_contigs_meta['id'].unique()
        return list(arr_contigs)

    def __repr__(self):
        return super().__repr__() 
    def __str__(self):
        return super().__repr__() 

    def replace_svid(self, to_replace, value):
        """
        replace_svid(to_replace, value)
        Renamed specified SV ID.

        Parameters
        ----------
        to_replace: int or str or List[int or str]
            SV ID which are replaced.
        value: int or str or List[int or str]
            Values of new SV ID.
        """
        if not isinstance(to_replace, list):
            to_replace = [to_replace]
        if not isinstance(value, list):
            value = [value]
        if len(to_replace) != len(value):
            raise ValueError('Two arguments should be the same length. {} vs {}'.format(len(to_replace), len(value)))
        
        set_table_list = set(self.table_list)
        set_table_list_header = set(self._odict_df_headers.keys())
        set_table_list_without_header = set_table_list - set_table_list_header
        for rep, val in zip(to_replace, value):
            for table_name in set_table_list_without_header:
                df_target = self._odict_alltables[table_name]
                df_target.loc[df_target['id'] == rep, 'id'] = val
                self._odict_alltables[table_name] = df_target
                if table_name in self._ls_infokeys:
                    self._odict_df_info[table_name.upper()] = df_target


    
    def add_info_table(self, table_name, table, number, type_, description, source=None, version=None):
        self._ls_infokeys += [table_name]
        self._odict_df_info[table_name.upper()] = table
        self._odict_alltables[table_name] = table
        df_meta = self.get_table('infos_meta')
        df_replace = df_meta.append({'id': table_name.upper(), 'number': number, 'type': type_, 'description': description, 'source': source, 'version': version},
                                    ignore_index=True)
        self._odict_df_headers['infos_meta'] = df_replace
        self._odict_alltables['infos_meta'] = df_replace # not beautiful code...
    
    def remove_info_table(self, table_name):
        del self._odict_df_info[table_name.upper()]
        del self._odict_alltables[table_name]
        df_replace = self.get_table('infos_meta')
        df_replace = df_replace.loc[df_replace['id'] != table_name.upper()]
        self._odict_df_headers['infos_meta'] = df_replace
        self._odict_alltables['infos_meta'] = df_replace
        self._ls_infokeys.remove(table_name)
        
    
    def drop_by_id(self, svid):
        """
        drop_by_id(svid)
        Remove SV records specified in "svid" argument.
        
        Paramters
        ---------
        svid: int or str or List[int or str]
            ID of SV record to be removed.
        inplace: bool, default False
            If False, return a copy. Otherwise, dropped SV record of the self and return None.
        
        Returns
        --------
        Vcf
            Return a removed Vcf instance.
        """
        if not isinstance(svid, list):
            svid = [svid]
        set_svid = set(svid)

        set_svid_all = set(self.ids)

        set_svid_preserved = set_svid_all - set_svid

        vcf_removed = self.filter_by_id(set_svid_preserved)

        return vcf_removed
    
    def copy(self):
        """
        copy()
        Return copy of the instance.
        """
        df_svpos = self.get_table('positions')
        df_filters = self.get_table('filters')
        odict_df_infos = OrderedDict([(k, self.get_table(k.lower())) for k, v in self._odict_df_info.items()])
        df_formats = self.get_table('formats')
        odict_df_headers = OrderedDict([(k, self.get_table(k)) for k,v in self._odict_df_headers.items()])
        metadata = self._metadata
        return Vcf(df_svpos, df_filters, odict_df_infos, df_formats, odict_df_headers, metadata)

    def to_vcf_like(self) -> pd.DataFrame:
        """
        to_vcf_like()
        Return a vcf-formatted DataFrame. Header information will not be reflected.
        """
        df_base_before_position_modification = self.get_table('positions')[['chrom1', 'pos1', 'id', 'ref', 'alt', 'qual', 'svtype', 'strand1']]
        def _modify_positions(x):
            svtype = x.name
            if svtype == 'DUP':
                x['pos1'] = x['pos1'] - 1
                return x
            elif svtype == 'INV':
                # if strand1 == '-', subtract 1 from pos1, otherwise subtract 0.
                arr_strand1 = x['strand1'].values 
                arr_num_of_subtraction = np.where(arr_strand1 == '+', 0, 1)
                x['pos1'] = x['pos1'] - arr_num_of_subtraction
                return x
            else:
                return x
                
        df_base = df_base_before_position_modification.groupby('svtype').apply(_modify_positions)
        df_base = df_base[['chrom1', 'pos1', 'id', 'ref', 'alt', 'qual']]
        df_base['qual'] = df_base['qual'].fillna('.')
        ser_id = df_base['id']

        ser_filter = pd.Series(["" for i in range(len(ser_id))], index=ser_id) 
        df_filter = self.get_table('filters')
        def _create_filter_field(x):
            out = ';'.join(x['filter'])
            return out
        ser_filter = df_filter.groupby('id').apply(_create_filter_field)
        df_filter = ser_filter.reset_index(name='filter')

        ser_vcfinfo = pd.Series(["" for i in range(len(ser_id))], index=ser_id)
        def _create_info_field(x, info):
            if x.iloc[1] > 0:
                return ',' + str(x.iloc[2])
            if type(x.iloc[2]) == bool:
                return ';' + info.upper()
            return ';' + info.upper() + '=' + str(x.iloc[2])
        for info in self._ls_infokeys:
            df_info = self.get_table(info)
            ser_be_appended = df_info.apply(_create_info_field, axis=1, **{'info':info})
            if ser_be_appended.empty:
                continue
            df_info_appended = df_info.copy()
            df_info_appended['info'] = ser_be_appended
            df_vcfinfo = df_info_appended.pivot(index='id', columns='value_idx', values='info')
            df_vcfinfo = df_vcfinfo.fillna('')
            ser_vcfinfo_to_append = df_vcfinfo.apply(''.join, axis=1)
            ser_vcfinfo.loc[ser_vcfinfo_to_append.index] = ser_vcfinfo.loc[ser_vcfinfo_to_append.index] + ser_vcfinfo_to_append
        ser_vcfinfo.replace("^;", "", regex=True, inplace=True)
        df_infofield = ser_vcfinfo.reset_index(name='info')

        df_format = self.get_table('formats')
        ls_samples = self.get_table('samples_meta')['id']

        def _create_format_field(x):
            arr_format_, ind_format_ = np.unique(x['format'], return_index=True)
            arr_format_ = arr_format_[np.argsort(ind_format_)]
            format_ = ':'.join(arr_format_)
            ls_sample_format_ = []
            for sample in ls_samples:
                ls_sample_values_ = []
                for a_format_ in arr_format_:
                    mask = (x['sample'] == sample) & (x['format'] == a_format_)
                    ls_sample_values_.append(','.join(x.loc[mask]['value'].astype(str)))
                ls_sample_format_.append(':'.join(ls_sample_values_))

            out_idx = ['format'] + list(ls_samples)
            return pd.Series([format_]+ls_sample_format_, index=out_idx)

        df_formatfield = df_format.groupby('id').apply(_create_format_field).reset_index()

        df_out = pd.merge(df_base, df_filter)
        df_out = pd.merge(df_out, df_infofield)
        df_out = pd.merge(df_out, df_formatfield)
        return df_out

    def to_vcf(self, path_or_buf = None, onlyinfo=False) -> str:
        """
        to_vcf()
        Return a vcf-formatted String. Header information will not be reflected.
        return csv file as str class.

        Parameters
        ----------
        path_or_buf: str, optional
            File path to save the VCF file.
        onlyinfo: bool
            if you only want "info", set this option to True
        
        Returns
        -------
        str
            return vcf file as a string.
        """

        def get_metadata():
            metadata = self._metadata
            out = ''
            if metadata is None:
                return out
            for key, value in metadata.items():
                if not isinstance(value, list):
                    value = [value]
                value = [str(s) for s in value]
                out += '##' + str(key) + '=' + ','.join(value) + '\n'
            return out
        
        def get_contig():
            df_contig = self.get_table('contigs_meta')
            if df_contig.empty:
                return ''
            ser_contig = '##contig=<ID=' + df_contig['id'].astype(str) + ',length=' + df_contig['length'].astype(str) + '>'
            out = '\n'.join(ser_contig)
            out += '\n'
            return out

        def get_info():
            str_info = ""
            for row in self.get_table("infos_meta").itertuples():
                if (row.number == None):
                    str_num = "."
                elif (row.number == -1):
                    str_num = "A"
                else:
                    str_num = str(row.number)
                str_info += "##INFO=<ID={},Number={},Type={},Description=\"{}\">".format(row.id, str_num, row.type,row.description)
                str_info += "\n"
            return str_info
        
        def get_format():
            df_format = self.get_table('formats_meta')
            if df_format.empty:
                return ''
            df_format['number'] = df_format['number'].fillna('.')
            ser_out = '##FORMAT=<ID=' + df_format['id'].astype(str) + ',Number=' + df_format['number'].astype(str) + \
            ',Type=' + df_format['type'].astype(str) + ',Description="' + df_format['description'].astype(str) + '">'
            out = '\n'.join(ser_out)
            out += '\n'
            return out
        
        def get_filter():
            df_filter = self.get_table('filters_meta')
            if df_filter.empty:
                return ''
            ser_out = '##FILTER=<ID=' + df_filter['id'].astype(str) + ',Description="' + df_filter['description'].astype(str) + '">'
            out = '\n'.join(ser_out)
            out += '\n'
            return out
        
        def get_alt():
            df_alt = self.get_table('alts_meta')
            if df_alt.empty:
                return ''
            ser_out = '##ALT=<ID=' + df_alt['id'].astype(str) + ',Description="' + df_alt['description'].astype(str) + '">'
            out = '\n'.join(ser_out)
            out += '\n'
            return out

        str_metadata = get_metadata()
        str_contig = get_contig()
        str_info = get_info()
        str_format = get_format()
        str_filter = get_filter()
        str_alt = get_alt()
        df_vcflike = self.to_vcf_like()
        str_table = df_vcflike.to_csv(sep='\t', header=False, index=False)
        ls_header = df_vcflike.columns.tolist()
        ls_header[0:9] = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
        str_header = "\t".join(ls_header)
        str_header += "\n"

        ls_vcf_data = [str_metadata, str_contig, str_info, str_format, str_filter, str_alt, str_header, str_table]

        #print(os.getcwd())

        if (onlyinfo):
            ret = str_info
        else:
            ret = "".join(ls_vcf_data)

        if (path_or_buf is not None):
            f = open(path_or_buf, 'w')
            f.write(ret)
            f.close()

        #return ret

    def to_bedpe_like(
        self,
        custom_infonames: Iterable[str] = [],
        add_filters: bool = False,
        add_formats: bool = False, 
        confidence_intervals: bool = False,
        unique_events: bool = False
    ) -> pd.DataFrame:
        """
        to_bedpe_like(custom_infonames=[], add_filters, add_formats, confidence_intervals: bool=False)
        Return a DataFrame in bedpe-like format.
        When specified, you can add INFOs, FILTERs, and FORMATs as additional columns.

        Parameters
        ---------------
        custom_infonames: list-like[str]
            The table names of INFOs to append.
        add_filters: bool, default False
            sth
        add_formats: bool, default False
            sth
        confidence_intervals: bool, default False
            Whether or not to consider confidence intervals of the breakpoints.  
            If True, confidence intervals for each breakpoint are represented by [start1, end1) and [start2, end2), respectively.
            Otherwise, breakpoints are represented by a single-nucleotide resolution.
        unique_events: bool, default False
        
        Returns 
        ---------------
        DataFrame
            A Dataframe in bedpe-like format.
            The columns include at least the following:  
            ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
             'name', 'score', 'strand1', 'strand2']
        """
        df_out = super().to_bedpe_like(confidence_intervals=confidence_intervals)
        if len(custom_infonames) != 0:
            df_out = self.append_infos(df_out, custom_infonames, left_on='name')
        if add_filters:
            df_out = self.append_filters(df_out, left_on='name')
        if add_formats:
            df_out = self.append_formats(df_out, left_on='name')
        if unique_events:
            # deprecated for now
            set_unique_ids = self._get_unique_events_ids()
            df_out.set_index('name', inplace=True)
            df_out = df_out.loc[set_unique_ids]
            df_out.reset_index(inplace=True)
        return df_out

    def append_infos(self, base_df,
        ls_tablenames,
        left_on: str = 'id',
        auto_fillna: bool = True) -> pd.DataFrame:
        """
        append_infos(base_df, ls_tablenames, left_on='id', auto_fillna=True)
        Append INFO tables to the right of the base_df, based on the SV id columns.
        If the name of the SV id column in base_df is not 'id', specify column name into left_on argument. 

        Parameters
        ---------------
        base_df: DataFrame
            The DataFrame to which the INFO tables are appended.
        ls_tablenames: list-like
            The list of INFO table names to be appended.
        left_on: str, default 'id'
            The name of SV id column of base_df
        auto_fillna: bool, default True
            If True, use the header information to handle missing values
            after merging DataFrames.
        
        Returns
        ---------------
        DataFrame
            A DataFrame which the INFO tables are added.
        
        """
        df = base_df.copy()
        df_infometa = self.get_table('infos_meta')
        for tablename in ls_tablenames:
            df_to_append_pre = self.get_table(tablename)
            df_to_append_pre['new_column_names'] = tablename + '_' + df_to_append_pre['value_idx'].astype(str)
            df_to_append = df_to_append_pre.pivot(index='id', columns='new_column_names', values=tablename)
            df = pd.merge(df, df_to_append, how='left', left_on=left_on, right_index=True)
            info_dtype = df_infometa.loc[df_infometa['id']==tablename.upper(), 'type'].iloc[0]
            len_info = df_to_append.shape[1]
            ls_ind_fancy = [tablename + '_' + str(i) for i in range(len_info)]
            if info_dtype == 'Integer':
                df[ls_ind_fancy] = df[ls_ind_fancy].fillna(0).astype(int)
            elif info_dtype == 'Flag':
                df[ls_ind_fancy] = df[ls_ind_fancy].fillna(False)
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

        if sq0 in self._ls_infokeys:
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
                set_out = self._filter_infos_flag(sq0, exclude=exclude) # defined in Bedpe
                return set_out

            if len(sq) == 3:
                set_out = self._filter_infos(sq[0], 0, sq[1], sq[2]) # defined in Bedpe
            else:
                set_out = self._filter_infos(*sq) # defined in Bedpe
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

        # is_locus?
        if sq0 in ['be1', 'be2', 'pos1', 'pos2']:
            split_locus = sq[1].split(':')
            chrom = split_locus[0]
            if chrom.startswith('!'):
                exclude_flag = True
                chrom = chrom[1:]
            else:
                exclude_flag = False
            if chrom not in self.contigs:
                raise ContigNotFoundError(chrom)
            if len(split_locus) == 1:
                st = None
                en = None
            elif len(split_locus) == 2:
                split_locus_coord = split_locus[1].split('-')
                if len(split_locus_coord) == 1:
                    st = int(split_locus_coord[0])
                    en = int(split_locus_coord[0]) + 1
                elif len(split_locus_coord) == 2:
                    st = split_locus_coord[0]
                    en = split_locus_coord[1]
                    if st == '':
                        st = None
                    else:
                        st = int(st)
                    if en == '':
                        en = None
                    else:
                        en = int(en)
            
            if sq0 in ['be1', 'pos1']:
                pos_num = 1
            elif sq0 in ['be2', 'pos2']:
                pos_num = 2

            args = [pos_num, chrom, st, en]
            
            if exclude_flag:
                return self._filter_by_positions_exclude(*args)
            return self._filter_by_positions(*args)

    def filter(self, ls_query, query_logic='and'):
        """
        filter(ls_query, query_logic)
        Filter Vcf object by the list of queries.
        Return object is also an instance of the Vcf object
        """
        ### != operation is dangerous
        if isinstance(ls_query, str):
            ls_query = [ls_query]
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
        """
        filter_by_id(arrlike_id)
        Filter Vcf object according to the list of SV ids.
        Return object is also an instance of the Vcf object

        Parameters
        ---------------
        arrlike_id: list-like
            SV ids which you would like to keep.
        
        Returns
        ---------------
        Vcf
            A Vcf object with the SV id specified in the arrlike_id argument.
            All records associated with SV ids that are not in the arrlike_id will be discarded.
        
        """
        out_svpos = self._filter_by_id('positions', arrlike_id)
        out_filters = self._filter_by_id('filters', arrlike_id)
        out_odict_df_info = OrderedDict([(k.upper(), self._filter_by_id(k, arrlike_id)) for k in self._ls_infokeys])
        out_formats = self._filter_by_id('formats', arrlike_id)
        out_odict_df_headers = self._odict_df_headers.copy()
        out_metadata = self._metadata
        return Vcf(out_svpos, out_filters, out_odict_df_info, out_formats, out_odict_df_headers, out_metadata)

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

    def _filter_formats(self, sample, item, item_idx=0, operator=None, threshold=None):
        df = self.get_table('formats')
        target_q = (df['sample'] == sample) & (df['format'] == item) & (df['value_idx'] == item_idx)
        df_target = df.loc[target_q]
        e = "df_target.loc[df_target['value'] {0} threshold]['id']".format(operator)
        return set(eval(e))
    
    def _filter_header(self, tablename):
        pass

    def annotate_bed(self, bed: Bed, annotation: str, suffix=['left', 'right'], description=None):
        df_svpos = self.get_table('positions')
        ls_left = []
        ls_right = []
        for idx, row in df_svpos.iterrows():
            svid = row['id']
            chrom1 = row['chrom1']
            pos1 = row['pos1']
            chrom2 = row['chrom2']
            pos2 = row['pos2']

            df_bp1 = bed.query(chrom1, pos1)
            df_bp2 = bed.query(chrom2, pos2)

            if not df_bp1.empty:
                ls_left.append([svid, 0, True])
            if not df_bp2.empty:
                ls_right.append([svid, 0, True])
        left_name = annotation + suffix[0]
        right_name = annotation + suffix[1]
        df_left = pd.DataFrame(ls_left, columns=('id', 'value_idx', left_name))
        df_right = pd.DataFrame(ls_right, columns=('id', 'value_idx', right_name))
        self.add_info_table(left_name, df_left, 0, type_="Flag", description=description)
        self.add_info_table(right_name, df_right, 0, type_="Flag", description=description)
    
    def breakend2breakpoint(self):
        """
        breakend2breakpoint()
        Transforms SV records whose svtype are BND, infer the SV type, and return a breakpoint-based Vcf object. 

        Returns
        --------
        Vcf
            SV records with svtype being BND were integrated into breakpoints, and svtype will be overwritten.
        """
        out = self.copy()
        if out._metadata['variantcaller'] == 'lumpy':
            ls_secondary = self.get_table('secondary')['id'].tolist()
            out = out.drop_by_id(ls_secondary)
            out.remove_info_table('secondary')
        df_svpos = out.get_table('positions')
        df_svtype = out.get_table('svtype')

        if out._metadata['variantcaller'] == 'delly':
            df_svpos.loc[df_svpos['svtype'] == 'BND', 'svtype'] = 'TRA'
            df_svtype.loc[df_svtype['svtype'] == 'BND', 'svtype'] = 'TRA'
            out._odict_alltables['positions'] = df_svpos
            out._odict_alltables['svtype'] = df_svtype
            out._odict_df_info['SVTYPE'] = df_svtype
            return out

        df_mateid = out.get_table('mateid')
        df_bnd = df_svpos[df_svpos['svtype'] == 'BND']
        ls_info_breakend_id = []
        breakpoint_id_num = 0
        if df_bnd.empty:
            return self
        arr_skip = np.array([])
        for idx, row in df_bnd.iterrows():
            svid = row['id']
            ser_mateid = df_mateid.loc[df_mateid['id'] == svid, 'mateid']
            if ser_mateid.empty:
                mateid = None
            else:
                mateid = ser_mateid.item()
            if np.isin(svid, arr_skip):
                continue
            if mateid is None:
                svtype = 'BND'
            elif row['chrom1'] != row['chrom2']:
                svtype = 'TRA'
            elif row['strand1'] == row['strand2']:
                svtype = 'INV'
            elif get_inslen_and_insseq_from_alt(row['alt'])[0] > abs(row['pos1'] - row['pos2']) * 0.5:
                svtype = 'INS'
            elif (row['pos1'] < row['pos2']) & (row['strand1'] == '-') & (row['strand2'] == '+'):
                svtype = 'DUP'
            elif (row['pos1'] > row['pos2']) & (row['strand1'] == '+') & (row['strand2'] == '-'):
                svtype = 'DUP'
            else:
                svtype = 'DEL'
            
            arr_skip = np.append(arr_skip, mateid)

            breakpoint_id = 'viola_breakpoint:' + str(breakpoint_id_num)
            if mateid is not None:
                ls_info_breakend_id += [[breakpoint_id, 0, svid], [breakpoint_id, 1, mateid]]
            else:
                ls_info_breakend_id += [[breakpoint_id, 0, svid]]


            df_svpos.loc[df_svpos['id'] == svid, ['id', 'svtype']] = [breakpoint_id, svtype]
            df_svtype.loc[df_svtype['id'] == svid, ['id', 'svtype']] = [breakpoint_id, svtype]

            out.replace_svid(svid, breakpoint_id)
            
            breakpoint_id_num += 1
        
        out._odict_alltables['positions'] = df_svpos
        out._odict_alltables['svtype'] = df_svtype
        out._odict_df_info['SVTYPE'] = df_svtype
        out.remove_info_table('mateid')

        df_info_breakend_id = pd.DataFrame(ls_info_breakend_id, columns=('id', 'value_idx', 'orgbeid'))
        out.add_info_table('orgbeid', df_info_breakend_id, type_='String', number=2, description='Breakend ID which were in original VCF file.', source='Python package, Viola-SV.')

        if out._metadata['variantcaller'] != 'lumpy':
            out = out.drop_by_id(list(arr_skip))

        return out
            
    def _get_unique_events_ids(self) -> Set[IntOrStr]:
        """
        now yields errors!!!
        ATOMAWASHI!!!
        """
        if 'mateid' not in self.table_list:
            print("Can't find mateid table")
            return
        df = self.get_table('mateid')
        df2 = df.reset_index().set_index('mateid_0').loc[df['id']]
        arr_mask = df2['index'].values > np.arange(df2.shape[0])
        set_to_subtract = set(df2.loc[arr_mask]['id'])
        set_all_ids = self.get_ids()
        set_result_ids = set_all_ids - set_to_subtract
        return set_result_ids

    def remove_duplicated_records(self):
        if 'event' not in self._ls_infokeys:
            return
        print(self.get_table('event'))
    
            