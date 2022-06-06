import os,sys
import numpy as np
import pandas as pd
import re
import pkgutil
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
    SVIDNotFoundError,
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
    patient_name
        Patient name.
    
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
        "_sig_criteria",
        "_patient_name"
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

    def __init__(self, df_svpos: pd.DataFrame, odict_df_info: 'OrderedDict[str, pd.DataFrame]', patient_name=None):
        if not isinstance(odict_df_info, OrderedDict):
            raise TypeError('the type of the argument "odict_df_info" should be collections.OrderedDict')
        self._df_svpos = df_svpos
        self._odict_df_info = odict_df_info
        self._ls_infokeys = [x.lower() for x in odict_df_info.keys()]
        self._patient_name = patient_name
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
    def patient_name(self):
        """
        Return the name of the patient.
        """
        return self._patient_name
    
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
    
    def copy(self):
        """
        copy()
        Return copy of the instance
        """
        df_svpos = self.get_table('positions')
        odict_df_infos = OrderedDict([(k, self.get_table(k.lower())) for k, v in self._odict_df_info.items()])
        patient_name = self.patient_name
        return Bedpe(df_svpos, odict_df_infos, patient_name)

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
        if table_name in self._ls_infokeys:
            self.remove_info_table(table_name)
        self._ls_infokeys += [table_name]
        self._odict_alltables[table_name] = df
        self._odict_df_info[table_name] = df

    def remove_info_table(self, table_name: str):
        """
        remove_info_table(table_name)
        Remove an INFO table from self.

        Parameters
        -----------
        table_name: str
            The name of the INFO table to be removed.
        """
        del self._odict_df_info[table_name]
        del self._odict_alltables[table_name]
        self._ls_infokeys.remove(table_name)
    
    def set_value_for_info_by_id(self, table_name, sv_id, value_idx=0, value=None):
        """
        set_value_for_info_by_id(table_name, sv_id, value_idx, value)
        Set value to the specified info table by sv_id. The value will be overwrited if it already exists.

        Parameters
        -------------
        table_name: str
            Name of the INFO table.
        sv_id: str or int
            Target SV ID
        value_idx: int, default 0
            0-origin index. This argument should be 0 in most cases unless multiple values are required such as CIPOS and CIEND.
        value: int or str
            INFO value to be set. 
        """
        if table_name not in self._ls_infokeys:
            raise InfoNotFoundError(table_name)
        if sv_id not in self.ids:
            raise SVIDNotFoundError(sv_id)
        df = self.get_table(table_name) 
        df.set_index('id' ,inplace=True)
        if df.empty:
            # Python 3.6 and 3.7 do not infer dtypes of new values when the df is empty.
            df = df.astype({'value_idx': int, table_name: type(value)})
        df.loc[sv_id] = [value_idx, value]
        df.reset_index(inplace=True)
        self.replace_table(table_name, df)

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
        desc_info = 'INFO='
        desc_doc = 'Documentation of Bedpe object ==> '
        doc_link = 'https://dermasugita.github.io/ViolaDocs/docs/html/reference/bedpe.html'
        out = desc_info + str_infokeys + '\n' + desc_doc + doc_link + '\n' + str_df_out
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
                df_svpos['start1'] = df_svpos['pos1'] + df_svpos['cipos_0'] - 1
                df_svpos['end1'] = df_svpos['pos1'] + df_svpos['cipos_1']
                df_svpos['start2'] = df_svpos['pos2'] + df_svpos['ciend_0'] - 1
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
    
    def to_bedpe(self,
        path_or_buf: str, 
        custom_infonames: Iterable[str] = [],
        confidence_intervals: bool = False):
        """
        to_beddpe(path_or_buf, custom_infonames, confidence_intervals)
        Return a BEDPE file.

        Parameters
        ----------
        path_or_buf: str, optional
            File path to save the BEDPE file.
        custom_infonames: list-like[str]
            The table names of INFOs to append.
        confidence_intervals: bool, default False
            Whether or not to consider confidence intervals of the breakpoints.  
            If True, confidence intervals for each breakpoint are represented by [start1, end1) and [start2, end2), respectively.
            Otherwise, breakpoints are represented by a single-nucleotide resolution.
        """
        df_bedpe = self.to_bedpe_like(custom_infonames=custom_infonames, confidence_intervals=confidence_intervals)
        df_bedpe.to_csv(path_or_buf, index=None, sep='\t')

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

        if sq0.lower() in self._ls_infokeys:
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
                elif len(sq) == 1:
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
                if set_query is None:
                    set_query = set()
                set_result = set_result & set_query
        elif query_logic == 'or':
            set_result = set()
            for query in ls_query:
                set_query = self._parse_filter_query(query)
                if set_query is None:
                    set_query = set()
                set_result = set_result | set_query
        else:
            ls_set_query = [self._parse_filter_query(q) for q in ls_query]
            ls_set_query = [set() if q is None else q for q in ls_set_query]
            pattern = re.compile('([^0-9]*)([0-9]+)([^0-9]*)')
            target = r"\1ls_set_query[\2]\3"
            expr = pattern.sub(target, query_logic)
            set_result = eval(expr)
        
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
        df = self.get_table(tablename.lower())
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
        return Bedpe(out_svpos, out_odict_df_info, self.patient_name)

    def _filter_pos_table(self, item, operator, threshold):
        df = self.get_table('positions')
        e = "df.loc[df[item] {0} threshold]['id']".format(operator)
        return set(eval(e))

    def _filter_infos(self, infoname, value_idx=0, operator=None, threshold=None):## returning result ids
        infoname = infoname.lower()
        df = self.get_table(infoname)
        value_idx = int(value_idx)
        df = df.loc[df['value_idx'] == value_idx]
        e = "df.loc[df[infoname] {0} threshold]['id']".format(operator)
        set_out = set(eval(e))
        return set_out

    def _filter_infos_flag(self, infoname, exclude=False):
        df = self.get_table(infoname)
        df = df.loc[df[infoname] == True]
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
        
    def classify_manual_svtype(self, definitions=None, ls_conditions=None, ls_names=None, ls_order=None, return_series=True):
        """
        classify_manual_svtype(definitions, ls_conditions, ls_names, ls_order=None)
        Classify SV records by user-defined criteria. A new INFO table named
        'manual_sv_type' will be created.

        Parameters
        ------------
        definitions: path_or_buf or str, default None
            Path to the file which specifies the definitions of custom SV classification. This argument is disabled when "ls_condition" is not None.
            If "default" is specified, the simple length-based SV classification will be employed.
            If "article" is specified, the same definition file which was used in the Viola publication will be reflected.
            Below is the links to each of definition file you can specify on this method.

            "default" -> https://github.com/dermasugita/Viola-SV/blob/master/examples/demo_sig/resources/definitions/sv_class_default.txt

            "article" -> https://github.com/dermasugita/Viola-SV/blob/master/examples/demo_sig/resources/definitions/sv_class_article.txt
        ls_conditions: List[callable] or List[str], default None
            List of definitions of custom SV classification. The data type of the elements in the list can be callable or SV ID (str).
            callable --> Functions that takes a self and returns a list of SV ID that satisfy the conditions of the SV class to be defined. 
            SV ID --> Lists of SV ID that satisfy the conditions of the SV class to be defined.
            This argument is disabled when "definitions" is not None.
        ls_names: List[str], default None
            List of the names of the custom SV class corresponding to the "ls_conditions". This argument is disabled when "definitions" is not None.
        return_series: bool, default True
            Return counts of each custom SV class as a pd.Series.
        
        Returns
        ---------
        pd.Series or None
        """
        set_ids_current = set(self.ids)
        obj = self
        ls_ids = []
        ls_result_names = []

        if definitions is not None:
            if isinstance(definitions, str):
                if definitions == "default":
                    d = os.path.dirname(sys.modules["viola"].__file__)
                    definitions = os.path.join(d, "data/sv_class_default.txt")
                    ls_conditions, ls_names = self._parse_signature_definition_file(open(definitions, 'r'))
                elif definitions == "article":
                    d = os.path.dirname(sys.modules["viola"].__file__)
                    definitions = os.path.join(d, "data/sv_class_article.txt")
                    ls_conditions, ls_names = self._parse_signature_definition_file(open(definitions, 'r'))
                else:
                    ls_conditions, ls_names = self._parse_signature_definition_file(open(definitions, 'r'))
            else:
                ls_conditions, ls_names = self._parse_signature_definition_file(definitions)

        for cond, name in zip(ls_conditions, ls_names):
            obj = obj.filter_by_id(set_ids_current)
            if callable(cond):
                ids = cond(obj)
            else:
                ids = cond
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
        """
        get_feature_count_as_series(feature, ls_order)
        Return counts of unique values as a pd.Series for the INFO specified in the "feature" argument.

        Parameters
        -----------
        feature: str, default 'manual_sv_type'
            The name of INFO to be counted
        ls_order: List[str], default None
            Order of the index (unique feature values) of the output Series.
        """
        ser_feature_counts = self.get_table(feature)[feature].value_counts()
        if ls_order is not None:
            pd_ind_reindex = pd.Index(ls_order)
            ser_feature_counts = ser_feature_counts.reindex(index=pd_ind_reindex, fill_value=0)
        return ser_feature_counts

    def _parse_signature_definition_file(self, infile):
        ls_query = []
        ls_conditions = []
        ls_names = []
        for line in infile:
            line = line.strip('\n')
            if line.startswith('#'): continue
            if line == '': continue
            if line.startswith('name'):
                sig_name = line[5:]
                ls_names.append(sig_name.strip('"\''))
                continue
            if line.startswith('logic'):
                query_logic = line[6:]
                result = self.filter(ls_query, query_logic=query_logic).ids
                ls_conditions.append(result)
                ls_query=[]
                continue

            pattern = re.compile('^[0-9]+ ')
            ls_query.append(pattern.sub('', line))
        return ls_conditions, ls_names


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
            str1 = [param["str1h"], param["str1w"]]
            str2 = [param["str2h"], param["str2w"]]
        elif mode == "reverse":
            chr1 = [param["chr1h"], param["chr2w"]]
            chr2 = [param["chr2h"], param["chr1w"]]
            str1 = [param["str1h"], param["str2w"]]
            str2 = [param["str2h"], param["str1w"]]

        proposition_chr = (chr1[0] == chr1[1]) and (chr2[0] == chr2[1])
        if str_missing: 
            proposition_str = ((str1[0] == str1[1]) or (str1[0]==".") or (str1[1]==".")) and ((str2[0] == str2[1]) or (str2[0]==".") or (str2[1]=="."))
        else:
            proposition_str = (str1[0] == str1[1]) and (str2[0] == str2[1])

        proposition = proposition_chr and proposition_str
        return proposition

    def _generate_distance_matrix_by_distance(self, multiobject, penalty_length=3e9, str_missing=True):
        """
        _generate_distance_matrix_by_distance(multiobject, penalty_length=3e9)
        Generate distance_matrix by simply calculating the distance between all conbinations of SV records.
        
        Parameters
        -----------
        multiobject: MultiBedpe or TmpVcfForMerge
            Multi-something object which includes all the samples to be merged.
        penalty_length: int or float, default 3e9
            The value which virtually gives constraints to the clustering model not to merge
            the certain pairs of SV records. 
        str_missing: bool, default True
            If True, all the missing strands are considered to be identical to the others.
        
        Returns
        --------
        array
            distance matrix.
        """
        positions_table = multiobject.get_table("positions")
        N = len(positions_table)
        distance_matrix = np.full((N,N), penalty_length)
        columns = ["chrom1", "chrom2", "pos1", "pos2", "strand1", "strand2"]
        for h in range(N):
            for w in range(N):
                param = {}
                for col in columns:
                    key_h = col[:3] + col[-1] + "h"
                    value_h = positions_table.at[positions_table.index[h], col]
                    key_w = col[:3] + col[-1] + "w"
                    value_w = positions_table.at[positions_table.index[w], col]
                    param[key_h] = value_h
                    param[key_w] = value_w

                if self._necessary_condition4merge(param = param, mode = "normal", str_missing = str_missing): 
                    if self._nonoverlap(param = param):
                        distance_matrix[h, w] = penalty_length
                    else:
                        distance_matrix[h, w] = max(np.abs(param["pos1h"] - param["pos1w"]), np.abs(param["pos2h"] - param["pos2w"]))                          

                elif self._necessary_condition4merge(param = param, mode = "reverse", str_missing = str_missing):
                    if self._nonoverlap(param = param):
                        distance_matrix[h, w] = penalty_length
                    else:
                        distance_matrix[h, w] = max(np.abs(param["pos1h"] - param["pos2w"]), np.abs(param["pos2h"] - param["pos1w"]))
        return distance_matrix

    def merge(self, ls_bedpe = [], ls_caller_names = None, threshold = 100, linkage = "complete", str_missing = True):
        """
        merge(ls_bedpe:list, ls_caller_names:list, threshold:float, linkage = "complete", str_missing=True)
        Return a merged bedpe object from mulitple  caller's bedpe objects in ls_bedpe

        Parameters
        ----------
        ls_bedpe:list
            A list of bedpe objects to be merged, which are the same order with ls_caller_names
        ls_caller_names:list
            A list of names of bedpe objects to be merged, which should have self's name as the first element
        threshold:float
            Two SVs whose diference of positions is under this threshold are cosidered to be identical.
        linkage:{‘complete’, ‘average’, ‘single’}, default=’complete’
            The linkage of hierarchical clustering.
            To keep the mutual distance of all SVs in each cluster below the threshold, 
            "complete" is recommended.
        str_missing:boolean, default="True"
            If True, all the missing strands are considered to be identical to the others. 

        Returns
        ----------
        A merged bedpe object
        """
        if self in ls_bedpe:
            pass
        else:
            ls_bedpe = [self] + ls_bedpe

        multibedpe = viola.MultiBedpe(ls_bedpe, ls_caller_names)
        distance_matrix = self._generate_distance_matrix_by_distance(multibedpe, penalty_length=3e9, str_missing=str_missing)
        hcl_clustering_model = AgglomerativeClustering(n_clusters=None, affinity="precomputed", linkage=linkage, distance_threshold=threshold)
        labels = hcl_clustering_model.fit_predict(X = distance_matrix)
        
        positions_table = multibedpe.get_table("positions")
        mergedid_dict = {labels[0]:0}
        ls_mergedid = []
        idx_head = 0
        for label in labels:
            if label in mergedid_dict:
                ls_mergedid.append(mergedid_dict[label])
            else:
                idx_head += 1
                mergedid_dict[label] = idx_head
                ls_mergedid.append(mergedid_dict[label])

        N = len(positions_table)
        value_idx = pd.Series(np.zeros(N, dtype=int))
        df_mergedid = pd.DataFrame({"id":positions_table["id"],"value_idx":value_idx, "mergedid":pd.Series(ls_mergedid)})
        
        originalid = multibedpe.get_table("global_id")["id"]
        df_originalid = pd.DataFrame({"id":positions_table["id"], "value_idx":value_idx, "originalid":originalid})
        
        ############## Edited by Sugita ##################
        df_id = multibedpe.get_table("global_id")
        df_patients = multibedpe.get_table("patients")
        df_id_patients = df_id.merge(df_patients, left_on="patient_id", right_on="id")
        caller = df_id_patients["patients"]
        df_caller = pd.DataFrame({"id":positions_table["id"], "value_idx":value_idx, "caller":caller})
        ############## /Edited by Sugita #################

        df_svpos = multibedpe._df_svpos
        odict_df_info = multibedpe._odict_df_info
    
        merged_bedpe = viola.Bedpe(df_svpos=df_svpos, odict_df_info=odict_df_info)
        merged_bedpe.add_info_table(table_name="mergedid", df=df_mergedid)
        merged_bedpe.add_info_table(table_name="originalid", df=df_originalid)
        merged_bedpe.add_info_table(table_name="caller", df=df_caller)
        
        return merged_bedpe