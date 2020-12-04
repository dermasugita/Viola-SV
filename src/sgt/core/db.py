import numpy as np
import pandas as pd
from typing import (
    List,
    Dict,
    Set,
    Iterable,
)
import sgt
from sgt.core.indexing import Indexer
from sgt.core.bed import Bed
from sgt._typing import (
    IntOrStr,
    StrOrIterableStr,
)
from sgt._exceptions import (
    TableNotFoundError,
)

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
    dict_df_info: dict[str, DataFrame]
        Dictionary of DataFrames which contain additional information on SV record (equivalent to INFO field of vcf).
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
        "_dict_df_info",
        "_ls_infokeys",
        "_dict_alltables",
        "_repr_config",
        "_sig_criteria"
    ]
    _internal_attrs_set = set(_internal_attrs)
    _repr_column_names = [
        "id",
        "bp1",
        "bp2",
        "strand",
        "qual",
        "svtype",
    ]
    _repr_column_names_set = set(_repr_column_names)

    def __init__(self, df_svpos: pd.DataFrame, dict_df_info: Dict[str, pd.DataFrame]):
        self._df_svpos = df_svpos
        self._dict_df_info = dict_df_info
        self._ls_infokeys = [x.lower() for x in dict_df_info.keys()]
        ls_keys = ['positions'] + self._ls_infokeys
        ls_values = [df_svpos] + list(dict_df_info.values())
        self._dict_alltables = {k: v for k, v in zip(ls_keys, ls_values)}
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
        return list(self._dict_alltables.keys())
    
    @property
    def ids(self):
        """
        Return all SV ids as list.
        """
        return list(self.get_ids())
    
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
        self._dict_alltables[table_name] = df
    
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
        ser_bp1 = df_svpos['chrom1'].astype(str) + ':' + df_svpos['pos1'].astype(str)
        ser_bp2 = df_svpos['chrom2'].astype(str) + ':' + df_svpos['pos2'].astype(str)
        ser_strand = df_svpos['strand1'] + df_svpos['strand2']
        ser_qual = df_svpos['qual']
        ser_svtype = df_svpos['svtype']
        ls_ser = [ser_id, ser_bp1, ser_bp2, ser_strand, ser_qual, ser_svtype]
        ls_key = ['id', 'bp1', 'bp2', 'strand', 'qual', 'svtype']
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
        table = self._dict_alltables[table_name]
        return table.copy()

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

    def filter(self,
        ls_query: StrOrIterableStr,
        query_logic: str = 'and'):
        """
        filter(ls_query, query_logic)
        Filter SgtSimple object by the list of queries.
        Return object is also an instance of the SgtSimple object
        """
        ### != operation is dangerous
        set_result = self._filter(ls_query, query_logic)
        out = self.filter_by_id(set_result)
        return out


    def _filter_by_id(self, tablename, arrlike_id):
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
            SV ids which you would like to keep.
        
        Returns
        ---------------
        Bedpe
            A Bedpe object with the SV id specified in the arrlike_id argument.
            All records associated with SV ids that are not in the arrlike_id will be discarded.
        
        """
        out_svpos = self._filter_by_id('positions', arrlike_id)
        out_dict_df_info = {k: self._filter_by_id(k, arrlike_id) for k in self._ls_infokeys}
        return Bedpe(out_svpos, out_dict_df_info)

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
    
    def annotate_bed(self, bed: Bed, annotation: str, suffix=['left', 'right']):
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
        self.add_info_table(left_name, df_left)
        self.add_info_table(right_name, df_right)

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
            set_ids = func(obj)
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
        ---------------
        positions_num:int
            1 for the first breakend and 2 for the other 
        chrom:str
            "chr1","chr2",...,"chrX","chrY" <- Vcf
            1,2,... <- Bedpe
        pos_min:int
        pos_sup:int

        Returns
        ---------------
        set
            set of ids which satisfies the argument
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
        Return set of ids except which are specified by the argument
        ---------------
        ex_positions_num:int
            1 or 2
        ex_chrom:str
            "chr1","chr2",...,<- Vcf
            1,2,... <- Bedpe
        ex_pos_min:int
        ex_pos_max:int

        Returns
        ---------------
        set
            set of ids except which satisfies the argument
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
    dict_df_info: dict[str, DataFrame]
        Dictionary of DataFrames which contain additional information on SV record (equivalent to INFO field of vcf).
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
    def __init__(self, df_svpos, df_filters, dict_df_info, df_formats, dict_df_headers = {}):
        self._df_svpos = df_svpos
        self._df_filters = df_filters
        self._dict_df_info = dict_df_info
        self._df_formats = df_formats
        self._dict_df_headers = dict_df_headers
        self._ls_infokeys = [ x.lower() for x in dict_df_headers['infos_meta']['id'].tolist()]
        ls_keys = ['positions', 'filters'] + self._ls_infokeys + ['formats'] + \
        list(dict_df_headers.keys())
        ls_values = [df_svpos, df_filters] + list(dict_df_info.values()) + [df_formats] + list(dict_df_headers.values())
        # self._dict_alltables is a {tablename: table} dictionary
        self._dict_alltables = {k: v for k, v in zip(ls_keys, ls_values)}
        self._repr_config = {
            'info': None,
        }

    def __repr__(self):
        return super().__repr__() 
    def __str__(self):
        return super().__repr__() 
    
    def create_info_table(self, table_name, table, number, type_, description, source=None, version=None):
        self._ls_infokeys += [table_name]
        self._dict_alltables[table_name] = table
        df_meta = self.get_table('infos_meta')
        df_replace = df_meta.append({'id': table_name.upper(), 'number': number, 'type': type_, 'description': description, 'source': source, 'version': version},
                                    ignore_index=True)
        self._dict_alltables['infos_meta'] = df_replace # not beautiful code...


    def to_vcf_like(self) -> pd.DataFrame:
        """
        to_vcf_like()
        Return a vcf-formatted DataFrame. Header information will not be reflected.
        """
        df_base = self.get_table('positions')[['chrom1', 'pos1', 'id', 'ref', 'alt', 'qual']]
        print(self._ls_infokeys)
        ser_id = df_base['id']
        ser_vcfinfo = pd.Series(["" for i in range(len(ser_id))], index=ser_id)
        def _create_info_field(x, info):
            if x.iloc[0] > 0:
                return ',' + str(x.iloc[1])
            if type(x.iloc[1]) == bool:
                return ';' + info.upper()
            return ';' + info.upper() + '=' + str(x.iloc[1])
        for info in self._ls_infokeys:
            df_info = self.get_table(info).set_index('id')
            ser_be_appended = df_info.apply(_create_info_field, axis=1, **{'info':info})
            ser_vcfinfo.loc[df_info.index] = ser_vcfinfo.loc[df_info.index] + ser_be_appended
        ser_vcfinfo.replace("^;", "", regex=True, inplace=True)
        df_infofield = ser_vcfinfo.reset_index(name='info')

        df_format = self.get_table('formats')
        ls_samples = self.get_table('samples_meta')['id']

        def _create_format_field(x):
            arr_format_ = np.unique(x['format'])
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

        df_out = pd.merge(df_base, df_infofield)
        df_out = pd.merge(df_out, df_formatfield)
        return df_out

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
        out_dict_df_info = {k: self._filter_by_id(k, arrlike_id) for k in self._ls_infokeys}
        out_formats = self._filter_by_id('formats', arrlike_id)
        out_dict_df_headers = self._dict_df_headers.copy()
        return Vcf(out_svpos, out_filters, out_dict_df_info, out_formats, out_dict_df_headers)

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
        self.create_info_table(left_name, df_left, 0, type_="Flag", description=description)
        self.create_info_table(right_name, df_right, 0, type_="Flag", description=description)
    
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

    def get_unique_events(self):
        set_result_ids = self._get_unique_events_ids()
        return self.filter_by_id(set_result_ids)
    
            

         