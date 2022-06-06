import os,sys
import numpy as np
import pandas as pd
import re
from intervaltree import Interval, IntervalTree
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
from viola.core.bedpe import Bedpe
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
    IllegalArgumentError,
    SVIDNotFoundError,
)

from sklearn.cluster import AgglomerativeClustering

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
    patient_name
        Patient name.
    contigs
        List of the contigs (chromosomes)

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
        "_sig_criteria",
        "_patient_name",
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
    def __init__(self, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers = {}, metadata = None, patient_name = None):
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
        self._patient_name = patient_name
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
        desc_doc = 'Documentation of Vcf object ==> '
        doc_link = 'https://dermasugita.github.io/ViolaDocs/docs/html/reference/vcf.html'
        out = desc_info + str_infokeys + '\n' + desc_doc + doc_link + '\n' + str_df_out
        if return_as_dataframe:
            return df_out
        return str(out)

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
        if table_name in self._ls_infokeys:
            self.remove_info_table(table_name)
        self._ls_infokeys += [table_name]
        self._odict_df_info[table_name.upper()] = table
        self._odict_alltables[table_name] = table
        df_meta = self.get_table('infos_meta')
        df_replace = pd.concat([df_meta, pd.DataFrame({'id': [table_name.upper()], 'number': [number], 'type': [type_], 'description': [description], 'source': [source], 'version': [version]})], ignore_index=True)
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
        value: int or str or bool
            INFO value to be set. For boolean INFO, set True or False. 
        """
        if table_name not in self._ls_infokeys:
            raise InfoNotFoundError(table_name)
        if sv_id not in self.ids:
            raise SVIDNotFoundError(sv_id)
        df = self.get_table(table_name) 
        df.set_index('id' ,inplace=True)

        info_dtype = self.infos_meta.set_index('id').at[table_name.upper(), 'type']
        if df.empty:
            # Python 3.6 and 3.7 do not infer dtypes of new values when the df is empty.
            df = df.astype({'value_idx': int, table_name: type(value)})
        if info_dtype == 'Flag' and not value:
            if df.empty:
                pass
            elif sv_id not in df.index:
                pass
            else:
                df.drop(sv_id, inplace=True)
        else:
            df.loc[sv_id] = [value_idx, value]

        df.reset_index(inplace=True)
        self.replace_table(table_name, df)

        
    
    def drop_by_id(self, svid):
        """
        drop_by_id(svid)
        Remove SV records specified in "svid" argument.
        
        Parameters
        -----------
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
        patient_name = self.patient_name
        return Vcf(df_svpos, df_filters, odict_df_infos, df_formats, odict_df_headers, metadata, patient_name)

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
                if self._metadata.get('variantcaller', None) != 'delly':
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
        to_vcf(path_or_buf)
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
            for key, values in metadata.items():
                if not isinstance(values, list):
                    values = [values]
                for value in values:
                    out += '##' + str(key) + '=' + str(value) + '\n'
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
                elif np.isnan(row.number):
                    str_num = "."
                elif (row.number == -1):
                    str_num = "A"
                else:
                    str_num = str(int(row.number))
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
        return df_out

    def to_bedpe(
        self,
        file_or_buf: str,
        custom_infonames: Iterable[str] = [],
        add_filters: bool = False,
        add_formats: bool = False, 
        confidence_intervals: bool = False,
    ):
        """
        to_bedpe_like(file_or_buf, custom_infonames=[], add_filters, add_formats, confidence_intervals: bool=False)
        Return a BEDPE file.

        Parameters
        ---------------
        file_or_buf: str
            File path to save the VCF file.
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
        """
        bedpe = self.to_bedpe_like(custom_infonames=custom_infonames, add_filters=add_filters, add_formats=add_formats, confidence_intervals=confidence_intervals)
        bedpe.to_csv(file_or_buf, index=None, sep='\t')
    
    def as_bedpe(self):
        """
        as_bedpe()
        Convert Vcf object into Bedpe object.

        Returns
        --------
        Bedpe

        Notes
        ------
        This process is lossy because only the "positions" table and INFO tables are inherited by Bedpe class.

        Examples
        ---------
        >>> import viola
        >>> vcf = viola.read_vcf('https://raw.githubusercontent.com/dermasugita/ViolaDocs/main/docs/html/_static/tutorial.manta.vcf', patient_name='patient1')
        >>> bedpe = vcf.as_bedpe()
        >>> bedpe
        INFO=imprecise,svtype,svlen,end,cipos,ciend,cigar,mateid,event,homlen,homseq,svinslen,svinsseq,left_svinsseq,right_svinsseq,contig,bnd_depth,mate_bnd_depth,somatic,somaticscore,junction_somaticscore,inv3,inv5
        Documentation of Bedpe object ==> https://dermasugita.github.io/ViolaDocs/docs/html/reference/bedpe.html
                id              be1              be2 strand  qual svtype
        0    test1    chr1:82550461    chr1:82554226     +-  None    DEL
        1    test2    chr1:22814217    chr1:92581132     --  None    INV
        2    test3    chr1:60567906    chr1:60675941     +-  None    DEL
        3    test4    chr1:69583190    chr1:69590948     +-  None    DEL
        4    test5  chr11:104534877  chr11:104536574     +-  None    DEL
        5  test6_1  chr11:111134697   chr17:26470495     +-  None    BND
        6  test6_2   chr17:26470495  chr11:111134697     -+  None    BND
        """
        df_svpos = self.get_table('positions')
        odict_df_info_view = self._odict_df_info
        odict_df_info = OrderedDict((k, v.copy()) for k, v in odict_df_info_view.items())
        bedpe = Bedpe(df_svpos, odict_df_info)
        return bedpe

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
        """
        append_formats(base_df, left_on='id')
        Append formats to the right of the base_df, based on the SV id columns.
        If the name of the SV id column in base_df is not 'id', specify column name into left_on argument. 

        Parameters
        ---------------
        base_df: DataFrame
            The DataFrame to which the INFO tables are appended.
        left_on: str, default 'id'
            The name of SV id column of base_df
        
        Returns
        ---------------
        DataFrame
            A DataFrame which the formats tables are added.
        """
        df_format = self.get_table('formats')
        df_format['format_id'] = df_format['sample'] + '_' + df_format['format'] + '_' + df_format['value_idx'].astype(str) 
        df_format.drop(['sample', 'format', 'value_idx'], axis=1, inplace=True)
        df_format = df_format.pivot(index='id', columns='format_id', values='value')
        df_out = pd.merge(base_df, df_format, how='left', left_on=left_on, right_index=True)
        return df_out

    def append_filters(self, base_df, left_on='id'):
        """
        append_filters(base_df, left_on='id')
        Append filters to the right of the base_df, based on the SV id columns.
        If the name of the SV id column in base_df is not 'id', specify column name into left_on argument. 

        Parameters
        ---------------
        base_df: DataFrame
            The DataFrame to which the INFO tables are appended.
        left_on: str, default 'id'
            The name of SV id column of base_df
        
        Returns
        ---------------
        DataFrame
            A DataFrame which the filters tables are added.
        
        """
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

        if sq0.lower() in self._ls_infokeys:
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
        out_patient_name = self.patient_name
        return Vcf(out_svpos, out_filters, out_odict_df_info, out_formats, out_odict_df_headers, out_metadata, out_patient_name)

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
        Converts a Vcf object into a breakpoint-based Vcf object by integrating the paired breakends (BND) and infering their SVTYPE. 

        Returns
        --------
        Vcf
            SV records with svtype being BND were integrated into breakpoints, and svtype INFO will be overwritten.
        """
        out = self.copy()
        # Lumpy gives "secondary" flag to breakends instead of "mateid".
        if out._metadata['variantcaller'] == 'lumpy':
            ls_secondary = self.get_table('secondary')['id'].tolist()
            out = out.drop_by_id(ls_secondary)
            out.remove_info_table('secondary')
        
        df_svpos = out.get_table('positions')
        df_svtype = out.get_table('svtype')

        # Breakends of Delly don't have mates.
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
                svlen = 0
            elif row['chrom1'] != row['chrom2']:
                svtype = 'TRA'
                svlen = 0
            elif row['strand1'] == row['strand2']:
                svtype = 'INV'
                svlen = abs(row['pos2'] - row['pos1'])
            elif get_inslen_and_insseq_from_alt(row['alt'])[0] > abs(row['pos1'] - row['pos2']) * 0.5:
                svtype = 'INS'
                svlen = get_inslen_and_insseq_from_alt(row['alt'])[0]
            elif (row['pos1'] < row['pos2']) & (row['strand1'] == '-') & (row['strand2'] == '+'):
                svtype = 'DUP'
                svlen = abs(row['pos2'] - row['pos1'])
            elif (row['pos1'] > row['pos2']) & (row['strand1'] == '+') & (row['strand2'] == '-'):
                svtype = 'DUP'
                svlen = abs(row['pos2'] - row['pos1'])
            else:
                svtype = 'DEL'
                svlen = -abs(row['pos2'] - row['pos1'])
            svlen = int(svlen)
            out.set_value_for_info_by_id('svlen', svid, 0, value=svlen)
            
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

        if out._metadata['variantcaller'] != 'lumpy': # It is enough to exclude only Lumpy's data because Delly's output have been already returned.
            # remove the SV records of mateid.
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


    def classify_manual_svtype(self, definitions=None, ls_conditions=None, ls_names=None, ls_order=None, return_series=True):
        """
        classify_manual_svtype(definitions, ls_conditions, ls_names, ls_order=None)
        Classify SV records by user-defined criteria. A new INFO table named
        'manual_sv_type' will be created.

        Parameters
        ------------
        definitions: path_or_buf, default None
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
        self.add_info_table('manual_sv_type', df_result, number=1, type_='String', description='Custom SV class defined by user')
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

    def _generate_distance_matrix_by_confidence_intervals(self, multiobject, penalty_length=3e9, str_missing=True):
        positions_table = multiobject.get_table("positions")
        N = len(positions_table)
        distance_matrix = np.full((N,N), penalty_length)
        df_cipos = multiobject.get_table('cipos').pivot(values='cipos', index='id', columns='value_idx').astype(int)
        df_ciend = multiobject.get_table('ciend').pivot(values='ciend', index='id', columns='value_idx').astype(int)
        positions_table.set_index('id', inplace=True)
        ind = positions_table.index.to_list()
        ser_pos1_chrom = positions_table['chrom1']
        ser_pos1_left = positions_table['pos1'].add(df_cipos[0])
        ser_pos1_right = positions_table['pos1'].add(df_cipos[1])
        ser_pos1_strand = positions_table['strand1']
        ser_pos2_chrom = positions_table['chrom2']
        ser_pos2_left = positions_table['pos2'].add(df_ciend[0])
        ser_pos2_right = positions_table['pos2'].add(df_ciend[1])
        ser_pos2_strand = positions_table['strand2']
        df_pos1_ci = pd.concat([ser_pos1_chrom, ser_pos1_left, ser_pos1_right, ser_pos1_strand], axis=1).loc[ind]
        df_pos2_ci = pd.concat([ser_pos2_chrom, ser_pos2_left, ser_pos2_right, ser_pos2_strand], axis=1).loc[ind]
        df_pos1_ci.columns = ['chrom', 'chromStart', 'chromEnd', 'strand']
        df_pos2_ci.columns = ['chrom', 'chromStart', 'chromEnd', 'strand']
        df_pos1_ci.reset_index(inplace=True)
        df_pos2_ci.reset_index(inplace=True)
        bed_pos1 = viola.IntervalTreeForMerge(df_pos1_ci, 0)
        bed_pos2 = viola.IntervalTreeForMerge(df_pos2_ci, 0)
        
        id_int = 0
        for idx in range(N):
            chrom1 = df_pos1_ci.loc[idx, 'chrom']
            start1 = df_pos1_ci.loc[idx, 'chromStart']
            end1 = df_pos1_ci.loc[idx, 'chromEnd'] + 1
            strand1 = df_pos1_ci.loc[idx, 'strand']
            chrom2 = df_pos2_ci.loc[idx, 'chrom']
            start2 = df_pos2_ci.loc[idx, 'chromStart']
            end2 = df_pos2_ci.loc[idx, 'chromEnd'] + 1
            strand2 = df_pos2_ci.loc[idx, 'strand']
            hit1_f = set(bed_pos1.query(chrom1, start1, end1).query('strand == @strand1').index)
            hit2_f = set(bed_pos2.query(chrom2, start2, end2).query('strand == @strand2').index)
            hit1_r = set(bed_pos1.query(chrom2, start2, end2).query('strand == @strand2').index)
            hit2_r = set(bed_pos2.query(chrom1, start1, end1).query('strand == @strand1').index)
            hit_result_f = hit1_f & hit2_f
            hit_result_r = hit1_r & hit2_r
            hit_result = hit_result_f | hit_result_r
            for j in hit_result:
                distance_matrix[id_int, j] = 0
                distance_matrix[j, id_int] = 0
            print(idx)
            print(id_int)
            id_int += 1
        return distance_matrix

    
    def merge(self, ls_vcf = [], ls_caller_names = None, mode = 'distance', threshold = 100, linkage = "complete", str_missing = True, integration = True):
        """
        merge(ls_vcf:list, ls_caller_names:list, threshold:float, linkage = "complete", str_missing=True, integration=True)
        Return a merged or integrated vcf object from mulitple  caller's bedpe objects in ls_bedpe

        Parameters
        ----------
        ls_vcf:list
            A list of vcf objects to be merged, which are the same order with ls_caller_names
        ls_caller_names:list
            A list of names of bedpe objects to be merged, which should have self's name as the first element
        mode: {'distance', 'confidence_intervals'}, default 'distance'
            The mode of the merging strategy. 

            * ``'distance'``: Merge SV records by representative SV positions, that is, coordinates of POS field or that of END in the INFO field. If multiple SV positions are within the distance specified in ``threshold`` each other, they will be merged.
            * ``'confidence_intervals'``: Merge SV records according to the confidence intervals reported by SV callers. If confidence intervals of multiple SV records share their genomic coordinates at least 1bp, the will be merged.

        threshold:float, default 100
            Two SVs whose diference of positions is under this threshold are cosidered to be identical.
            This argument is enabled only when ``mode='distance'``.
        linkage:{‘complete’, ‘average’, ‘single’}, default ’complete’
            The linkage of hierarchical clustering.
            To keep the mutual distance of all SVs in each cluster below the threshold, 
            "complete" is recommended.
            This argument is enabled only when ``mode='distance'``.
        str_missing:bool, default True
            If True, all the missing strands are considered to be the same with the others.
        integration:bool, default True 
            If True, vcf objects in ls_vcf will be merged and integrated with priority as ls_caller_names.

        Returns
        ----------
        Vcf
            A merged  vcf object or an integrated vcf object
            
        """
        if self in ls_vcf:
            pass
        else:
            ls_vcf = [self] + ls_vcf
        
        if ls_caller_names is None:
            ls_caller_names = [vcf._metadata["variantcaller"] for vcf in ls_vcf] 

        multivcf = viola.TmpVcfForMerge(ls_vcf, ls_caller_names)
        if mode == 'distance':
            distance_matrix = self._generate_distance_matrix_by_distance(multivcf, penalty_length=3e9, str_missing=str_missing)
            hcl_clustering_model = AgglomerativeClustering(n_clusters=None, affinity="precomputed", linkage=linkage, distance_threshold=threshold)
            labels = hcl_clustering_model.fit_predict(X = distance_matrix)
        elif mode == 'confidence_intervals':
            distance_matrix = self._generate_distance_matrix_by_confidence_intervals(multivcf)
            hcl_clustering_model = AgglomerativeClustering(n_clusters=None, affinity="precomputed", linkage='complete', distance_threshold=100)
            labels = hcl_clustering_model.fit_predict(X = distance_matrix)
        else:
            raise IllegalArgumentError(mode)
        
        
        positions_table = multivcf.get_table("positions")
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
        
        originalid = multivcf.get_table("global_id")["id"]
        df_originalid = pd.DataFrame({"id":positions_table["id"], "value_idx":value_idx, "originalid":originalid})
        
        ############## Edited by Sugita ##################
        df_id = multivcf.get_table("global_id")
        df_patients = multivcf.get_table("patients")
        df_id_patients = df_id.merge(df_patients, left_on="patient_id", right_on="id")
        caller = df_id_patients["patients"]
        df_caller = pd.DataFrame({"id":positions_table["id"], "value_idx":value_idx, "caller":caller})
        ############## /Edited by Sugita #################

        df_pos = multivcf._df_svpos
        df_filters = multivcf._df_filters
        odict_df_info = multivcf._odict_df_info
        df_formats = multivcf._df_formats
        odict_df_headers = multivcf._odict_df_headers
        metadata = OrderedDict({'fileformat': 'VCFv4.2'}) # add by sugita
        args = [df_pos, df_filters, odict_df_info, df_formats, odict_df_headers, metadata] # edited by sugita
        merged_vcf = viola.Vcf(*args)
    
        merged_vcf.add_info_table(table_name="mergedid", table=df_mergedid, number=1, type_="String", description="ID of breakpoints.")
        merged_vcf.add_info_table(table_name="originalid", table=df_originalid, number=1, type_="String", description="The SV-caller-derived ID before merging.")
        merged_vcf.add_info_table(table_name="caller", table=df_caller, number=1, type_="String", description="The name of SV caller which identified the SV record.")

        if integration is False:
            return merged_vcf
        
        if integration is True:
            intergrated_vcf = merged_vcf.integrate(merged_vcf = merged_vcf, priority=ls_caller_names)    
            return intergrated_vcf

    def integrate(self, merged_vcf, priority):
        """
        integrate(merged_vcf:Vcf object, priority:List)
        Return an integrated Vcf object 

        Parameters
        ----------
        merged_vcf:Vcf object
            A Vcf object to be integrated
        priority:list
            Have caller names in list, 
            and shows the priority at which id is adopted in a mergedid cluster.

        Returns
        ----------
        An integrated vcf object
            
        """
        vcf = merged_vcf.copy()
        prior_dict = {}
        for i, c in enumerate(priority):
            prior_dict[c] = i
        mergedid = vcf.get_table("mergedid")["mergedid"]
        globalid = vcf.get_table("positions")["id"]
        caller = vcf.get_table("caller")["caller"]
        priority_num = [prior_dict[c] for c in caller]
        intg_df = pd.DataFrame({"mergedid":mergedid, "supportingid":globalid, "supportingcaller":caller, "priority":priority_num}, 
                                columns=['mergedid', 'supportingid', 'supportingcaller', "priority"])
        
        id_set = set()
        array_dict = {}
        info_str_lst = ["supportingid", "supportingcaller"]
        info_int_lst = ["supportingidcount", "supportingcallercount"]
        info_list = info_str_lst + info_int_lst
        for info in info_list:
            array_dict[info] = np.empty((0, 3))
        
        for bp in intg_df.groupby("mergedid"):
            priority_block = next(iter(bp[1].groupby("priority")))[1]
            id = priority_block.at[priority_block.index[0], "supportingid"]
            id_set.add(id)
            info_dict = {}
            info_dict["supportingid"] = bp[1]["supportingid"].values.reshape(-1,1)
            info_dict["supportingcaller"] = np.array(list(set(bp[1]["supportingcaller"].values))).reshape(-1,1)
            info_dict["supportingidcount"] = np.array(len(bp[1]["supportingid"].values)).reshape(-1,1)
            info_dict["supportingcallercount"] = np.array(len(set(bp[1]["supportingcaller"].values))).reshape(-1,1)
            for info in info_list:
                N = info_dict[info].shape[0]
                block = np.concatenate([np.array([id]*N).reshape(-1,1), np.arange(0, N, dtype=int).reshape(-1,1), info_dict[info]], axis=1)
                array_dict[info] = np.concatenate([array_dict[info], block])
        
        df = {}
        for info in info_list:
            df[info] = pd.DataFrame(data=array_dict[info], columns=["id", "value_idx", info]).astype({'value_idx': int})
        for info in info_int_lst:
            df[info] = df[info].astype({info:int})

        integrated_vcf = vcf.drop_by_id(list(set(globalid) - id_set))
        integrated_vcf.add_info_table(table_name="supportingid", table=df["supportingid"], number=None, 
                                    type_="String", description="IDs of original SV records supporting the merged SV record.")
        integrated_vcf.add_info_table(table_name="supportingcaller", table=df["supportingcaller"], number=None, 
                                    type_="String", description="SV callers supporting the variant.")
        integrated_vcf.add_info_table(table_name="supportingidcount", table=df["supportingidcount"], number=1, 
                                    type_="Integer", description="Number of original SV records supporting the merged SV record.")
        integrated_vcf.add_info_table(table_name="supportingcallercount", table=df["supportingcallercount"], number=1, 
                                    type_="Integer", description="Count of SV callers supporting the variant.")
        return integrated_vcf
        
        




            