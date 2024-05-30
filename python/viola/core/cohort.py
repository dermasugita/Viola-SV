import pandas as pd
import sys, os
from collections import OrderedDict
from viola.core.bedpe import Bedpe
from viola.core.vcf import Vcf
from typing import (
    List,
    Optional,
)

class MultiBedpe(Bedpe):
    """
    A database-like object that contains information of multiple BEDPE files.
    In this class, main keys in most tables are "global id" instead of using
    "SV id" from SV callers. "global id" is unique ID of all the SV record
    across all the samples.
    """
    _internal_attrs = [
        "_df_id",
        "_df_patients",
        "_df_svpos",
        "_odict_df_info",
        "_ls_patients",
        "_ls_infokeys",
        "_odict_alltables",
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
    def __init__(
        self,
        ls_bedpe: List[Bedpe] = None, 
        ls_patient_names: List[str] = None, 
        direct_tables: Optional[List[pd.DataFrame]] = None
        ):
        if direct_tables is None:
            df_id, df_patients, df_svpos, odict_df_info = self.__init__from_ls_bedpe(ls_bedpe, ls_patient_names)
            self.__init__common(df_id, df_patients, df_svpos, odict_df_info)
        else:
            self.__init__common(*direct_tables)
    
    def __init__from_ls_bedpe(self, ls_bedpe, ls_patient_names):
        ls_df_id = []
        ls_df_svpos = []
        dict_ls_df_info = dict() 
        ls_patient_id = [i for i in range(len(ls_patient_names))]
        df_patients = pd.DataFrame({'id': ls_patient_id, 'patients': ls_patient_names})
        for bedpe, patient_id, patient_name in zip(ls_bedpe, ls_patient_id, ls_patient_names):
            df_svpos = bedpe.get_table('positions')
            df_id = df_svpos[['id']].copy()
            df_id['patient_id'] = patient_id
            df_id['global_id'] = str(patient_name) + '_' + df_id['id'].astype(str)
            df_id = df_id[['global_id', 'patient_id', 'id']]
            ls_df_id.append(df_id)

            df_svpos['id'] = str(patient_name) + '_' + df_svpos['id'].astype(str)
            ls_df_svpos.append(df_svpos)

            for key, value in bedpe._odict_df_info.items():
                value = value.copy()
                value['id'] = str(patient_name) + '_' + value['id'].astype(str)
                if dict_ls_df_info.get(key) is None:
                    dict_ls_df_info[key] = [value]
                else:
                    dict_ls_df_info[key].append(value)
        df_concat_id = pd.concat(ls_df_id, ignore_index=True)
        df_concat_svpos = pd.concat(ls_df_svpos, ignore_index=True)
        odict_df_info = OrderedDict()
        for key, value in dict_ls_df_info.items():
            odict_df_info[key] = pd.concat(value)
        
        return (df_concat_id, df_patients, df_concat_svpos, odict_df_info)
    
    def __init__common(self, df_id, df_patients, df_svpos, odict_df_info):
        self._df_id = df_id
        self._df_patients = df_patients
        self._ls_patients = df_patients['patients'].to_list()
        self._df_svpos = df_svpos 
        self._odict_df_info = odict_df_info
        self._ls_infokeys = [x.lower() for x in odict_df_info.keys()]
        ls_keys = ['global_id', 'patients', 'positions'] + self._ls_infokeys
        ls_values = [df_id, df_patients, df_svpos] + list(odict_df_info.values())
        self._odict_alltables = OrderedDict([(k, v) for k, v in zip(ls_keys, ls_values)])
        self._repr_config = {
            'info': None,
        }
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
        desc_doc = 'Documentation of MultiBedpe object ==> '
        doc_link = 'https://dermasugita.github.io/ViolaDocs/docs/html/reference/multi_bedpe.html'
        out = desc_info + str_infokeys + '\n' + desc_doc + doc_link + '\n' + str_df_out
        if return_as_dataframe:
            return df_out
        return str(out)
        

    def filter_by_id(self, arrlike_id):
        """
        filter_by_id(arrlike_id)
        Filter MultiBedpe object according to the list of SV ids.
        Return object is also an instance of the MultiBedpe object

        Parameters
        ---------------
        arrlike_id: list-like
            Global ids which you would like to keep.
        
        Returns
        ---------------
        MultiBedpe
            A MultiBedpe object with the SV id specified in the arrlike_id argument.
            All records associated with SV ids that are not in the arrlike_id will be discarded.
        
        """
        df_global_id = self.get_table('global_id')
        out_global_id = df_global_id.loc[df_global_id['global_id'].isin(arrlike_id)].reset_index(drop=True)
        out_patients = self.get_table('patients')
        out_svpos = self._filter_by_id('positions', arrlike_id)
        out_odict_df_info = OrderedDict([(k, self._filter_by_id(k, arrlike_id)) for k in self._ls_infokeys])
        return MultiBedpe(direct_tables=[out_global_id, out_patients, out_svpos, out_odict_df_info])
    

    def classify_manual_svtype(self, definitions=None, ls_conditions=None, ls_names=None, ls_order=None, return_data_frame=True, exclude_empty_cases=False):
        """
        classify_manual_svtype(definitions, ls_conditions, ls_names, ls_order=None, exclude_empty_cases=False)
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
        exclude_empty_cases: bool, default False
            If True, samples which have no SV record will be excluded.
        
        Returns
        ---------
        pd.DataFrame or None
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
        if return_data_frame:
            if ls_order is None:
                pd_ind_reindex = pd.Index(ls_names + ['others'])
            else:
                pd_ind_reindex = pd.Index(ls_order)
            df_feature_counts = self.get_feature_count_as_data_frame(ls_order=pd_ind_reindex, exclude_empty_cases=exclude_empty_cases)
            return df_feature_counts

    def get_feature_count_as_data_frame(self, feature='manual_sv_type', ls_order=None, exclude_empty_cases=False):
        df_feature = self.get_table(feature)
        df_id = self.get_table('global_id')
        df_patients = self.get_table('patients')
        df_merged = pd.merge(df_feature, df_id, left_on='id', right_on='global_id')
        df_merged = df_merged.merge(df_patients, left_on='patient_id', right_on='id')
        df_feature_counts = df_merged.pivot_table('global_id', index='patients', columns=feature, aggfunc='count', fill_value=0)
        if not exclude_empty_cases:
            df_feature_counts = df_feature_counts.reindex(self._ls_patients, fill_value=0)
        if ls_order is not None:
            pd_ind_reindex = pd.Index(ls_order, name=feature)
            df_feature_counts = df_feature_counts.reindex(columns=pd_ind_reindex, fill_value=0)
        return df_feature_counts

class MultiVcf(Vcf):
    """
    A database-like object that contains information of multiple Vcf files.
    In this class, main keys in most tables are "global id" instead of using
    "SV id" from SV callers. "global id" is unique ID of all the SV record
    across all the samples.
    """
    _internal_attrs = [
        "_df_id",
        "_df_patients",
        "_df_svpos",
        "_odict_df_info",
        "_ls_patients",
        "_ls_infokeys",
        "_odict_alltables",
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
    def __init__(
        self,
        ls_vcf: List[Vcf] = None, 
        ls_patient_names: List[str] = None, 
        direct_tables: Optional[List[pd.DataFrame]] = None
        ):
        if direct_tables is None:
            df_id, df_patients, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers  = self.__init__from_ls_vcf(ls_vcf, ls_patient_names)
            self.__init__common(df_id, df_patients, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers)
        else:
            self.__init__common(*direct_tables)
    
    def __init__from_ls_vcf(self, ls_vcf, ls_patient_names):
        ls_df_id = []
        ls_df_svpos = []
        ls_df_filters = []
        odict_ls_df_info = OrderedDict() 
        ls_df_formats = []
        odict_ls_df_headers = OrderedDict()

        # Header Integration
        for vcf, patient_name in zip(ls_vcf, ls_patient_names):
            for key, value in vcf._odict_df_headers.items():
                value = value.copy()
                if odict_ls_df_headers.get(key) is None:
                    odict_ls_df_headers[key] = [value]
                else:
                    odict_ls_df_headers[key].append(value)

        odict_df_headers = OrderedDict()
        for key, value in odict_ls_df_headers.items():
            for idx, df in enumerate(value):
                if idx == 0:
                    df_merged = df
                    continue
                on = list(df_merged.columns)
                df_merged = df_merged.merge(df, how='outer', on=on)
            
            odict_df_headers[key] = df_merged
        # /Header Integration

        ls_patient_id = [i for i in range(len(ls_patient_names))]
        df_patients = pd.DataFrame({'id': ls_patient_id, 'patients': ls_patient_names})
        for vcf, patient_id, patient_name in zip(ls_vcf, ls_patient_id, ls_patient_names):
            df_svpos = vcf.get_table('positions')
            df_filters = vcf.get_table('filters')
            df_formats = vcf.get_table('formats')
            df_id = df_svpos[['id']].copy()
            df_id['patient_id'] = patient_id
            df_id['global_id'] = str(patient_name) + '_' + df_id['id'].astype(str)
            df_id = df_id[['global_id', 'patient_id', 'id']]
            ls_df_id.append(df_id)

            df_svpos['id'] = str(patient_name) + '_' + df_svpos['id'].astype(str)
            ls_df_svpos.append(df_svpos)

            df_filters['id'] = str(patient_name) + '_' + df_filters['id'].astype(str)
            ls_df_filters.append(df_filters)

            df_formats['id'] = str(patient_name) + '_' + df_formats['id'].astype(str)
            ls_df_formats.append(df_formats)

            for info in odict_df_headers['infos_meta'].id:
                df_info_ = vcf._odict_df_info.get(info, None)
                if df_info_ is None:
                    df_info = pd.DataFrame(columns=('id', 'value_idx', info.lower()))
                else:
                    df_info = df_info_.copy()
                    df_info['id'] = str(patient_name) + '_' + df_info['id'].astype(str)
                if odict_ls_df_info.get(info) is None:
                    odict_ls_df_info[info] = [df_info]
                else:
                    odict_ls_df_info[info].append(df_info)

            
        df_concat_id = pd.concat(ls_df_id, ignore_index=True)
        df_concat_svpos = pd.concat(ls_df_svpos, ignore_index=True)
        df_concat_filters = pd.concat(ls_df_filters, ignore_index=True)
        df_concat_formats = pd.concat(ls_df_formats, ignore_index=True)
        odict_df_info = OrderedDict()

        for key, value in odict_ls_df_info.items():
            odict_df_info[key] = pd.concat(value)
        
        return (df_concat_id, df_patients, df_concat_svpos, df_concat_filters, odict_df_info, df_concat_formats, odict_df_headers)

    
    def __init__common(self, df_id, df_patients, df_svpos, df_filters, odict_df_info, df_formats, odict_df_headers = {}):
        self._df_id = df_id
        self._df_patients = df_patients
        self._df_svpos = df_svpos
        self._df_filters = df_filters
        self._odict_df_info = odict_df_info
        self._df_formats = df_formats
        self._odict_df_headers = odict_df_headers
        self._ls_patients = df_patients['patients'].to_list()
        self._ls_infokeys = [ x.lower() for x in odict_df_headers['infos_meta']['id'].tolist()]
        ls_keys = ['global_id', 'patients', 'positions', 'filters'] + self._ls_infokeys + ['formats'] + \
        list(odict_df_headers.keys())
        ls_values = [df_id, df_patients, df_svpos, df_filters] + list(odict_df_info.values()) + [df_formats] + list(odict_df_headers.values())
        self._odict_alltables = OrderedDict([(k, v) for k, v in zip(ls_keys, ls_values)])
        self._repr_config = {
            'info': None,
        }

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
        desc_doc = 'Documentation of MultiBedpe object ==> '
        doc_link = 'https://dermasugita.github.io/ViolaDocs/docs/html/reference/multi_vcf.html'
        out = desc_info + str_infokeys + '\n' + desc_doc + doc_link + '\n' + str_df_out
        if return_as_dataframe:
            return df_out
        return str(out)
        

    def filter_by_id(self, arrlike_id):
        """
        filter_by_id(arrlike_id)
        Filter MultiBedpe object according to the list of SV ids.
        Return object is also an instance of the MultiBedpe object

        Parameters
        ---------------
        arrlike_id: list-like
            Global ids which you would like to keep.
        
        Returns
        ---------------
        MultiBedpe
            A MultiBedpe object with the SV id specified in the arrlike_id argument.
            All records associated with SV ids that are not in the arrlike_id will be discarded.
        
        """
        df_global_id = self.get_table('global_id')
        out_global_id = df_global_id.loc[df_global_id['global_id'].isin(arrlike_id)].reset_index(drop=True)
        out_patients = self.get_table('patients')
        out_svpos = self._filter_by_id('positions', arrlike_id)
        out_filters = self._filter_by_id('filters', arrlike_id)
        out_odict_df_info = OrderedDict([(k, self._filter_by_id(k, arrlike_id)) for k in self._ls_infokeys])
        out_formats = self._filter_by_id('formats', arrlike_id)
        out_odict_df_headers = self._odict_df_headers.copy()
        return MultiVcf(direct_tables=[out_global_id, out_patients, out_svpos, out_filters, out_odict_df_info, out_formats, out_odict_df_headers])
    

    def classify_manual_svtype(self, definitions=None, ls_conditions=None, ls_names=None, ls_order=None, return_data_frame=True, exclude_empty_cases=False):
        """
        classify_manual_svtype(definitions, ls_conditions, ls_names, ls_order=None, exclude_empty_cases=False)
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
        exclude_empty_cases: bool, default False
            If True, samples which have no SV record will be excluded.
        
        Returns
        ---------
        pd.DataFrame or None
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
        if return_data_frame:
            if ls_order is None:
                pd_ind_reindex = pd.Index(ls_names + ['others'])
            else:
                pd_ind_reindex = pd.Index(ls_order)
            df_feature_counts = self.get_feature_count_as_data_frame(ls_order=pd_ind_reindex, exclude_empty_cases=exclude_empty_cases)
            return df_feature_counts
    
    def as_bedpe_multi(self):
        """
        as_bedpe_multi()
        Convert MultiVcf object into MultiBedpe object.
        
        Returns
        --------
        MultiBedpe

        Notes
        ------
        This process is lossy because only the "global_id", "patients", "positions" and INFO tables are inherited by MultiBedpe class.

        Examples
        ---------
        The file tree of this script example. You can't run this code itself, sorry!

        .. code-block::

            ├── this_script.py
            └── vcf
                ├── patient1.vcf
                └── patient2.vcf
        
        Here is the example code.

        >>> import viola
        >>> multi_vcf = viola.read_vcf_multi('./vcf', variant_caller='manta')
        >>> multi_bedpe = multi_vcf.as_bedpe()
        >>> multi_bedpe
        INFO=imprecise,svtype,svlen,end,cipos,ciend,cigar,mateid,event,homlen,homseq,svinslen,svinsseq,left_svinsseq,right_svinsseq,contig,bnd_depth,mate_bnd_depth,somatic,somaticscore,junction_somaticscore,inv3,inv5
        Documentation of MultiBedpe object ==> https://dermasugita.github.io/ViolaDocs/docs/html/reference/multi_bedpe.html
                id              be1              be2 strand  qual svtype
        0    patient1_test1    chr1:82550461    chr1:82554226     +-  None    DEL
        1    patient1_test2    chr1:22814217    chr1:92581132     --  None    INV
        2    patient2_test1    chr1:60567906    chr1:60675941     +-  None    DEL
        3    patient2_test2    chr1:69583190    chr1:69590948     +-  None    DEL
        5  patient2_test3_1  chr11:111134697   chr17:26470495     +-  None    BND
        6  patient2_test3_2   chr17:26470495  chr11:111134697     -+  None    BND
        """
        df_id = self.get_table('global_id')
        df_patients = self.get_table('patients')
        df_svpos = self.get_table('positions')
        odict_df_info_view = self._odict_df_info
        odict_df_info = OrderedDict((k, v.copy()) for k, v in odict_df_info_view.items())
        multi_bedpe = MultiBedpe(direct_tables=[df_id, df_patients, df_svpos, odict_df_info])
        return multi_bedpe

    def as_bedpe(self):
        """
        as_bedpe()
        The same as as_bedpe_multi()
        """
        return self.as_bedpe_multi()

    
    def get_feature_count_as_data_frame(self, feature='manual_sv_type', ls_order=None, exclude_empty_cases=False):
        df_feature = self.get_table(feature)
        df_id = self.get_table('global_id')
        df_patients = self.get_table('patients')
        df_merged = pd.merge(df_feature, df_id, left_on='id', right_on='global_id')
        df_merged = df_merged.merge(df_patients, left_on='patient_id', right_on='id')
        df_feature_counts = df_merged.pivot_table('global_id', index='patients', columns=feature, aggfunc='count', fill_value=0)
        if not exclude_empty_cases:
            df_feature_counts = df_feature_counts.reindex(self._ls_patients, fill_value=0)
        if ls_order is not None:
            pd_ind_reindex = pd.Index(ls_order, name=feature)
            df_feature_counts = df_feature_counts.reindex(columns=pd_ind_reindex, fill_value=0)
        return df_feature_counts