import vcf
import numpy as np
import pandas as pd
import re
import itertools
import urllib.request
from collections import OrderedDict
import warnings
from typing import (
    Union,
    Optional,
)
from io import StringIO
from viola.core.bedpe import Bedpe
from viola.core.vcf import Vcf
from viola.core.bed import Bed
from viola.utils.utils import is_url
from viola.io._vcf_parser import (
    read_vcf_manta,
    read_vcf_delly,
    read_vcf_lumpy,
    read_vcf_gridss,
)
pd.set_option('display.max_columns', 10)
pd.set_option('display.max_colwidth', 30)
pd.set_option('display.width', 1000)


class _VcfReader():
    _singular_metadata = ['fileformat', 'fileDate', 'reference', 'variantcaller']
    def __init__(self, variant_caller, patient_name):
        self.variant_caller = variant_caller
        self.patient_name = patient_name
        self.metadata = OrderedDict()
        self._svtype = ''
        self.odict_odict_headers = OrderedDict()
        self.odict_filters = OrderedDict({'id': [], 'filter': []})
        self.odict_odict_infos = OrderedDict()
        self.odict_formats = OrderedDict({'id': [], 'sample': [], \
            'format': [], 'value_idx': [], 'value': []})
        self.odict_svpos = OrderedDict({
            'id': [],
            'chrom1': [],
            'pos1': [],
            'chrom2': [],
            'pos2': [],
            'strand1': [],
            'strand2': [],
            'qual': [],
            'ref': [],
            'alt': [],
            'svtype': [],
        })
        self.re_header_split = re.compile(r"^(?P<key>[^=]*)=(?P<items>.*)")
        self.re_header_items = re.compile(r'([^<=,\s]*)=([^,\s]*|".*")[,>]')

    def _refine_items(self, key, dtype, value):
        value = value.replace('"', '')
        if dtype is None:
            if value.isdigit():
                value = int(value)
            elif key == 'GT' and value == '.':
                value = '.'
            elif value == ".":
                value = None
            elif value.replace('.', '', 1).isdigit():
                value = float(value)
            elif key == 'Number' and value == 'A':
                value = -1
            elif key == 'Number' and value == 'G':
                value = -2
        else:
            if dtype == 'Integer':
                value = int(value)
            elif dtype == 'Float':
                value = float(value)
            elif dtype == 'String':
                value = str(value)
            elif dtype == 'Flag':
                value = bool(value)
        return value

    def _vcf_header_parser(self, line):
        header_content = line[2:]
        header_split = self.re_header_split.match(header_content)
        header_key = header_split.group('key')
        header_items = header_split.group('items')
        if not header_items.startswith('<'):
            if header_key in self._singular_metadata:
                self.metadata[header_key] = header_items
            else:
                if header_key not in self.metadata:
                    self.metadata[header_key] = []
                self.metadata[header_key].append(header_items)
            return
        ls_tuple_header_items = self.re_header_items.findall(header_items)

        if header_key == 'contig':
            header_key = 'contigs_meta'
        elif header_key == 'INFO':
            header_key = 'infos_meta'
        elif header_key == 'FILTER':
            header_key = 'filters_meta'
        elif header_key == 'FORMAT':
            header_key = 'formats_meta'
        elif header_key == 'ALT':
            header_key = 'alts_meta'


        odict_header_items = \
            OrderedDict([(k.lower(), [self._refine_items(k, None, v)]) for k, v in ls_tuple_header_items])
        
        if self.variant_caller == 'delly':
            if header_key == 'infos_meta' and odict_header_items['id'][0] == 'SVLEN':
                odict_header_items.update({'id': ['SVLENORG', 'SVLEN']})
                odict_header_items['number'].append(None)
                odict_header_items['type'].append('Integer')
                odict_header_items['description'].append('Length of SV')
                if 'source' in odict_header_items:
                    odict_header_items['source'].append(None)
                if 'version' in odict_header_items:
                    odict_header_items['version'].append(None)
        elif self.variant_caller == 'lumpy':
            if header_key == 'infos_meta' and odict_header_items['id'][0] == 'SU':
                odict_header_items.update({'id': ['SU', 'SUORG']})
                odict_header_items['number'].append(None)
                odict_header_items['type'].append('Integer')
                odict_header_items['description'].append('Original value of Lumpy SU')
                if 'source' in odict_header_items:
                    odict_header_items['source'].append(None)
                if 'version' in odict_header_items:
                    odict_header_items['version'].append(None)

        
        if self.odict_odict_headers.get(header_key) is None:
            self.odict_odict_headers[header_key] = \
                odict_header_items
        else:
            for k, v in odict_header_items.items():
                if self.odict_odict_headers[header_key].get(k) is None:
                    self.odict_odict_headers[header_key][k] = v
                else:
                    self.odict_odict_headers[header_key][k].extend(v)
        
    def _header_df_constructor(self):
        odict_df_headers = OrderedDict()
        for k, v in self.odict_odict_headers.items():
            if k == 'infos_meta':
                df = pd.DataFrame(v, columns=('id', 'number', 'type', 'description', 'source', 'version'))
            else:
                df = pd.DataFrame(v)
            
            if 'number' in df.columns:
                df['number'] = df['number'].astype(object) # backward compatibility
            odict_df_headers[k] = df
        if self.variant_caller == 'lumpy':
            odict_df_headers['contigs_meta'] = pd.DataFrame(columns=('id', 'length'))
            odict_df_headers['filters_meta'] = pd.DataFrame(columns=('id', 'description'))
        self.odict_df_headers = odict_df_headers
    
    def _sample_extractor(self, header_line):
        header_line = header_line.replace('\n', '')
        ls_header_line = header_line.split('\t')
        ls_samples = ls_header_line[9:]
        self.ls_samples = ls_samples
        self.odict_odict_headers['samples_meta'] = OrderedDict({'id': ls_samples})
    
    def _filter_parser(self, ls_line):
        svid = ls_line[2]
        ls_filter_ = ls_line[6].split(';')
        ls_filter_ = ['PASS' if filter_ == '.' else filter_ for filter_ in ls_filter_]
        len_filter_ = len(ls_filter_)
        if self.variant_caller == 'lumpy' and self._svtype == 'INV':
            svid1 = str(svid) + '_1'
            svid2 = str(svid) + '_2'
            ls_svid = [svid1, svid2] * len_filter_
            ls_filter_ = ls_filter_ * 2
        else:
            ls_svid = [svid] * len_filter_
        self.odict_filters['id'].extend(ls_svid)
        self.odict_filters['filter'].extend(ls_filter_)
        
    def _get_sv_type(self, ls_info):
        for info in ls_info:
            ls_key_item = info.split('=')
            k = ls_key_item[0]
            if k == 'SVTYPE':
                self._svtype = ls_key_item[1]
                return ls_key_item[1]
    
    def _info_parser(self, ls_line):
        odict_infos = OrderedDict()
        odict_infos_out = OrderedDict()
        svid = ls_line[2]
        info_line = ls_line[7]
        ls_info = info_line.split(';')
        if self.variant_caller == 'lumpy':
            is_inv = self._get_sv_type(ls_info) == 'INV'
        
        def _append_info_to_odict(svid, ls_key_item, k, info_dtype):
            if len(ls_key_item) == 1:
                ls_v = [True]
                odict_infos = OrderedDict({'id': [svid], \
                    'value_idx': [0], k: ls_v})
            else:
                v = ls_key_item[1]
                ls_v = v.split(',')
                ls_v = [self._refine_items(None, info_dtype, item) for item in ls_v]
                len_v = len(ls_v)
                ls_svid = [svid] * len_v
                odict_infos = OrderedDict({'id': ls_svid, \
                    'value_idx': list(range(len_v)), k: ls_v})
            odict_infos_out[k] = ls_v
            if self.odict_odict_infos.get(k) is None:
                self.odict_odict_infos[k] = odict_infos
            else:
                for table_key, table_value in odict_infos.items():
                    self.odict_odict_infos[k][table_key].extend(table_value)
        
        for info in ls_info:
            ls_key_item = info.split('=')
            k = ls_key_item[0]
            info_dtype_idx = self.odict_odict_headers['infos_meta']['id'].index(k)
            info_dtype = self.odict_odict_headers['infos_meta']['type'][info_dtype_idx]
            k = k.lower()
            if self.variant_caller == 'delly' and k == 'svlen':
                k = 'svlenorg'
            if self.variant_caller == 'lumpy' and is_inv:
                if k == 'strands':
                    ls_v = ls_key_item[1].split(',')
                    for idx_v, each_v in enumerate(ls_v):
                        svid_ = str(svid) + '_' + str(idx_v + 1)
                        su = each_v.split(':')[1]
                        k_su = 'su'
                        ls_key_item = [k, each_v]
                        ls_key_item_su = ['su', su]
                        info_dtype_su = 'Integer'
                        _append_info_to_odict(svid_, ls_key_item, k, info_dtype)
                        _append_info_to_odict(svid_, ls_key_item_su, k_su, info_dtype_su)
                else:
                    if k == 'su':
                        k = 'suorg'
                        ls_key_item[0] = k
                    svid1 = str(svid) + '_1'
                    svid2 = str(svid) + '_2'
                    _append_info_to_odict(svid1, ls_key_item, k, info_dtype)
                    _append_info_to_odict(svid2, ls_key_item, k, info_dtype)
                continue
                    
            _append_info_to_odict(svid, ls_key_item, k, info_dtype)
        
        if self.variant_caller == 'delly':
            svtype = odict_infos_out['svtype'][0]
            pos1 = int(ls_line[1])
            end = odict_infos_out['end'][0]
            if svtype == 'BND':
                svlen = 0
            elif svtype == 'TRA':
                svlen = 0
            elif svtype == 'DEL':
                svlen = pos1 - end
            elif svtype == 'DUP':
                svlen = end - pos1
            elif svtype == 'INV':
                svlen = end - pos1
            else:
                svlen = 0
            if 'svlen' not in self.odict_odict_infos:
                self.odict_odict_infos['svlen'] = OrderedDict({
                    'id': [],
                    'value_idx': [],
                    'svlen': []
                })
            self.odict_odict_infos['svlen']['id'].extend([svid])
            self.odict_odict_infos['svlen']['value_idx'].extend([0])
            self.odict_odict_infos['svlen']['svlen'].extend([svlen])
            odict_infos_out['svlen'] = [svlen]
        elif self.variant_caller == 'lumpy' and is_inv:
            odict_infos_out['event'] = [svid]
            svid1 = str(svid) + '_1'
            svid2 = str(svid) + '_2'
            if 'event' not in self.odict_odict_infos:
                self.odict_odict_infos['event'] = OrderedDict({
                    'id': [svid1, svid2],
                    'value_idx': [0, 0],
                    'event': [svid, svid]
                })
            else:
                self.odict_odict_infos['event']['id'].extend([svid1, svid2])
                self.odict_odict_infos['event']['value_idx'].extend([0, 0])
                self.odict_odict_infos['event']['event'].extend([svid, svid])
        elif self.variant_caller == 'gridss' and odict_infos.get('cipos') is None:
            odict_infos_out['cipos'] = [0, 0]
            if self.odict_odict_infos.get('cipos') is None:
                self.odict_odict_infos['cipos'] = OrderedDict({
                    'id': [svid, svid],
                    'value_idx': [0, 1],
                    'cipos': [0, 0]
                })
            else:
                self.odict_odict_infos['cipos']['id'].extend([svid, svid])
                self.odict_odict_infos['cipos']['value_idx'].extend([0, 1])
                self.odict_odict_infos['cipos']['cipos'].extend([0, 0])
        
        return odict_infos_out
    
    def _format_parser(self, ls_line):
        svid = ls_line[2]
        format_ = ls_line[8]
        ls_formats_names = format_.split(':')
        idx_sample = 9
        for sample in self.ls_samples:
            each_format = ls_line[idx_sample]
            ls_each_format = each_format.split(':')
            for formats_name, each_format_values in zip(ls_formats_names, ls_each_format):
                formats_dtype_idx = self.odict_odict_headers['formats_meta']['id'].index(formats_name)
                formats_dtype = self.odict_odict_headers['formats_meta']['type'][formats_dtype_idx]
                ls_each_format_values = each_format_values.split(',')
                ls_each_format_values = [self._refine_items(formats_name, formats_dtype, item) for item in ls_each_format_values]
                len_each_format_values = len(ls_each_format_values)
                if self.variant_caller == 'lumpy' and self._svtype == 'INV':
                    svid1 = str(svid) + '_1'
                    svid2 = str(svid) + '_2'
                    ls_svid = [svid1, svid2] * len_each_format_values
                    ls_sample = [sample] * len_each_format_values * 2
                    ls_format_names = [formats_name] * len_each_format_values * 2
                    ls_value_idx = list(range(len_each_format_values)) * 2
                    ls_each_format_values *= 2
                else:
                    ls_svid = [svid] * len_each_format_values
                    ls_sample = [sample] * len_each_format_values
                    ls_format_names = [formats_name] * len_each_format_values
                    ls_value_idx = list(range(len_each_format_values))
                self.odict_formats['id'].extend(ls_svid)
                self.odict_formats['sample'].extend(ls_sample)
                self.odict_formats['format'].extend(ls_format_names)
                self.odict_formats['value_idx'].extend(ls_value_idx)
                self.odict_formats['value'].extend(ls_each_format_values)
            idx_sample += 1
    
    def _alt_parser_for_bnd(self, alt):
        if re.search('[\[\]]', alt) is not None:
            if (alt[0] == '[') or (alt[0] == ']'):
                strand1 = '-'
            else:
                strand1 = '+'
            if re.search('\[', alt) is not None:
                strand2 = '-'
            else:
                strand2 = '+'
            ls_alt = re.split('[\[\]]', alt)
            chrom2, pos2 = ls_alt[1].split(':')
            pos2 = int(pos2)
        elif (alt[0] == '.') or (alt[-1] == '.'):
            chrom2 = None
            pos2 = None
            if alt[0] == '.':
                strand1 = '-'
            else:
                strand1 = '+'
            strand2 = None
        
        return chrom2, pos2, strand1, strand2

    
    def _positions_parser(self, ls_line, odict_infos):
        chrom1 = ls_line[0]
        pos1 = int(ls_line[1])
        svid = ls_line[2]
        ref = ls_line[3]
        alt = ls_line[4]
        qual = self._refine_items(None, None, ls_line[5])

        svtype = odict_infos['svtype'][0]
        if svtype == 'DEL':
            chrom2 = chrom1
            pos2 = int(odict_infos['end'][0]) + 1
            strand1 = '+'
            strand2 = '-'
        elif svtype == 'DUP':
            pos1 += 1
            chrom2 = chrom1
            pos2 = int(odict_infos['end'][0])
            strand1 = '-'
            strand2 = '+'
        elif svtype == 'INV':
            chrom2 = chrom1
            pos2 = int(odict_infos['end'][0])
            if self.variant_caller == 'manta':
                if odict_infos.get('inv3', False):
                    strand1 = '+'
                    strand2 = '+'
                else:
                    pos1 += 1
                    pos2 += 1
                    strand1 = '-'
                    strand2 = '-'
            elif self.variant_caller == 'delly':
                if odict_infos['ct'][0] == '3to3':
                    strand1, strand2 = '+', '+'
                else:
                    strand1, strand2 = '-', '-'

        elif svtype == 'INS':
            chrom2 = chrom1
            pos2 = int(pos1)
            strand1 = '+'
            strand2 = '-'
        elif svtype == 'BND':
            chrom2, pos2, strand1, strand2 = self._alt_parser_for_bnd(alt)

        keys = ['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', \
            'strand2', 'qual', 'ref', 'alt', 'svtype']
        
        if self.variant_caller == 'lumpy' and svtype == 'INV':
            svid1 = str(svid) + '_1'
            svid2 = str(svid) + '_2'
            items1 = [svid1, chrom1, pos1, chrom2, pos2, '+', \
                '+', qual, ref, alt, svtype]
            items2 = [svid2, chrom1, pos1+1, chrom2, pos2+1, '-', \
                '-', qual, ref, alt, svtype]
            for key, item in zip(keys, items1):
                self.odict_svpos[key].append(item)
            for key, item in zip(keys, items2):
                self.odict_svpos[key].append(item)

        else:
            items = [svid, chrom1, pos1, chrom2, pos2, strand1, \
                strand2, qual, ref, alt, svtype]

            for key, item in zip(keys, items):
                self.odict_svpos[key].append(item)
    
    def _main_line_parser(self, line):
        line = line.replace('\n', '')
        ls_line = line.split('\t')
        odict_infos = self._info_parser(ls_line)
        self._filter_parser(ls_line)
        self._format_parser(ls_line)
        self._positions_parser(ls_line, odict_infos)

    
    def _filter_df_constructor(self):
        df_filters = pd.DataFrame(self.odict_filters)
        self.df_filters = df_filters

    def _info_df_constructor(self):
        odict_df_infos = OrderedDict()
        ls_info_names = self.odict_odict_headers['infos_meta']['id']
        ls_info_dtype = self.odict_odict_headers['infos_meta']['type']
        for k, dtype in zip(ls_info_names, ls_info_dtype):
            k = k.lower()
            if self.odict_odict_infos.get(k) is None:
                odict_df_infos[k] = pd.DataFrame(columns=('id', 'value_idx', k))
            else:
                odict_info = self.odict_odict_infos[k]
                df_info = pd.DataFrame(odict_info)
                if dtype == 'Integer':
                    df_info[k] = df_info[k].astype(int)
                odict_df_infos[k] = df_info
        #### Generate CIEND table by merging MATEID and CIPOS
        if self.variant_caller == 'lumpy':
            pass
        else:
            if 'mateid' in odict_df_infos:
                df_mateid = odict_df_infos['mateid']
                df_cipos = odict_df_infos['cipos']
                df_merged = pd.merge(df_cipos, df_mateid, on='id')
                df_merged = df_merged[['mateid', 'value_idx_x', 'cipos']]
                df_merged.columns = ['id', 'value_idx', 'ciend']
                df_ciend = odict_df_infos['ciend']
                if df_ciend.empty:
                    odict_df_infos['ciend'] = df_merged
                else:
                    odict_df_infos['ciend'] = pd.concat([df_ciend, df_merged], ignore_index=True)
        #### /Generate CIEND table by merging MATEID and CIPOS

        self.odict_df_infos = odict_df_infos
    
    def _format_df_constructor(self):
        df_formats = pd.DataFrame(self.odict_formats)
        self.df_formats = df_formats
    
    def _position_df_constructor(self):
        df_svpos = pd.DataFrame(self.odict_svpos)
        self.df_svpos = df_svpos
    

        
def read_vcf2(filepath_or_buffer, variant_caller, patient_name=None):
    if patient_name is None:
        warnings.warn(
            'Passing NoneType to the "patient_name" argument is deprecated.',
            DeprecationWarning
        )
    reader = _VcfReader(variant_caller, patient_name)
    if isinstance(filepath_or_buffer, StringIO):
        f = filepath_or_buffer
    else:
        f = open(filepath_or_buffer, 'r')

    for line in f:
        if line.startswith('##'):
            reader._vcf_header_parser(line)
            continue
        elif line.startswith('#'):
            reader._sample_extractor(line)
            continue
        reader._main_line_parser(line)
        
    reader._header_df_constructor()
    reader._filter_df_constructor()
    reader._info_df_constructor()
    reader._format_df_constructor()
    reader._position_df_constructor()
    reader.metadata['variantcaller'] = variant_caller

    ############## Generate contigs_meta ####################
    if variant_caller == 'lumpy':
        arr_contigs = np.append(reader.df_svpos['chrom1'].values,  reader.df_svpos['chrom2'].values)
        arr_contigs = np.unique(arr_contigs)
        arr_contigs_length = np.repeat('.', len(arr_contigs))
        df_contigs_meta = pd.DataFrame({'id': arr_contigs, 'length': arr_contigs_length})
        reader.odict_df_headers['contigs_meta'] = df_contigs_meta
    
    reader.odict_df_infos = OrderedDict([(k.upper(), v) for k, v in reader.odict_df_infos.items()])
   

    args = [
        reader.df_svpos,
        reader.df_filters,
        reader.odict_df_infos,
        reader.df_formats,
        reader.odict_df_headers,
        reader.metadata,
        patient_name
    ]

    return Vcf(*args)



def read_vcf(filepath_or_buffer: Union[str, StringIO], variant_caller: str = "manta", patient_name = None):
    """
    read_vcf(filepath_or_buffer, variant_callser = "manta", patient_name = None)
    Read vcf file of SV and return Vcf object.

    Parameters
    ---------------
    filepath_or_buffer: str or StringIO
        String path to the vcf file. StringIO is also acceptable.
        (Wether URL is acceptable or not hasn't been tested.)
        (Acceptable types should be extended in the future)
    variant_caller: str
        Let this function know which SV caller was used to create vcf file.
    patient name: str or None, default None
    
    Returns
    ---------------
    A Vcf object
    """
    if patient_name is None:
        warnings.warn(
            'Passing NoneType to the "patient_name" argument is deprecated.',
            DeprecationWarning
        )
    # read vcf files using PyVcf package
    if isinstance(filepath_or_buffer, str) and is_url(filepath_or_buffer):
        b = StringIO(urllib.request.urlopen(filepath_or_buffer).read().decode('utf-8'))
        vcf_reader = vcf.Reader(b)
    elif isinstance(filepath_or_buffer, str):
        vcf_reader = vcf.Reader(open(filepath_or_buffer, 'r'))
    elif isinstance(filepath_or_buffer, StringIO):
        vcf_reader = vcf.Reader(filepath_or_buffer)
    else:
        vcf_reader = vcf.Reader(filepath_or_buffer)
        #raise TypeError("should be file or buffer")

    if variant_caller == 'manta':
        return read_vcf_manta(vcf_reader, patient_name=patient_name)
    elif variant_caller == 'delly':
        return read_vcf_delly(vcf_reader, patient_name)
    elif variant_caller == 'lumpy':
        return read_vcf_lumpy(vcf_reader, patient_name)
    elif variant_caller == 'gridss':
        return read_vcf_gridss(vcf_reader, patient_name)

    # obtain header informations
    odict_contigs = vcf_reader.contigs
    df_contigs_meta = pd.DataFrame(odict_contigs, index=('id', 'length')).T.reset_index(drop=True)

    # obtain alteration class informations
    odict_alts = vcf_reader.alts
    df_alts_meta = pd.DataFrame(odict_alts, index=('id', 'description')).T.reset_index(drop=True)

    # obtain info field information
    odict_infos = vcf_reader.infos
    df_infos_meta = pd.DataFrame(odict_infos, index=('id', 'number', 'type', 'description', 'source', 'version')).T.reset_index(drop=True)

    # obtain FORMAT informations
    odict_formats = vcf_reader.formats
    df_formats_meta = pd.DataFrame(odict_formats, index=('id', 'number', 'type','description')).T.reset_index(drop=True)

    # obtain FILTER informations
    odict_filters = vcf_reader.filters
    df_filters_meta = pd.DataFrame(odict_filters, index=('id', 'description')).T.reset_index(drop=True)

    # obtain SAMPLE informations
    ls_samples = vcf_reader.samples
    df_samples = pd.DataFrame(ls_samples, columns=['id'])

    odict_df_headers = OrderedDict(
        contigs_meta = df_contigs_meta,
        alts_meta = df_alts_meta, 
        infos_meta = df_infos_meta, 
        formats_meta = df_formats_meta, 
        filters_meta = df_filters_meta, 
        samples_meta = df_samples)

    
    ls_pos = []
    ls_filters = []
    dict_infos = {k: [] for k in df_infos_meta.id}
    ls_df_formats = []

    
    for record in vcf_reader:
        dict_row = {}
        # ex: Record(CHROM=chr1, POS=9163435, REF=G, ALT=[<DUP:TANDEM>])
        row_ID = record.ID
        row_CHROM1 = record.CHROM
        row_POS1 = record.POS
        row_REF = record.REF
        row_ALT = record.ALT[0] # This operation is safe when parsing manta vcf files, but potential cause of information loss in another callers.
        row_QUAL = record.QUAL
         
        
        #####FILTER
        ls_filter = record.FILTER
        if ls_filter is None:
            ls_filter = ['PASS']
        elif len(ls_filter) == 0:
            ls_filter = ['PASS'] 
        row_FILTER = ls_filter
        for filter_ in row_FILTER:
            ls_filters.append({'id': row_ID, 'filter': filter_})
        #####/FILTER 

        #####INFO
        row_INFO = record.INFO
        for info in df_infos_meta.id:
            values = row_INFO.get(info, 'none')
            if values == 'none':
                continue
            if not isinstance(values, list):
                values = [values]
            ls_keys = ['id', 'value_idx', info.lower()] 
            for idx, value in enumerate(values):
                ls_value = [row_ID, idx, value]
                dict_a_info = {k: v for k, v in zip(ls_keys, ls_value)}
                dict_infos[info].append(dict_a_info)
        #####/INFO

        ###POS
        if isinstance(row_ALT, vcf.model._Breakend):
            row_CHROM2 = row_ALT.chr
            row_POS2 = row_ALT.pos
            row_STRAND1 = '-' if row_ALT.orientation else '+'
            row_STRAND2 = '-' if row_ALT.remoteOrientation else '+'
        elif row_INFO['SVTYPE'] == 'INV': # manta-specific operation
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
            # manta
            if variant_caller == "manta":
                if row_INFO.get('INV3', False):
                    row_STRANDs = '+'
                else:
                    row_POS1 += 1
                    row_POS2 += 1
                    row_STRANDs = '-'
            # /manta
            # delly
            elif variant_caller == "delly":
                if row_INFO['CT'] == '3to3':
                    row_STRANDs = '+'
                else:
                    row_POS1 += 1
                    row_POS2 += 1
                    row_STRANDs = '-'
            # /delly
            row_STRAND1, row_STRAND2 = row_STRANDs, row_STRANDs
        elif row_INFO['SVTYPE'] == 'DEL':
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END'] + 1
            row_STRAND1 = '+'
            row_STRAND2 = '-'
        elif row_INFO['SVTYPE'] == 'DUP':
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
            row_POS1 += 1
            row_STRAND1 = '-'
            row_STRAND2 = '+'
        else:
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
            row_STRAND1 = '.'
            row_STRAND2 = '.'

        row_SVTYPE = row_INFO['SVTYPE']
        ls_pos.append({
                'id': row_ID, 
                'chrom1': row_CHROM1, 
                'pos1': row_POS1, 
                'chrom2': row_CHROM2, 
                'pos2': row_POS2, 
                'strand1': row_STRAND1, 
                'strand2': row_STRAND2,
                'ref': row_REF,
                'alt': row_ALT,
                'qual': row_QUAL,
                'svtype': row_SVTYPE
            })
        ###/POS

        ####FORMAT
        format_ = record.FORMAT

        # operation below is limited for manta
        ls_formats = []
        for a_sample in record.samples:
            for a_format in format_.split(':'):
                values = eval('a_sample.data.' + str(a_format))
                if not isinstance(values, list):
                    values = [values]
                for value_idx in range(len(values)):
                    ls_formats.append([row_ID, a_sample.sample, a_format, value_idx, values[value_idx]])
        df_formats_each_record = pd.DataFrame(ls_formats)
        ls_df_formats.append(df_formats_each_record)
        # end of the manta-limited operation

        ####/FORMAT

    ###FILTER
    df_filters = pd.DataFrame(ls_filters)
    ###/FILTER
    
    ###INFO
    ls_df_infos = []
    for info in df_infos_meta.id:
        if len(dict_infos[info]) == 0:
            df_a_info = pd.DataFrame(columns=('id', 'value_idx', info.lower()))
        else:
            df_a_info = pd.DataFrame(dict_infos[info])
        ls_df_infos.append(df_a_info)
    odict_df_infos = OrderedDict([(k, v) for k, v in zip(df_infos_meta.id, ls_df_infos)])
    ###/INFO

    ###POS
    df_pos = pd.DataFrame(ls_pos)
    df_pos = df_pos[['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype']]
    ###/POS

    ###FORMAT
    df_formats = pd.concat(ls_df_formats, ignore_index=True)
    columns = ['id', 'sample', 'format', 'value_idx', 'value']
    df_formats.columns = columns
    ###/FORMAT
   
    args = [df_pos, df_filters, odict_df_infos, df_formats, odict_df_headers, patient_name]
    return Vcf(*args)


def _read_bedpe_empty(df_bedpe, patient_name):
    ls_header = list(df_bedpe.columns)
    ls_header_required = ls_header[:10]
    ls_header_option = ls_header[10:]
    df_svpos = pd.DataFrame(columns=('id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype'))
    df_svlen = pd.DataFrame(columns=('id', 'value_idx', 'svlen'))
    df_svtype = pd.DataFrame(columns=('id', 'value_idx', 'svtype'))
    df_cipos = pd.DataFrame(columns=('id', 'value_idx', 'svtype'))
    df_ciend = pd.DataFrame(columns=('id', 'value_idx', 'svtype'))
    ls_df_infos = []
    for info in ls_header_option:
        df_info = pd.DataFrame(columns=('id', 'value_idx', info))
        ls_df_infos.append(df_info)
    ls_df_infos = [df_svlen, df_svtype, df_cipos, df_ciend] + ls_df_infos   
    ls_infokeys = ['svlen', 'svtype', 'cipos', 'ciend'] + ls_header_option
    odict_df_infos = OrderedDict([(k, v) for k, v in zip(ls_infokeys, ls_df_infos)])
    args = [df_svpos, odict_df_infos, patient_name]
    return Bedpe(*args)
    

def read_bedpe(filepath,
    header_info_path = None,
    svtype_col_name: Optional[str] = None,
    patient_name = None):
    """
    read_bedpe(filepath, header_info_path, svtype_col_name)
    Read a BEDPE file of SV and return Bedpe object.

    Parameters
    ---------------
    filepath: str or file-like object
        Acceptable type is equivalent to that of pandas.read_csv().
    header_info_path:
        Haven't been coded yet.
    svtype_col_name: str or None, default None
        If the bedpe file has a svtype column, please pass the column name to this argument.
    patient_name: str or None, default None
        The patient name
    
    Returns
    ---------------
    A Bedpe object
    """
    if patient_name is None:
        warnings.warn(
            'Passing NoneType to the "patient_name" argument is deprecated.',
            DeprecationWarning
        )
    df_bedpe = pd.read_csv(filepath, sep='\t')
    if df_bedpe.shape[0] == 0:
        return _read_bedpe_empty(df_bedpe, patient_name)
    ls_header = list(df_bedpe.columns)
    ls_header_option = ls_header[10:]
    ls_new_header = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2'] + ls_header_option
    df_bedpe.columns = ls_new_header
    df_bedpe['pos1'] = (df_bedpe['start1'] + df_bedpe['end1']) // 2 + 1
    df_bedpe['pos2'] = (df_bedpe['start2'] + df_bedpe['end2']) // 2 + 1
    df_bedpe['chrom1'] = prepend_chr(df_bedpe['chrom1'])
    df_bedpe['chrom2'] = prepend_chr(df_bedpe['chrom2'])

    if svtype_col_name is None:
        df_bedpe = infer_svtype_from_position(df_bedpe)
        df_svpos = df_bedpe[['name', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'score', 'svtype']].copy()
    else:
        df_svpos = df_bedpe[['name', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'score', svtype_col_name]].rename(columns={svtype_col_name: 'svtype'}).copy()
        df_bedpe['svtype'] = df_svpos['svtype']

    df_svpos = df_svpos.rename(columns={'name': 'id', 'score': 'qual'})
    df_svpos['ref'] = 'N'
    df_svpos = create_alt_field_from_position(df_svpos)
    df_svpos = df_svpos[['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype']]

    ## below: construct INFO tables

    ### svlen table
    def _add_svlen(x):
        x['value_idx'] = 0
        if x.name == 'BND':
            x['svlen'] = 0
            return x
        elif x.name == 'TRA':
            x['svlen'] = 0
            return x
        elif x.name == 'DEL':
            x['svlen'] = x['pos1'] - x['pos2'] + 1
            return x
        elif x.name == 'DUP':
            x['svlen'] = x['pos2'] - x['pos1'] + 1
            return x
        elif x.name == 'INV':
            x['svlen'] = x['pos2'] - x['pos1']
            return x
        else:
            x['svlen'] = abs(x['pos2'] - x['pos1'])
            return x
    df_svlen = df_svpos.groupby('svtype').apply(_add_svlen)[['id', 'value_idx', 'svlen']]

    ### svtype table
    df_svtype = df_svpos[['id', 'svtype']].copy()
    df_svtype['value_idx'] = 0
    df_svtype = df_svtype[['id', 'value_idx', 'svtype']]

    ### cipos and ciend
    df_ci = df_bedpe.copy()
    df_ci[0] = df_ci['start1'] - (df_ci['pos1'] - 1)
    df_ci[1] = df_ci['end1'] - df_ci['pos1']
    df_cipos = df_ci[['name', 0, 1]].rename(columns={'name': 'id'}).set_index('id')
    df_cipos = df_cipos.stack()
    df_cipos = df_cipos.reset_index().rename(columns={'level_1': 'value_idx', 0: 'cipos'})
    
    df_ci = df_bedpe.copy()
    df_ci[0] = df_ci['start2'] - (df_ci['pos2'] - 1)
    df_ci[1] = df_ci['end2'] - df_ci['pos2']
    df_ciend = df_ci[['name', 0, 1]].rename(columns={'name': 'id'}).set_index('id')
    df_ciend = df_ciend.stack()
    df_ciend = df_ciend.reset_index().rename(columns={'level_1': 'value_idx', 0: 'ciend'})

    ls_df_infos = []
    for info in ls_header_option:
        df_info = df_bedpe[['name', info]].rename(columns={'name': 'id'}).copy()
        df_info['value_idx'] = 0
        df_info = df_info[['id', 'value_idx', info]]
        ls_df_infos.append(df_info)
    ls_df_infos = [df_svlen, df_svtype, df_cipos, df_ciend] + ls_df_infos   
    ls_infokeys = ['svlen', 'svtype', 'cipos', 'ciend'] + ls_header_option
    odict_df_infos = OrderedDict([(k, v) for k, v in zip(ls_infokeys, ls_df_infos)])

    args = [df_svpos, odict_df_infos, patient_name]
    return Bedpe(*args)

def prepend_chr(ser):
    """
    prepend_chr(ser)
    Add "chr" to the beginning of every element in the Series.
    If the "chr" prefix has already exist, nothing will be appended but
    the type of Series elements will be changed into str.

    Parameters
    ----------
    ser: pd.Series
        A Series of the chromosome number.
    
    Returns
    -------
    A Series of chromosome numbers in "chr" notation.
    """
    dict_regex = {"^([0-9]+|[XY]|MT)":r"chr\1", r"^(M)":r"chrMT"}
    return ser.astype(str).replace(regex=dict_regex)

def infer_svtype_from_position(position_table):
    df = position_table.copy()
    df['svtype'] = '.'
    mask_bnd = df['chrom1'] != df['chrom2']
    mask_del = (df['chrom1'] == df['chrom2']) & \
                (df['strand1'] == '+') & \
                (df['strand2'] == '-')
    mask_dup = (df['chrom1'] == df['chrom2']) & \
                (df['strand1'] == '-') & \
                (df['strand2'] == '+')
    mask_inv = (df['chrom1'] == df['chrom2']) & \
                (df['strand1'] != '.') & \
                (df['strand1'] == df['strand2'])
    ls_mask = [mask_bnd, mask_del, mask_dup, mask_inv]
    ls_svtype = ['BND', 'DEL', 'DUP', 'INV']
    for mask, svtype in zip(ls_mask, ls_svtype):
        df.loc[mask, 'svtype'] = svtype
    return df

def create_alt_field_from_position(position_table):
    '''
    svtype column is required
    '''
    def _f(x):
        if x.name == '.':
            x['alt'] = '.'
            return x
        elif (x.name != 'BND') & (x.name != 'TRA'):
            x['alt'] = "<{}>".format(x.name)
            return x
        ls_alt = []
        for idx, row in x.iterrows():
            strand1 = row['strand1']
            strand2 = row['strand2']
            chrom2 = row['chrom2']
            pos2 = row['pos2']
            ref = row['ref']
            if (strand1 == '+') & (strand2 == '-'):
                alt = '{0}[{1}:{2}['.format(ref, chrom2, pos2)
            elif (strand1 == '+') & (strand2 == '+'):
                alt = '{0}]{1}:{2}]'.format(ref, chrom2, pos2)
            elif (strand1 == '-') & (strand2 == '+'):
                alt = ']{0}:{1}]{2}'.format(chrom2, pos2, ref)
            else:
                alt = '[{0}:{1}[{2}'.format(chrom2, pos2, ref)
            ls_alt.append(alt)
        x['alt'] = ls_alt
        return x
    df = position_table.copy()
    df = df.groupby('svtype').apply(_f)
    return df

def read_bed(filepath_or_buffer):
    """
    read_bed(filepath_or_buffer)
    Read a BED file into Bed class.
    The input file should conform to the BED format defined by UCSC.
    https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    """
    header=None
    with open(filepath_or_buffer, 'r') as f:
        data = []
        for line in f:
            if line.startswith('track'):
                header = line
            else:
                data.append(line)
    df = pd.read_csv(
        StringIO(''.join(data)),
        sep='\t',
        header=None
    )
    df_columns_origin = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart',
                  'thickEnd', 'itemRgb', 'blockCount', 'blockSize', 'blockStarts']
    df_columns = df_columns_origin[:df.shape[1]]
    df.columns = df_columns

    return Bed(df, header)


    