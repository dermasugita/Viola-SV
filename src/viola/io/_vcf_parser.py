import vcf
import pandas as pd
import numpy as np
from collections import OrderedDict
from viola.core.vcf import Vcf

################################################
# Manta
################################################

def read_vcf_manta(vcf_reader, patient_name):
    metadata = vcf_reader.metadata
    metadata['variantcaller'] = 'manta'
    # obtain header informations
    odict_contigs = vcf_reader.contigs
    df_contigs_meta = pd.DataFrame(odict_contigs, index=('id', 'length')).T.reset_index(drop=True)
    df_contigs_meta['length'] = df_contigs_meta['length'].astype(int)

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
        row_SVTYPE = row_INFO['SVTYPE']
        for info in df_infos_meta.id:
            values = row_INFO.get(info, 'none')
            if values == 'none':
                if info == 'CIPOS':
                    dict_infos['CIPOS'].append({'id': row_ID, 'value_idx': 0, 'cipos': 0})
                    dict_infos['CIPOS'].append({'id': row_ID, 'value_idx': 1, 'cipos': 0})
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
            if row_POS2 is None:
                row_STRAND2 = None
            else:
                row_STRAND2 = '-' if row_ALT.remoteOrientation else '+'
        elif row_SVTYPE == 'INV': # manta-specific operation
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
            # Manta
            if row_INFO.get('INV3', False):
                row_STRANDs = '+'
            else:
                row_POS1 += 1
                row_POS2 += 1
                row_STRANDs = '-'
            # /Manta
            row_STRAND1, row_STRAND2 = row_STRANDs, row_STRANDs
        elif row_SVTYPE == 'DEL':
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END'] + 1
            row_STRAND1 = '+'
            row_STRAND2 = '-'
        elif row_SVTYPE == 'DUP':
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

        row_SVTYPE = row_SVTYPE
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
    df_filters = pd.DataFrame(ls_filters, columns=['id', 'filter'])
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

    #### Generate CIEND table by merging MATEID and CIPOS
    if 'MATEID' in odict_df_infos:
        df_mateid = odict_df_infos['MATEID']
        df_cipos = odict_df_infos['CIPOS']
        df_merged = pd.merge(df_cipos, df_mateid, on='id')
        df_merged = df_merged[['mateid', 'value_idx_x', 'cipos']]
        df_merged.columns = ['id', 'value_idx', 'ciend']
        df_ciend = odict_df_infos['CIEND']
        odict_df_infos['CIEND'] = pd.concat([df_ciend, df_merged], ignore_index=True)
    #### /Generate CIEND table by merging MATEID and CIPOS

    ###/INFO

    ###POS
    df_pos = pd.DataFrame(ls_pos)
    if df_pos.empty:
        df_pos = pd.DataFrame([], columns=['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype'])
    else:
        df_pos = df_pos[['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype']]
    ###/POS

    ###FORMAT
    if len(ls_df_formats) == 0:
        df_formats = pd.DataFrame([], columns=['id', 'sample', 'format', 'value_idx', 'value'])
    else:
        df_formats = pd.concat(ls_df_formats, ignore_index=True)
        columns = ['id', 'sample', 'format', 'value_idx', 'value']
        df_formats.columns = columns
    ###/FORMAT
   
    args = [df_pos, df_filters, odict_df_infos, df_formats, odict_df_headers, metadata, patient_name]
    return Vcf(*args)

################################################
# /Manta
################################################

################################################
# Delly
################################################

def read_vcf_delly(vcf_reader, patient_name):
    """
    Replace table name SVLEN into SVLENORG to avoid name conflict.
    """
    metadata = vcf_reader.metadata
    metadata['variantcaller'] = 'delly'
    # obtain header informations
    odict_contigs = vcf_reader.contigs
    df_contigs_meta = pd.DataFrame(odict_contigs, index=('id', 'length')).T.reset_index(drop=True)
    df_contigs_meta['length'] = df_contigs_meta['length'].astype(int)

    # obtain alteration class informations
    odict_alts = vcf_reader.alts
    df_alts_meta = pd.DataFrame(odict_alts, index=('id', 'description')).T.reset_index(drop=True)

    # obtain info field information
    odict_infos = vcf_reader.infos
    df_infos_meta = pd.DataFrame(odict_infos, index=('id', 'number', 'type', 'description', 'source', 'version')).T.reset_index(drop=True)
    ## Delly
    df_infos_meta.loc[df_infos_meta['id'] == 'SVLEN', 'id'] = 'SVLENORG'
    df_infos_meta = pd.concat([df_infos_meta, pd.DataFrame({'id':['SVLEN'], 'number':[None], 'type':['Integer'], 'description':['Length of SV'], 'source':[None], 'version':[None]})], ignore_index=True)
    ## /Delly

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
        row_SVTYPE = row_INFO['SVTYPE']

        ##### delly: get svlen
        if row_SVTYPE == 'BND':
            svlen = 0
        elif row_SVTYPE == 'TRA':
            svlen = 0
        elif row_SVTYPE == 'DEL':
            svlen = row_POS1 - row_INFO['END']
        elif row_SVTYPE == 'DUP':
            svlen = row_INFO['END'] - row_POS1
        elif row_SVTYPE == 'INV':
            svlen = row_INFO['END'] - row_POS1
        else:
            svlen = 0
        ##### /delly: get svlen

        for info in df_infos_meta.id:
            if info == 'SVLENORG':
                values = row_INFO.get('SVLEN', 'none')
            elif info == 'SVLEN':
                values = svlen
            else:
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
            if row_POS2 is None:
                row_STRAND2 = None
            else:
                row_STRAND2 = '-' if row_ALT.remoteOrientation else '+'
        elif row_SVTYPE == 'INV': # delly-specific operation
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
            # delly
            if row_INFO['CT'] == '3to3':
                row_STRANDs = '+'
            else:
                row_STRANDs = '-'
            # /delly
            row_STRAND1, row_STRAND2 = row_STRANDs, row_STRANDs
        elif row_SVTYPE == 'DEL':
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END'] + 1
            row_STRAND1 = '+'
            row_STRAND2 = '-'
        elif row_SVTYPE == 'DUP':
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
    df_filters = pd.DataFrame(ls_filters, columns=['id', 'filter'])
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
    if df_pos.empty:
        df_pos = pd.DataFrame([], columns=['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype'])
    else:
        df_pos = df_pos[['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype']]
    ###/POS

    ###FORMAT
    if len(ls_df_formats) == 0:
        df_formats = pd.DataFrame([], columns=['id', 'sample', 'format', 'value_idx', 'value'])
    else:
        df_formats = pd.concat(ls_df_formats, ignore_index=True)
        columns = ['id', 'sample', 'format', 'value_idx', 'value']
        df_formats.columns = columns
    ###/FORMAT
   
    args = [df_pos, df_filters, odict_df_infos, df_formats, odict_df_headers, metadata, patient_name]
    return Vcf(*args)

################################################
# /Delly
################################################

################################################
# Lumpy
################################################

def read_vcf_lumpy(vcf_reader, patient_name):
    """
    split INV lines.
    INFO=SU is modified according to the INFO=STRANDS.
    """
    metadata = vcf_reader.metadata
    metadata['variantcaller'] = 'lumpy'
    # obtain header informations
    ################# This is empty DataFrame! ############################
    odict_contigs = vcf_reader.contigs
    df_contigs_meta = pd.DataFrame(odict_contigs, index=('id', 'length')).T.reset_index(drop=True)
    df_contigs_meta['length'] = df_contigs_meta['length'].astype(int)

    # obtain alteration class informations
    odict_alts = vcf_reader.alts
    df_alts_meta = pd.DataFrame(odict_alts, index=('id', 'description')).T.reset_index(drop=True)

    # obtain info field information
    odict_infos = vcf_reader.infos
    df_infos_meta = pd.DataFrame(odict_infos, index=('id', 'number', 'type', 'description', 'source', 'version')).T.reset_index(drop=True)
    df_suorg = pd.DataFrame([['SUORG', None, 'Integer', 'Original value of Lumpy SU', None, None]], columns = ('id', 'number', 'type', 'description', 'source', 'version'))
    df_infos_meta = pd.concat([df_infos_meta, df_suorg], ignore_index=True)

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
        row_SVTYPE = record.INFO['SVTYPE']
        if row_SVTYPE == 'INV':
            is_inv = True
        else:
            is_inv = False
        
        if is_inv:
            row_ID = record.ID
            row_ID1 = str(record.ID) + '_1'
            row_ID2 = str(record.ID) + '_2'
        else:
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
            if is_inv:
                ls_filters.append({'id': row_ID1, 'filter': filter_})
                ls_filters.append({'id': row_ID2, 'filter': filter_})
            else:
                ls_filters.append({'id': row_ID, 'filter': filter_})

        #####/FILTER 

        #####INFO
        row_INFO = record.INFO
        for info in df_infos_meta.id:
            values = row_INFO.get(info, 'none')
            if values == 'none':
                if info == 'EVENT' and is_inv:
                    values = row_ID
                else:
                    continue
            if not isinstance(values, list):
                values = [values]
            ls_keys = ['id', 'value_idx', info.lower()] 
            if is_inv:
                for idx, value in enumerate(values):
                    if info == 'SU':
                        value = row_INFO['STRANDS'][0].split(':')[1]
                        value = int(value)
                    elif info == 'STRANDS':
                        if idx == 1: continue
                    ls_value1 = [row_ID1, idx, value]
                    dict_a_info1 = {k: v for k, v in zip(ls_keys, ls_value1)}
                    dict_infos[info].append(dict_a_info1)
                for idx, value in enumerate(values):
                    if info == 'SU':
                        value = row_INFO['STRANDS'][1].split(':')[1]
                        value = int(value)
                    elif info == 'STRANDS':
                        if idx == 0: continue
                        idx = 0
                    ls_value2 = [row_ID2, idx, value]
                    dict_a_info2 = {k: v for k, v in zip(ls_keys, ls_value2)}
                    dict_infos[info].append(dict_a_info2)
            else:
                for idx, value in enumerate(values):
                    ls_value = [row_ID, idx, value]
                    dict_a_info = {k: v for k, v in zip(ls_keys, ls_value)}
                    dict_infos[info].append(dict_a_info)
        
        if is_inv:
            suorg = row_INFO['SU'][0]
            dict_infos['SUORG'].append({'id': row_ID1, 'value_idx': 0, 'suorg': suorg})
            dict_infos['SUORG'].append({'id': row_ID2, 'value_idx': 0, 'suorg': suorg})
        #####/INFO

        ###POS
        if isinstance(row_ALT, vcf.model._Breakend):
            row_CHROM2 = row_ALT.chr
            row_POS2 = row_ALT.pos
            row_STRAND1 = '-' if row_ALT.orientation else '+'
            if row_POS2 is None:
                row_STRAND2 = None
            else:
                row_STRAND2 = '-' if row_ALT.remoteOrientation else '+'
        elif row_SVTYPE == 'INV': # manta-specific operation
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
        elif row_SVTYPE == 'DEL':
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END'] + 1
            row_STRAND1 = '+'
            row_STRAND2 = '-'
        elif row_SVTYPE == 'DUP':
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

        row_SVTYPE = row_SVTYPE
        if is_inv:
            ls_pos.append({'id': row_ID1, 'chrom1': row_CHROM1, 'pos1': row_POS1, 
                    'chrom2': row_CHROM2, 'pos2': row_POS2, 'strand1': '+', 
                    'strand2': '+', 'ref': row_REF, 'alt': row_ALT,
                    'qual': row_QUAL, 'svtype': row_SVTYPE
                })
            ls_pos.append({'id': row_ID2, 'chrom1': row_CHROM1, 'pos1': row_POS1+1, 
                    'chrom2': row_CHROM2, 'pos2': row_POS2+1, 'strand1': '-', 
                    'strand2': '-', 'ref': row_REF, 'alt': row_ALT,
                    'qual': row_QUAL, 'svtype': row_SVTYPE
                })
        else:
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
                    if is_inv:
                        ls_formats.append([row_ID1, a_sample.sample, a_format, value_idx, values[value_idx]])
                        ls_formats.append([row_ID2, a_sample.sample, a_format, value_idx, values[value_idx]])
                    else:
                        ls_formats.append([row_ID, a_sample.sample, a_format, value_idx, values[value_idx]])
        df_formats_each_record = pd.DataFrame(ls_formats)
        ls_df_formats.append(df_formats_each_record)
        # end of the manta-limited operation

        ####/FORMAT

    ###FILTER
    df_filters = pd.DataFrame(ls_filters, columns=['id', 'filter'])
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
    if df_pos.empty:
        df_pos = pd.DataFrame([], columns=['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype'])
    else:
        df_pos = df_pos[['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype']]
    ###/POS

    ###FORMAT
    if len(ls_df_formats) == 0:
        df_formats = pd.DataFrame([], columns=['id', 'sample', 'format', 'value_idx', 'value'])
    else:
        df_formats = pd.concat(ls_df_formats, ignore_index=True)
        columns = ['id', 'sample', 'format', 'value_idx', 'value']
        df_formats.columns = columns
    ###/FORMAT

    ############## Generate contigs_meta ####################
    arr_contigs = np.append(df_pos['chrom1'].values,  df_pos['chrom2'].values)
    arr_contigs = np.unique(arr_contigs)
    arr_contigs_length = np.repeat('.', len(arr_contigs))
    df_contigs_meta = pd.DataFrame({'id': arr_contigs, 'length': arr_contigs_length})
    odict_df_headers['contigs_meta'] = df_contigs_meta
   
    args = [df_pos, df_filters, odict_df_infos, df_formats, odict_df_headers, metadata, patient_name]
    return Vcf(*args)

################################################
# /Lumpy
################################################

################################################
# Gridss
################################################

def read_vcf_gridss(vcf_reader, patient_name):
    metadata = vcf_reader.metadata
    metadata['variantcaller'] = 'gridss'
    # obtain header informations
    odict_contigs = vcf_reader.contigs
    df_contigs_meta = pd.DataFrame(odict_contigs, index=('id', 'length')).T.reset_index(drop=True)
    df_contigs_meta['length'] = df_contigs_meta['length'].astype(int)

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
        row_SVTYPE = row_INFO['SVTYPE']
        for info in df_infos_meta.id:
            values = row_INFO.get(info, 'none')
            if values == 'none':
                if info == 'CIPOS':
                    dict_infos['CIPOS'].append({'id': row_ID, 'value_idx': 0, 'cipos': 0})
                    dict_infos['CIPOS'].append({'id': row_ID, 'value_idx': 1, 'cipos': 0})
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
            if row_POS2 is None:
                row_STRAND2 = None
            else:
                row_STRAND2 = '-' if row_ALT.remoteOrientation else '+'
        elif row_SVTYPE == 'INV': # manta-specific operation
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
            # Manta
            if row_INFO.get('INV3', False):
                row_STRANDs = '+'
            else:
                row_POS1 += 1
                row_POS2 += 1
                row_STRANDs = '-'
            # /Manta
            row_STRAND1, row_STRAND2 = row_STRANDs, row_STRANDs
        elif row_SVTYPE == 'DEL':
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END'] + 1
            row_STRAND1 = '+'
            row_STRAND2 = '-'
        elif row_SVTYPE == 'DUP':
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

        row_SVTYPE = row_SVTYPE
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
    df_filters = pd.DataFrame(ls_filters, columns=['id', 'filter'])
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

    #### Generate CIEND table by merging MATEID and CIPOS
    if 'MATEID' in odict_df_infos:
        df_mateid = odict_df_infos['MATEID']
        df_cipos = odict_df_infos['CIPOS']
        df_merged = pd.merge(df_cipos, df_mateid, on='id')
        df_merged = df_merged[['mateid', 'value_idx_x', 'cipos']]
        df_merged.columns = ['id', 'value_idx', 'ciend']
        odict_df_infos['CIEND'] = df_merged
    #### /Generate CIEND table by merging MATEID and CIPOS
    
    ###/INFO

    ###POS
    df_pos = pd.DataFrame(ls_pos)
    if df_pos.empty:
        df_pos = pd.DataFrame([], columns=['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype'])
    else:
        df_pos = df_pos[['id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype']]
    ###/POS

    ###FORMAT
    if len(ls_df_formats) == 0:
        df_formats = pd.DataFrame([], columns=['id', 'sample', 'format', 'value_idx', 'value'])
    else:
        df_formats = pd.concat(ls_df_formats, ignore_index=True)
        columns = ['id', 'sample', 'format', 'value_idx', 'value']
        df_formats.columns = columns
    ###/FORMAT
   
    args = [df_pos, df_filters, odict_df_infos, df_formats, odict_df_headers, metadata, patient_name]
    return Vcf(*args)

################################################
# /Gridss
################################################
