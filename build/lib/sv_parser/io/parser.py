import vcf
import pandas as pd
import sqlite3
pd.set_option('display.max_columns', 10)
pd.set_option('display.max_colwidth', 30)
pd.set_option('display.width', 1000) 

def read_vcf(filepath, variant_caller="manta"):
    # read vcf files using PyVcf package
    vcf_reader = vcf.Reader(open(filepath, 'r'))

    # obtain header informations
    odict_contigs = vcf_reader.contigs
    df_contigs_meta = pd.DataFrame(odict_contigs, index=('ID', 'Length')).T.reset_index(drop=True)

    # obtain alteration class informations
    odict_alts = vcf_reader.alts
    df_alts_meta = pd.DataFrame(odict_alts, index=('ID', 'Description')).T.reset_index(drop=True)

    # obtain info field information
    odict_infos = vcf_reader.infos
    df_infos_meta = pd.DataFrame(odict_infos, index=('ID', 'Number', 'Type', 'Description', 'Source', 'Version')).T.reset_index(drop=True)

    # obtain FORMAT informations
    odict_formats = vcf_reader.formats
    df_formats_meta = pd.DataFrame(odict_formats, index=('ID', 'Number', 'Type','Description')).T.reset_index(drop=True)

    # obtain FILTER informations
    odict_filters = vcf_reader.filters
    df_filters = pd.DataFrame(odict_filters, index=('ID', 'Description')).T.reset_index(drop=True)

    ls_df_headers = [df_contigs_meta, df_alts_meta, df_infos_meta, df_formats_meta, df_filters]

    
    ls_pos = []
    ls_filters = []
    dict_infos = {k: [] for k in df_infos_meta.ID}
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
        if len(ls_filter) == 0:
            ls_filter = ['PASS'] 
        row_FILTER = ls_filter
        for filter_ in row_FILTER:
            ls_filters.append({'ID': row_ID, 'FILTER': filter_})
        #####/FILTER 

        #####INFO
        row_INFO = record.INFO
        for info in df_infos_meta.ID:
            values = row_INFO.get(info, 'none')
            if values == 'none':
                continue
            if isinstance(values, list):
                ls_keys = ['ID', 'INFO'] + ['VALUE' + str(i+1) for i in range(len(values))]
                ls_all_values = [row_ID, info] + [v for v in values]
                dict_a_info = {k: v for k, v in zip(ls_keys, ls_all_values)}
            else:
                dict_a_info = {'ID': row_ID, 'INFO': info, 'VALUE1': values}
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
            row_STRANDs = '-' if row_INFO.get('INV3', False) else '+'
            row_STRAND1, row_STRAND2 = row_STRANDs, row_STRANDs
        else:
            row_CHROM2 = row_CHROM2
            row_POS2 = row_INFO['END']
            row_STRAND1 = '+'
            row_STRAND2 = '-'
        row_SVCLASS = row_INFO['SVTYPE']
        ls_pos.append({
                'ID': row_ID, 
                'CHROM1': row_CHROM1, 
                'POS1': row_POS1, 
                'CHROM2': row_CHROM2, 
                'POS2': row_POS2, 
                'STRAND1': row_STRAND1, 
                'STRAND2': row_STRAND2,
                'REF': row_REF,
                'ALT': row_ALT,
                'QUAL': row_QUAL,
                'SVCLASS': row_SVCLASS
            })
        ###/POS

        


        ####FORMAT
        format_ = record.FORMAT

        # operation below is limited for manta
        ls_formats = []
        for a_sample in record.samples:
            for a_format in format_.split(':'):
                values = eval('a_sample.data.' + str(a_format))
                ls_formats.append([row_ID, a_sample.sample, a_format] + values)
        df_formats_each_record = pd.DataFrame(ls_formats)
        ls_df_formats.append(df_formats_each_record)
        # end of the manta-limited operation

        ####/FORMAT

        

    ###FILTER
    df_filters = pd.DataFrame(ls_filters)
    ###/FILTER
    
    ###INFO
    ls_df_infos = []
    for info in df_infos_meta.ID:
        df_a_info = pd.DataFrame(dict_infos[info])
        ls_df_infos.append(df_a_info)
    dict_df_infos = {k: v for k, v in zip(df_infos_meta.ID, ls_df_infos)}
    ###/INFO

    ###POS
    df_pos = pd.DataFrame(ls_pos)
    ###/POS

    ###FORMAT
    df_formats = pd.concat(ls_df_formats, ignore_index=True)
    columns = ['ID', 'SAMPLE', 'FORMAT'] + ['VALUE' + str(i+1) for i in range(df_formats.shape[1] - 3)]
    df_formats.columns = columns
    ###/FORMAT
    
    return([df_pos, df_filters, dict_df_infos, df_formats, ls_df_headers])


