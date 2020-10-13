import vcf
import pandas as pd
from sgt.core.db import Sgt_core, Sgt_simple
pd.set_option('display.max_columns', 10)
pd.set_option('display.max_colwidth', 30)
pd.set_option('display.width', 1000) 


def read_vcf(filepath, variant_caller="manta"):
    # read vcf files using PyVcf package
    vcf_reader = vcf.Reader(open(filepath, 'r'))

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

    dict_df_headers = {'contigs_meta': df_contigs_meta, 'alts_meta': df_alts_meta, 'infos_meta': df_infos_meta, 'formats_meta': df_formats_meta, 'filters_meta': df_filters_meta, 'samples_meta': df_samples}

    
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
            if isinstance(values, list):
                ls_keys = ['id'] + [info.lower() + '_' + str(i) for i in range(len(values))]
                ls_all_values = [row_ID] + [v for v in values]
                dict_a_info = {k: v for k, v in zip(ls_keys, ls_all_values)}
            else:
                column_name = info.lower() + '_0'
                dict_a_info = {'id': row_ID, column_name: values}
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
        elif row_INFO['SVTYPE'] == 'DEL':
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
            row_STRAND1 = '+'
            row_STRAND2 = '-'
        elif row_INFO['SVTYPE'] == 'DUP':
            row_CHROM2 = row_CHROM1
            row_POS2 = row_INFO['END']
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
            df_a_info = pd.DataFrame(columns=('id', info + '_0'))
        else:
            df_a_info = pd.DataFrame(dict_infos[info])
        ls_df_infos.append(df_a_info)
    dict_df_infos = {k: v for k, v in zip(df_infos_meta.id, ls_df_infos)}
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
   
    args = [df_pos, df_filters, dict_df_infos, df_formats, dict_df_headers]
    return Sgt_core(*args)

def _read_bedpe_empty(df_bedpe):
    ls_header = list(df_bedpe.columns)
    ls_header_required = ls_header[:10]
    ls_header_option = ls_header[10:]
    df_svpos = pd.DataFrame(columns=('id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'ref', 'alt', 'qual', 'svtype'))
    df_svlen = pd.DataFrame(columns=('id', 'svlen_0'))
    df_svtype = pd.DataFrame(columns=('id', 'svtype_0'))
    ls_df_infos = []
    for info in ls_header_option:
        df_info = pd.DataFrame(columns=('id', info + '_0'))
        ls_df_infos.append(df_info)
    ls_df_infos = [df_svlen, df_svtype] + ls_df_infos   
    ls_infokeys = ['svlen', 'svtype'] + ls_header_option
    dict_df_infos = {k: v for k, v in zip(ls_infokeys, ls_df_infos)}
    args = [df_svpos, dict_df_infos]
    return Sgt_simple(*args)
    

def read_bedpe(filepath, header_info_path=None, svtype_col_name=''):
    df_bedpe = pd.read_csv(filepath, sep='\t')
    if df_bedpe.shape[0] == 0:
        return _read_bedpe_empty(df_bedpe)
    ls_header = list(df_bedpe.columns)
    ls_header_required = ls_header[:10]
    ls_header_option = ls_header[10:]
    ls_new_header = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2'] + ls_header_option
    df_bedpe.columns = ls_new_header
    df_bedpe['pos1'] = (df_bedpe['start1'] + df_bedpe['end1']) // 2
    df_bedpe['pos2'] = (df_bedpe['start2'] + df_bedpe['end2']) // 2

    if svtype_col_name == '':
        df_bedpe = infer_svtype_from_position(df_bedpe)
        df_svpos = df_bedpe[['name', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'score', 'svtype']].copy()
    else:
        df_svpos = df_bedpe[['name', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2', 'score', svtype_col_name]].rename(columns={svtype_col_name: 'svtype'}).copy()
        df_bedpe['svtype'] = df_svpos['svtype']

    df_svpos = df_svpos.rename(columns={'name': 'id', 'score': 'qual'})
    df_svpos['ref'] = 'N'
    df_svpos['alt'] = '.'

    def _add_svlen(x):
        if x.name == 'BND':
            x['svlen_0'] = 0
            return x
        elif x.name == 'TRA':
            x['svlen_0'] = 0
            return x
        elif x.name == 'DEL':
            x['svlen_0'] = x['pos1'] - x['pos2']
            return x
        elif x.name == 'DUP':
            x['svlen_0'] = x['pos2'] - x['pos1']
            return x
        elif x.name == 'INV':
            x['svlen_0'] = x['pos2'] - x['pos1']
            return x
        else:
            x['svlen_0'] = 0
            return x

    df_bedpe = df_bedpe.groupby('svtype').apply(_add_svlen)

    df_bedpe['cipos_0'] = df_bedpe['start1'] - df_bedpe['pos1']
    df_bedpe['cipos_1'] = df_bedpe['end1'] - df_bedpe['pos1'] - 1
    df_bedpe['ciend_0'] = df_bedpe['start2'] - df_bedpe['pos2']
    df_bedpe['ciend_1'] = df_bedpe['end2'] - df_bedpe['pos2'] - 1

    df_svlen = df_bedpe[['name', 'svlen_0']].rename(columns={'name': 'id'})
    df_svtype = df_svpos[['id', 'svtype']].rename(columns={'svtype': 'svtype_0'})
    df_cipos = df_bedpe[['name', 'cipos_0', 'cipos_1']].rename(columns={'name': 'id'})
    df_ciend = df_bedpe[['name', 'ciend_0', 'ciend_1']].rename(columns={'name': 'id'})

    ls_df_infos = []
    for info in ls_header_option:
        df_info = df_bedpe[['name', info]].rename(columns={info: info + '_0', 'name': 'id'})
        ls_df_infos.append(df_info)
    ls_df_infos = [df_svlen, df_svtype, df_cipos, df_ciend] + ls_df_infos   
    ls_infokeys = ['svlen', 'svtype', 'cipos', 'ciend'] + ls_header_option
    dict_df_infos = {k: v for k, v in zip(ls_infokeys, ls_df_infos)}

    args = [df_svpos, dict_df_infos]
    return Sgt_simple(*args)

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
            return x
        elif (x.name != 'BND') & (x.name != 'TRA'):
            x['alt'] = "<{}>".format(x.name)
            return x
        ls_out = []
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
            row['alt'] = alt
            ls_out.append(row)
        df_out = pd.DataFrame(ls_out)
        return df_out
    df = position_table.copy()
    df = df.groupby('svtype').apply(_f)
    return df
    



