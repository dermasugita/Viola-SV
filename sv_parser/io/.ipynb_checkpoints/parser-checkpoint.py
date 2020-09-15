import vcf

def read_vcf(
    filepath,
    variant_caller="manta"
):
    vcf_reader = vcf.Reader(open(filepath, 'r'))
    vcf_infos = vcf_reader.infos
    ls_info_keys = ['INFO.' + k for k in vcf_infos.keys()]
    ls_format = [a_sample + '.' + a_format for a_sample in vcf_reader.samples for a_format in vcf_reader.formats.keys()]
    ls_columns = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER'] + ls_info_keys + ls_format
    
    ls_all = []
    
    for record in vcf_reader:
        dict_row = {}
        # ex: Record(CHROM=chr1, POS=9163435, REF=G, ALT=[<DUP:TANDEM>])
        dict_row['CHROM'] = record.CHROM
        dict_row['POS'] = record.POS
        dict_row['REF'] = record.REF
        dict_row['ALT'] = record.ALT[0]
        dict_row['QUAL'] = record.QUAL
         
        ls_filter = record.FILTER
        if len(ls_filter) == 0:
            ls_filter = ['PASS']
        dict_row['FILTER'] = ls_filter
        
        dict_info = record.INFO
        dict_info = {'INFO.' + k: v for k, v in dict_info.items()}
        format_ = record.FORMAT
        
        dict_format = {}
        for a_sample in record.samples:
            for a_format in format_.split(':'):
                key = a_sample.sample + '.' + str(a_format)
                value = eval('a_sample.data.' + str(a_format))
                dict_format[key] = value
        
        dict_row = dict(dict_row, **dict_info, **dict_format)
        
        ls_all.append(dict_row)
        
    return(pd.DataFrame(ls_all, columns=ls_columns))