import os
from viola.io.parser import read_bedpe, read_vcf
from viola.core.cohort import MultiBedpe, MultiVcf

def read_vcf_multi(dir_path: str,
    variant_caller: str = 'manta',
    as_breakpoint: bool = False,
    exclude_empty_cases: bool = False,
    file_extension: str = 'vcf',
    escape_dot_files: bool = True):
    """
    read_vcf_multi(dir_path, variant_caller, as_breakpoint, file_extension, escape_dot_files)
    Read VCF files in a specified directory at the same time and return as MultiBedpe object.
    
    Parameters
    ----------
    dir_path: str
        Path to the directory containing the VCF files.
    variant_caller: str
        Let this function know which SV caller was used to create vcf file.
        Only "manta" is supported for now.
    as_breakpoint: bool, default False
        Convert SVTYPE=BND record into breakpoint-wise SV. The SVTYPE is predicted and can be DEL, DUP, INV, INS, TRA or BND.
    exclude_empty_cases: bool, default False
        If True, skip reading empty VCF files.
    file_extension: str or None, default 'vcf'
        File extension of BEDPE files. If you want to load files with no extension, specify None.
    escape_dot_files: bool, default True
        If True, avoid reading hidden files in the directory.
    """
    ls_vcf = []
    ls_names = []
    for f in os.listdir(dir_path):
        if escape_dot_files and f.startswith('.'): continue
        if (file_extension is not None) and (f.split('.')[-1] != file_extension): continue
        abspath = os.path.abspath(os.path.join(dir_path, f))
        vcf = read_vcf(abspath, variant_caller=variant_caller)
        if exclude_empty_cases & (vcf.sv_count == 0):
            continue
        if as_breakpoint:
            vcf = vcf.breakend2breakpoint()
        ls_vcf.append(vcf)
        patient_id = f.replace('.vcf', '')
        ls_names.append(patient_id)
    multi_vcf = MultiVcf(ls_vcf, ls_names)
    return multi_vcf



def read_bedpe_multi(dir_path: str,
    svtype_col_name: str = None,
    exclude_empty_cases: bool = False,
    file_extension: str = 'bedpe',
    escape_dot_files: bool = True):
    """
    read_bedpe_multi(dir_path, svtype_col_name, file_extension, escape_dot_files)
    Read BEDPE files in a specified directory at the same time and return as MultiBedpe object.
    
    Parameters
    ----------
    dir_path: str
        Path to the directory containing the BEDPE files.
    svtype_col_name: str, default 'svclass'
        If the bedpe file has a svtype column, please pass the column name to this argument.
    exclude_empty_cases: bool, default False
        If True, skip reading empty BEDPE files.
    file_extension: str or None, default 'bedpe'
        File extension of BEDPE files. If you want to load files with no extension, specify None.
    escape_dot_files: bool, default True
        If True, avoid reading hidden files in the directory.
    """
    ls_bedpe = []
    ls_names = []
    for f in os.listdir(dir_path):
        if escape_dot_files and f.startswith('.'): continue
        if (file_extension is not None) and (f.split('.')[-1] != file_extension): continue
        abspath = os.path.abspath(os.path.join(dir_path, f))
        bedpe = read_bedpe(abspath, svtype_col_name=svtype_col_name)
        if exclude_empty_cases & (bedpe.sv_count == 0): continue
        ls_bedpe.append(bedpe)
        patient_id = f.replace('.bedpe', '')
        ls_names.append(patient_id)
    multi_bedpe = MultiBedpe(ls_bedpe, ls_names)
    return multi_bedpe