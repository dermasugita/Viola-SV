import os
from sgt.io.parser import read_bedpe
from sgt.core.cohort import MultiBedpe

def read_bedpe_multi(dir_path: str,
    svtype_col_name: str = 'svclass',
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
        ls_bedpe.append(read_bedpe(abspath, svtype_col_name=svtype_col_name))
        patient_id = f.replace('.bedpe', '')
        ls_names.append(patient_id)
    multi_bedpe = MultiBedpe(ls_bedpe, ls_names)
    return multi_bedpe