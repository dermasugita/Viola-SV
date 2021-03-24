import click
import viola
from io import StringIO, TextIOWrapper
import os

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version='1.0.0')
@click.option('--input-dir', default=None, help='The directory of input files. When specified, the --files argument is disabled.')
@click.option('--input-files', default=None, help='The input files separeted by comma. When specified, the --dir argument is disabled.')
@click.option('--input-files-id', default=None, help='The sample ID of input files separeted by comma.')
@click.option('--format', 'format_', default='vcf', type=click.Choice(['vcf', 'bedpe']), help='File format. vcf or bedpe.')
@click.option('--caller', default='manta', type=click.Choice(['manta', 'delly', 'lumpy', 'gridss']), help='The name of SV caller by which the input VCF was generated. This option can be specified when --format=vcf.')
@click.option('--svtype-col-name', default=None, help='Name of the column of BEDPE files that indicate SV type. If not specified, SV type will be infered. This option can be specified when --format=bedpe')
@click.option('--as-breakpoint', is_flag=True, help='Convert SVTYPE=BND records into breakpoint-wise SV records and infer its SVTYPE. This option is used when --format=vcf')
@click.option('--definitions', default=None, help='Path to the definition file of custom SV class.')
@click.argument('output', type=click.File('w'))

def generate_feature_matrix(input_dir, input_files, input_files_id, format_, caller, svtype_col_name, as_breakpoint, definitions, output):
    """
    Generate feature matrix from VCF or BEDPE files.
    """
    if format_ == 'bedpe':
        if (input_dir is None) & (input_files is None):
            return
        elif (input_files is None):
            data = viola.read_bedpe_multi(input_dir, svtype_col_name=svtype_col_name)
        elif (input_dir is None):
            ls_input = input_files.split(',')
            ls_bedpe = [viola.read_bedpe(path, svtype_col_name=svtype_col_name) for path in ls_input]
            if input_files_id is None:
                ls_names = range(len(ls_bedpe))
            else:
                ls_names = input_files_id.split(',')
            data = viola.MultiBedpe(ls_bedpe, ls_names)
        else:
            return
    else:
        if (input_dir is None) & (input_files is None):
            return
        elif (input_files is None):
            data = viola.read_vcf_multi(input_dir, variant_caller=caller, as_breakpoint=as_breakpoint)
        elif (input_dir is None):
            ls_input = input_files.split(',')
            if as_breakpoint:
                ls_vcf = [viola.read_vcf(path, variant_caller=caller).breakend2breakpoint() for path in ls_input]
            else:
                ls_vcf = [viola.read_vcf(path, variant_caller=caller) for path in ls_input]

            if input_files_id is None:
                ls_names = range(len(ls_vcf))
            else:
                ls_names = input_files_id.split(',')
            data = viola.MultiBedpe(ls_vcf, ls_names)
        else:
            return

    result = data.classify_manual_svtype(definitions=definitions)
    result.to_csv(output, sep='\t')
    