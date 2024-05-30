import viola
import click
import sys
from io import StringIO, TextIOWrapper

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version='1.0.0')
@click.option('--caller', default='manta', help='The name of SV caller by which the input VCF was generated.\n[manta, delly, lumpy, gridss] could be acceptable (default, manta).')
@click.option('-i', '--info', help='The names of INFO fields to return. To specify multiple INFO, separate them by commas. ex. --info SVTYPE,SVLEN,END')
@click.option('-f','--filter', 'filter_', is_flag=True, help='If specified, FILTER field of the VCF files is included in output BEDPE.')
@click.option('-m', '--format', 'format_', is_flag=True, help='If specified, FORMAT field of the VCF files is included in output BEDPE.')
@click.argument('vcf', default=sys.stdin, type=click.File('r'))
def vcf2bedpe(caller, info, filter_, format_, vcf):
   """
   Convert a VCF file into a BEDPE file.

   A VCF argument is the path to the input VCF file.
   """
   if isinstance(vcf, TextIOWrapper):
      vcf_obj = viola.read_vcf(StringIO(vcf.read()), variant_caller=caller)
   else:
      vcf_obj = viola.read_vcf(vcf, variant_caller=caller)
   if info is not None:
      ls_info = info.split(',')
      ls_info_lower = [i.lower() for i in ls_info]
   else:
      ls_info_lower = []
   bedpe_like = vcf_obj.to_bedpe_like(custom_infonames=ls_info_lower, add_filters=filter_, add_formats=format_)
   click.echo(bedpe_like.to_string(index=None))
