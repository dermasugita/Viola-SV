import click
from viola.cli.vcf2bedpe import vcf2bedpe
@click.group()
def viola():
   pass

viola.add_command(vcf2bedpe)

if __name__ == '__main__':
   viola()