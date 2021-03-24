import click
from viola.cli.vcf2bedpe import vcf2bedpe
from viola.cli.matrix_generator import generate_feature_matrix
from viola.cli.signature_extractor import extract_signature
@click.group()
def viola():
   pass

viola.add_command(vcf2bedpe)
viola.add_command(generate_feature_matrix)
viola.add_command(extract_signature)

if __name__ == '__main__':
   viola()