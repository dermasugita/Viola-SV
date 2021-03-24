from viola.cli.viola import viola
from click.testing import CliRunner
import sys, os
HERE = os.path.abspath(os.path.dirname(__file__))
bedpe_path = os.path.join(HERE, 'data/bedpe')
vcf_path = os.path.join(HERE, 'data/vcf')
bedpe_definition_path = os.path.join(HERE, 'data/example_definition_bedpe.txt')
vcf_definition_path = os.path.join(HERE, 'data/example_definition_vcf.txt')
out_path = os.path.join(HERE, 'data/out')

def test_generate_matrix_bedpe_dir():
    output = os.path.join(out_path, 'out1.tsv')
    runner = CliRunner()
    result = runner.invoke(viola, ['generate-feature-matrix', '--input-dir', bedpe_path, '--format', 'bedpe', '--definitions', bedpe_definition_path, output])
    assert result.exit_code == 0

def test_generate_matrix_bedpe_files():
    output = os.path.join(out_path, 'out2.tsv')
    path1 = os.path.join(bedpe_path, 'bedpe1.bedpe')
    path2 = os.path.join(bedpe_path, 'bedpe2.bedpe')
    path = ','.join([path1, path2])
    runner = CliRunner()
    result = runner.invoke(viola, ['generate-feature-matrix', '--input-files', path, '--format', 'bedpe', '--definitions', bedpe_definition_path, output])
    assert result.exit_code == 0

def test_generate_matrix_bedpe_files2():
    output = os.path.join(out_path, 'out3.tsv')
    path1 = os.path.join(bedpe_path, 'bedpe1.bedpe')
    path2 = os.path.join(bedpe_path, 'bedpe2.bedpe')
    path = ','.join([path1, path2])
    runner = CliRunner()
    result = runner.invoke(viola, ['generate-feature-matrix', '--input-files', path, '--input-files-id', 'bedpe1,bedpe2', '--format', 'bedpe', '--definitions', bedpe_definition_path, output])
    assert result.exit_code == 0

def test_generate_matrix_vcf_dir():
    output = os.path.join(out_path, 'out4.tsv')
    runner = CliRunner()
    result = runner.invoke(viola, ['generate-feature-matrix', '--input-dir', vcf_path, '--format', 'vcf', '--definitions', vcf_definition_path, output])
    assert result.exit_code == 0

def test_generate_matrix_vcf_dir_as_bp():
    output = os.path.join(out_path, 'out5.tsv')
    runner = CliRunner()
    result = runner.invoke(viola, ['generate-feature-matrix', '--input-dir', vcf_path, '--format', 'vcf', '--definitions', vcf_definition_path, '--as-breakpoint', output])
    assert result.exit_code == 0

def test_generate_matrix_vcf_files():
    path1 = os.path.join(vcf_path, 'manta1.vcf')
    path2 = os.path.join(vcf_path, 'manta2.vcf')
    path = ','.join([path1, path2])
    output = os.path.join(out_path, 'out6.tsv')
    runner = CliRunner()
    result = runner.invoke(viola, ['generate-feature-matrix', '--input-files', path, '--format', 'vcf', '--definitions', vcf_definition_path, '--as-breakpoint', output])
    assert result.exit_code == 0

def test_generate_matrix_vcf_files2():
    path1 = os.path.join(vcf_path, 'manta1.vcf')
    path2 = os.path.join(vcf_path, 'manta2.vcf')
    path = ','.join([path1, path2])
    output = os.path.join(out_path, 'out7.tsv')
    runner = CliRunner()
    result = runner.invoke(viola, ['generate-feature-matrix', '--input-files', path, '--input-files-id', 'manta1,manta2', '--format', 'vcf', '--definitions', vcf_definition_path, '--as-breakpoint', output])
    assert result.exit_code == 0