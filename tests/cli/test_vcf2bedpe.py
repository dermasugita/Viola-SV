from viola.cli.vcf2bedpe import vcf2bedpe
from click.testing import CliRunner
import sys, os
HERE = os.path.abspath(os.path.dirname(__file__))

def test_vcf2bedpe_manta():
  path = os.path.join(HERE, '../io/data/test.manta.vcf')
  runner = CliRunner()
  result = runner.invoke(vcf2bedpe, [path])
  assert result.exit_code == 0


def test_vcf2bedpe_delly():
  path = os.path.join(HERE, '../io/data/test.delly.vcf')
  runner = CliRunner()
  result = runner.invoke(vcf2bedpe, ['--caller','delly', path])
  assert result.exit_code == 0

def test_vcf2bedpe_lumpy():
  path = os.path.join(HERE, '../io/data/test.lumpy.vcf')
  runner = CliRunner()
  result = runner.invoke(vcf2bedpe, ['--caller=lumpy', path])
  assert result.exit_code == 0

def test_vcf2bedpe_gridss():
  path = os.path.join(HERE, '../io/data/test.gridss.vcf')
  runner = CliRunner()
  result = runner.invoke(vcf2bedpe, ['--caller=gridss', path])
  assert result.exit_code == 0
