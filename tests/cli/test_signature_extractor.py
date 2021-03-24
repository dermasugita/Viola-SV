from viola.cli.viola import viola
from click.testing import CliRunner
import pandas as pd
import sys, os
HERE = os.path.abspath(os.path.dirname(__file__))
infile = os.path.join(HERE, 'data/feature_matrix.csv')
outfile = os.path.join(HERE, 'data/out/test_signature.tsv')

def test_signature_extractor():
    runner = CliRunner()
    result = runner.invoke(viola, ['extract-signature', infile, outfile])
    assert result.exit_code == 0
