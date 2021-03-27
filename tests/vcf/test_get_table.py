import viola
import pytest
from viola._exceptions import (
    TableNotFoundError,
)
import os
HERE = os.path.abspath(os.path.dirname(__file__))

def test_get_table():
    data = viola.read_vcf(os.path.join(HERE, 'data/test.manta.vcf'))

    # get positions table
    positions = data.get_table('positions')

    ## test for TableNotFoundError
    with pytest.raises(TableNotFoundError):
        data.get_table("notexists")