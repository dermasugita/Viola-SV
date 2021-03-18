import viola
import pytest
from viola._exceptions import (
    TableNotFoundError,
)

def test_get_table():
    data = viola.read_vcf('resources/vcf/manta1.vcf')

    # get positions table
    positions = data.get_table('positions')

    ## test for TableNotFoundError
    with pytest.raises(TableNotFoundError):
        data.get_table("notexists")