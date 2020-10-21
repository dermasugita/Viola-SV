import sgt
import pytest
from sgt._exceptions import (
    TableNotFoundError,
)

def test_get_table():
    data = sgt.read_vcf('resources/vcf/manta1.vcf')

    ## test for TableNotFoundError
    with pytest.raises(TableNotFoundError):
        data.get_table("notexists")