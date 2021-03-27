import viola
import pytest
from viola._exceptions import (
    TableNotFoundError,
)
<<<<<<< HEAD
import os
HERE = os.path.abspath(os.path.dirname(__file__))

def test_get_table():
    data = viola.read_vcf(os.path.join(HERE, 'data/test.manta.vcf'))
=======

def test_get_table():
    data = viola.read_vcf('resources/vcf/manta1.vcf')
>>>>>>> 5b88c4e09abac3da91babb308491bb51425a6705

    # get positions table
    positions = data.get_table('positions')

    ## test for TableNotFoundError
    with pytest.raises(TableNotFoundError):
        data.get_table("notexists")