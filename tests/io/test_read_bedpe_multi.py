import viola
import pandas as pd
import os
HERE = os.path.abspath(os.path.dirname(__file__))

def test_read_bedpe_multi_exclude_empty():
    bedpe = viola.read_bedpe_multi(os.path.join(HERE, 'data/multibedpe'), exclude_empty_cases=True)
    assert set(bedpe._ls_patients) == {'test1', 'test2'}

def test_read_bedpe_multi_keep_empty():
    bedpe = viola.read_bedpe_multi(os.path.join(HERE, 'data/multibedpe'), exclude_empty_cases=False)
    assert set(bedpe._ls_patients) == {'test1', 'test2', 'empty'}