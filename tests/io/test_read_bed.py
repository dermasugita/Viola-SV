import viola
import os
HERE = os.path.abspath(os.path.dirname(__file__))
f_path = os.path.join(HERE, 'data/test.bed')
f_path_gzip = os.path.join(HERE, 'data/compressed/test.bed.gz')

def test_read_bed():
    bed = viola.read_bed(f_path)
    print(bed)

def test_read_bed_gzip():
    bed = viola.read_bed(f_path_gzip)
    print(bed)