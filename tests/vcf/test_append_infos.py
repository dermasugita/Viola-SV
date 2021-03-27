import pytest
import viola
import sys, os
HERE = os.path.abspath(os.path.dirname(__file__))

class TestAppendInfos:
<<<<<<< HEAD
    manta_path = os.path.join(HERE, 'data/test.manta.vcf')
=======
    manta_path = os.path.join(HERE, 'data/manta1.inv.vcf')
>>>>>>> 5b88c4e09abac3da91babb308491bb51425a6705
    reader = viola.read_vcf(manta_path)
    positions = reader.get_table('positions')

    def test_append_svlen(self):
        if 'svlen' in self.reader.table_list:
            result = viola.Vcf.append_infos(self.reader, self.positions, ls_tablenames=['svlen'])
            assert result.columns[-1] == 'svlen_0'
            assert result['svlen_0'].notnull().any()