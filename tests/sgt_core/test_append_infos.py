import pytest
import sgt
import sys, os
HERE = os.path.abspath(os.path.dirname(__file__))

class TestAppendInfos:
    manta_path = os.path.join(HERE, 'data/manta1.inv.vcf')
    reader = sgt.read_vcf(manta_path)
    positions = reader.get_table('positions')

    def test_append_svlen(self):
        if 'svlen' in self.reader.table_list:
            result = sgt.Vcf.append_infos(self.reader, self.positions, ls_tablenames=['svlen'])
            assert result.columns[-1] == 'svlen_0'
            assert result['svlen_0'].notnull().any()