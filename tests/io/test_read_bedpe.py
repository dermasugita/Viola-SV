import sgt
import pandas as pd
from io import StringIO

class TestReadBedpe:
    data = """chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2
chr1\t10\t13\tchr2\t20\t30\ttest1\t10\t+\t+
chr3\t100\t130\tchr3\t210\t230\ttest2\t30\t-\t+
"""
    b = StringIO(data)
    obj = sgt.read_bedpe(b)
    def test_read_bedpe(self):
        b = StringIO(self.data)
        sgt.read_bedpe(b)
    
    def test_svpos(self):
        expected_data = """id\tchrom1\tpos1\tchrom2\tpos2\tstrand1\tstrand2\tqual\tsvtype\tref\talt
test1\tchr1\t11\tchr2\t25\t+\t+\t10\tBND\tN\tN]chr2:25]
test2\tchr3\t115\tchr3\t220\t-\t+\t30\tDUP\tN\t<DUP>
        """
        df_svpos = self.obj.get_table('positions') 
        df_expected = pd.read_csv(StringIO(expected_data), sep="\t")
        pd.testing.assert_frame_equal(df_svpos, df_expected)
    
    def test_create_alt_field_from_position(self):
        test_data = """id\tchrom1\tpos1\tchrom2\tpos2\tstrand1\tstrand2\tsvtype\tref
test1\tchr1\t10\tchr2\t10\t+\t-\tBND\tN
test2\tchr1\t10\tchr1\t10\t+\t-\tDEL\tN
"""
        b = StringIO(test_data)
        df_svpos = pd.read_csv(b, sep="\t")
        result = sgt.io.parser.create_alt_field_from_position(df_svpos)
        