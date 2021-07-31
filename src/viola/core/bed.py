from intervaltree import Interval, IntervalTree
import pandas as pd
class Bed(object):
    """
    A class for optimizing the processing of BED files.
    Interval tree algorithm is employed for storing region information.
    This allows for fast retrieval of matching BED rows for queries in the form of genomic coordinates.

    See Also
    ---------
    read_bed: Read a BED file into Bed class.
    """
    def __init__(self, df, header):
        self._df = df
        self._header = header
        self._tree = self._intervaltree_init(df)
    
    def _intervaltree_init(self, df):
        dict_tree = dict()
        for idx, row in df.iterrows():
            chrom = row['chrom']
            st = row['chromStart']
            en = row['chromEnd'] + 1
            if dict_tree.get(chrom) is None:
                dict_tree[chrom] = IntervalTree()
            dict_tree[chrom][st:en] = row
        
        return dict_tree
    
    def query(self, chrom, st, en=None) -> pd.DataFrame:
        if self._tree.get(chrom) is None:
            return pd.DataFrame(columns = self._df.columns)
        tree = self._tree[chrom]
        if en is None:
            set_out = tree[st]
        else:
            en += 1
            set_out = tree[st:en]
        if len(set_out) == 0:
            return pd.DataFrame(columns = self._df.columns)
        ls_out = list(set_out)
        ls_ser_out = [x.data for x in ls_out]
        return pd.DataFrame(ls_ser_out)
