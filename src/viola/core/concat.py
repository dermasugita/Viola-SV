import pandas as pd
import viola
from viola.core.bedpe import Bedpe
from viola.core.vcf import Vcf
from viola.core.cohort import MultiBedpe, MultiVcf

def concat_bedpe(bedpe1, bedpe2):
    # checking type
    if not isinstance(bedpe1, Bedpe) | isinstance(bedpe1, MultiBedpe):
        raise TypeError('Input values should be Bedpe or MultiBedpe class.')
    if not isinstance(bedpe2, Bedpe) | isinstance(bedpe2, MultiBedpe):
        raise TypeError('Input values should be Bedpe or MultiBedpe class.')
    is_multi = isinstance(bedpe1, MultiBedpe)

def concat(ls_objs):
    """
    concat(ls_objs)
    """
    return concat_bedpe(ls_objs[0], ls_objs[1])
    