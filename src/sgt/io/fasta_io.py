from collections import OrderedDict
import requests
from io import StringIO
from Bio import SeqIO
from sgt.core.fasta import Fasta
from sgt.utils.utils import is_url

def read_fasta(path_or_url):
    if is_url(path_or_url):
        b = StringIO(requests.get(path_or_url).text)
    else:
        b = path_or_url
    fasta_import = SeqIO.parse(b, "fasta")
    fasta = Fasta()
    for seq in fasta_import:
        fasta[seq.id] = seq.seq
    return fasta
