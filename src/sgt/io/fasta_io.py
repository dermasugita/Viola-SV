from collections import OrderedDict
from sgt.core.fasta import Fasta
def read_fasta(path):
    from Bio import SeqIO
    fasta_import = SeqIO.parse(path, "fasta")
    fasta = Fasta()
    for seq in fasta_import:
        fasta[seq.id] = seq.seq
    return fasta
