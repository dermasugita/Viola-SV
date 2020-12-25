from collections import OrderedDict
def read_fasta(path):
    from Bio import SeqIO
    fasta_import = SeqIO.parse(path, "fasta")
    fasta = OrderedDict()
    for seq in fasta_import:
        fasta[seq.id] = seq.seq
    return fasta
