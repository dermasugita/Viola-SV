from collections import OrderedDict
import urllib
# import requests ミスった!! 保留。
import gzip
from io import StringIO, BytesIO
from Bio import SeqIO
from viola.core.fasta import Fasta
from viola.utils.utils import is_url

def read_fasta(path_or_url):
    if is_url(path_or_url):
        if path_or_url.split('.')[-1] == 'gz':
            response = urllib.request.urlopen(path_or_url)
            bytes_response = response.read()
            text = gzip.open(BytesIO(bytes_response), 'rt')
        else:
            pass
            #text = StringIO(requests.get(path_or_url).text)
    else:
        text = path_or_url
    fasta_import = SeqIO.parse(text, "fasta")
    fasta = Fasta()
    for seq in fasta_import:
        fasta[seq.id] = seq.seq
    return fasta
