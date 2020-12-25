import sgt
import Bio
import os
HERE = os.path.abspath(os.path.dirname(__file__))
fasta_path = os.path.join(HERE, '../../../SV_parser_test/jupyters/fasta/hg19.fa')
bedpepath1 = os.path.join(HERE, '../../resources/pcawg/10ad692b-4c3d-42de-9b5e-4968441388b3.pcawg_consensus_1.6.161116.somatic.sv.bedpe')
bedpepath2 = os.path.join(HERE, '../../resources/pcawg/0b19bee7-5281-4915-9d98-c20eb3e84ecf.pcawg_consensus_1.6.161116.somatic.sv.bedpe')

def test_get_microhomology():
    fasta = sgt.read_fasta(fasta_path)
    bedpe = sgt.read_bedpe(bedpepath2)
    bedpe.get_microhomology(fasta)
    print(bedpe.to_bedpe_like(custom_infonames=['homlen', 'homseq']))