import os
import viola
HERE = os.path.abspath(os.path.dirname(__file__))
filepath = os.path.join(HERE, '../../resources/vcf/manta1.inv.vcf')
bedpepath1 = os.path.join(HERE, '../../resources/pcawg/10ad692b-4c3d-42de-9b5e-4968441388b3.pcawg_consensus_1.6.161116.somatic.sv.bedpe')
bedpepath2 = os.path.join(HERE, '../../resources/pcawg/0b19bee7-5281-4915-9d98-c20eb3e84ecf.pcawg_consensus_1.6.161116.somatic.sv.bedpe')
bedpath = os.path.join(HERE, '../../resources/bed/fragile_site.hg19.bed')
def test_hello():
    bedpe1 = viola.read_bedpe(bedpepath1)
    bedpe2 = viola.read_bedpe(bedpepath2)
    viola.MultiBedpe([bedpe1, bedpe2], ['p1', 'p2'])
    bed = viola.read_bed(bedpath)
    bedpe1.annotate_bed(bed, 'fragile')
    bedpe2.annotate_bed(bed, 'fragile')
    test = viola.get_id_by_boolean_info(bedpe2, info='fragileleft', svtype='DEL')
    ls_func = [
        fragile,
        del1,
        del2,
        del3,
        del4,
        dup1,
        dup2,
        dup3,
        dup4
    ]
    ls_names = [
        'fragile',
        'DEL::<50k',
        'DEL::<500k',
        'DEL::<5M',
        'DEL::>=5M',
        'DUP::<50k',
        'DUP::<500k',
        'DUP::<5M',
        'DUP::>=5M'
    ]
    print(bedpe1.classify_manual_svtype(ls_func, ls_names))

def fragile(x):
    set1 = viola.get_id_by_boolean_info(x, info='fragileleft', svtype='DEL')
    set2 = viola.get_id_by_boolean_info(x, info='fragileright', svtype='DEL')
    return set1 | set2

def del1(x):
    return viola.get_id_by_slicing_info(x, info='svlen', st=-50000, svtype='DEL')
def del2(x):
    return viola.get_id_by_slicing_info(x, info='svlen', st=-500000, en=-50000, svtype='DEL')
def del3(x):
    return viola.get_id_by_slicing_info(x, info='svlen', st=-5000000, en=-500000, svtype='DEL')
def del4(x):
    return viola.get_id_by_slicing_info(x, info='svlen', en=-5000000, svtype='DEL')

def dup1(x):
    return viola.get_id_by_slicing_info(x, info='svlen', en=50000, svtype='DUP')
def dup2(x):
    return viola.get_id_by_slicing_info(x, info='svlen', st=50000, en=500000, svtype='DUP')
def dup3(x):
    return viola.get_id_by_slicing_info(x, info='svlen', st=500000, en=5000000, svtype='DUP')
def dup4(x):
    return viola.get_id_by_slicing_info(x, info='svlen', st=5000000, svtype='DUP')

