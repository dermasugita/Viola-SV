import viola
from viola.testing import assert_vcf_equal
from io import StringIO
import os
HEADER ="""##fileformat=VCFv4.2
##source=LUMPY
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">
##INFO=<ID=SECONDARY,Number=0,Type=Flag,Description="Secondary breakend in a multi-line variants">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="Confidence interval (95%) around POS for imprecise variants">
##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="Confidence interval (95%) around END for imprecise variants">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">
##INFO=<ID=BD,Number=.,Type=Integer,Description="Amount of BED evidence supporting the variant across all samples">
##INFO=<ID=EV,Number=.,Type=String,Description="Type of LUMPY evidence contributing to the variant call">
##INFO=<ID=PRPOS,Number=.,Type=String,Description="LUMPY probability curve of the POS breakend">
##INFO=<ID=PREND,Number=.,Type=String,Description="LUMPY probability curve of the END breakend">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting the variant">
##FORMAT=<ID=BD,Number=1,Type=Integer,Description="Amount of BED evidence supporting the variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1T	1N
"""
HEADER_expected ="""##fileformat=VCFv4.2
##source=LUMPY
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="Confidence interval (95%) around POS for imprecise variants">
##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="Confidence interval (95%) around END for imprecise variants">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">
##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">
##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">
##INFO=<ID=BD,Number=.,Type=Integer,Description="Amount of BED evidence supporting the variant across all samples">
##INFO=<ID=EV,Number=.,Type=String,Description="Type of LUMPY evidence contributing to the variant call">
##INFO=<ID=PRPOS,Number=.,Type=String,Description="LUMPY probability curve of the POS breakend">
##INFO=<ID=PREND,Number=.,Type=String,Description="LUMPY probability curve of the END breakend">
##INFO=<ID=ORGBEID,Number=2,Type=String,Description="Breakend ID which were in original VCF file.",Source="Python package, Viola-SV.">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP:TANDEM,Description="Tandem duplication">
##ALT=<ID=INS,Description="Insertion of novel sequence">
##ALT=<ID=CNV,Description="Copy number variable region">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SU,Number=1,Type=Integer,Description="Number of pieces of evidence supporting the variant">
##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-end reads supporting the variant">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Number of split reads supporting the variant">
##FORMAT=<ID=BD,Number=1,Type=Integer,Description="Amount of BED evidence supporting the variant">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1T	1N
"""

body = """chr10	39984191	43	N	<INV>	.	.	SVTYPE=INV;STRANDS=++:4,--:2;SVLEN=176;END=39984367;CIPOS=-8,11;CIEND=-17,31;CIPOS95=-3,6;CIEND95=-1,15;IMPRECISE;SU=6;PE=6;SR=0	GT:SU:PE:SR	./.:3:3:0	./.:3:3:0
chr10	121422696	788	N	<DEL>	.	.	SVTYPE=DEL;STRANDS=+-:8;SVLEN=-21;END=121422717;CIPOS=-9,8;CIEND=-9,8;CIPOS95=0,0;CIEND95=0,0;SU=8;PE=0;SR=8	GT:SU:PE:SR	./.:8:0:8	./.:0:0:0
chr14	45812810	7294_1	N	[chr14:65332314[N	.	.	SVTYPE=BND;STRANDS=--:10;EVENT=7294;MATEID=7294_2;CIPOS=-7,6;CIEND=-7,6;CIPOS95=0,0;CIEND95=0,0;SU=10;PE=0;SR=10	GT:SU:PE:SR	./.:10:0:10	./.:0:0:0
chr14	65332314	7294_2	N	[chr14:45812810[N	.	.	SVTYPE=BND;STRANDS=--:10;SECONDARY;EVENT=7294;MATEID=7294_1;CIPOS=-7,6;CIEND=-7,6;CIPOS95=0,0;CIEND95=0,0;SU=10;PE=0;SR=10	GT:SU:PE:SR	./.:10:0:10	./.:0:0:0
chr15	60882465	8801	N	<DUP>	.	.	SVTYPE=DUP;STRANDS=-+:41;SVLEN=3099777;END=63982242;CIPOS=-3,1;CIEND=-3,2;CIPOS95=0,0;CIEND95=0,0;SU=41;PE=19;SR=22	GT:SU:PE:SR	./.:41:19:22	./.:0:0:0
"""

body_expected = """chr10	39984191	43	N	<INV>	.	.	SVTYPE=INV;STRANDS=++:4,--:2;SVLEN=176;END=39984367;CIPOS=-8,11;CIEND=-17,31;CIPOS95=-3,6;CIEND95=-1,15;IMPRECISE;SU=6;PE=6;SR=0	GT:SU:PE:SR	./.:3:3:0	./.:3:3:0
chr10	121422696	788	N	<DEL>	.	.	SVTYPE=DEL;STRANDS=+-:8;SVLEN=-21;END=121422717;CIPOS=-9,8;CIEND=-9,8;CIPOS95=0,0;CIEND95=0,0;SU=8;PE=0;SR=8	GT:SU:PE:SR	./.:8:0:8	./.:0:0:0
chr14	45812810	viola_breakpoint:0	N	[chr14:65332314[N	.	.	SVTYPE=BND;STRANDS=--:10;EVENT=7294;MATEID=7294_2;CIPOS=-7,6;CIEND=-7,6;CIPOS95=0,0;CIEND95=0,0;SU=10;PE=0;SR=10;ORGBEID=7294_1,7294_2	GT:SU:PE:SR	./.:10:0:10	./.:0:0:0
chr15	60882465	8801	N	<DUP>	.	.	SVTYPE=DUP;STRANDS=-+:41;SVLEN=3099777;END=63982242;CIPOS=-3,1;CIEND=-3,2;CIPOS95=0,0;CIEND95=0,0;SU=41;PE=19;SR=22	GT:SU:PE:SR	./.:41:19:22	./.:0:0:0
"""


def test_breakend2breakpoint():
    vcf = viola.read_vcf(StringIO(HEADER + body), variant_caller='lumpy')
    vcf_expected = viola.read_vcf(StringIO(HEADER_expected + body_expected), variant_caller='lumpy')
    ex_svpos = vcf_expected.get_table('positions')
    ex_svpos.loc[3, 'svtype'] = 'INV'
    ex_svtype = vcf_expected.get_table('svtype')
    ex_svtype.iloc[3, 2] = 'INV'
    ex_infos_meta = vcf_expected.get_table('infos_meta')
    ex_infos_meta.index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 17]
    vcf_expected._odict_alltables['positions'] = ex_svpos
    vcf_expected._odict_alltables['svtype'] = ex_svtype
    vcf_expected._odict_alltables['infos_meta'] = ex_infos_meta
    vcf_result = vcf.breakend2breakpoint()
    assert_vcf_equal(vcf_result, vcf_expected)