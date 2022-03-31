import viola
from viola.testing import assert_vcf_equal
from io import StringIO
import os
import pandas as pd
HEADER = """##fileformat=VCFv4.1
##contig=<ID=chr1,length=195471971>
##contig=<ID=chr2,length=182113224>
##contig=<ID=chr11,length=122082543>
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=CIGAR,Number=A,Type=String,Description="CIGAR alignment for each alternate indel allele">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical homology at event breakpoints">
##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description="Length of insertion">
##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description="Sequence of insertion">
##INFO=<ID=LEFT_SVINSSEQ,Number=.,Type=String,Description="Known left side of insertion for an insertion of unknown length">
##INFO=<ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description="Known right side of insertion for an insertion of unknown length">
##INFO=<ID=CONTIG,Number=1,Type=String,Description="Assembled contig sequence">
##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description="Read depth at local translocation breakend">
##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description="Read depth at remote translocation mate breakend">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
##INFO=<ID=SOMATICSCORE,Number=1,Type=Integer,Description="Somatic variant quality score">
##INFO=<ID=JUNCTION_SOMATICSCORE,Number=1,Type=Integer,Description="If the SV junctino is part of an EVENT (ie. a multi-adjacency variant), this field provides the SOMATICSCORE value for the adjacency in question only">
##INFO=<ID=INV3,Number=0,Type=Flag,Description="Inversion breakends open 3' of reported location">
##INFO=<ID=INV5,Number=0,Type=Flag,Description="Inversion breakends open 5' of reported location">
##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999">
##FILTER=<ID=MaxDepth,Description="Normal sample site depth is greater than 3x the median chromosome depth near one or both variant breakends">
##FILTER=<ID=MinSomaticScore,Description="Somatic score is less than 30">
##FILTER=<ID=MaxMQ0Frac,Description="For a small variant (<1000 bases) in the normal sample, the fraction of reads with MAPQ0 around either breakend exceeds 0.4">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mouse1_N	mouse1_T
"""

body = """chr1	82550461	test1	G	<DEL>	.	MinSomaticScore	SVTYPE=DEL;SVLEN=-3764;END=82554225;CIPOS=-51,52;CIEND=-51,52;SOMATIC;SOMATICSCORE=10	PR:SR	21,0:10,0	43,4:15,3
chr1	22814216	test2	T	<INV>	.	MinSomaticScore;MaxMQ0Frac	SVTYPE=INV;SVLEN=69766915;END=92581131;CIPOS=-51,51;CIEND=-89,90;SOMATIC;SOMATICSCORE=11;INV5	PR	24,0	35,5
chr8	69735694	test4_1	A	A[chr11:30018803[	.	MinSomaticScore	IMPRECISE;SVTYPE=BND;CIPOS=-100,100;MATEID=test4_2;BND_DEPTH=49;MATE_BND_DEPTH=35;SOMATIC;SOMATICSCORE=12	PR	30,1	68,9
chr11	46689527	test3	C	<INV>	.	MinSomaticScore	IMPRECISE;SVTYPE=INV;SVLEN=61679;END=46751206;CIPOS=-64,64;CIEND=-65,65;SOMATIC;SOMATICSCORE=13;INV3	PR	17,1	55,19
chr11	30018803	test4_2	G	]chr8:69735694]G	.	MinSomaticScore	IMPRECISE;SVTYPE=BND;CIPOS=-100,100;MATEID=test4_1;BND_DEPTH=35;MATE_BND_DEPTH=49;SOMATIC;SOMATICSCORE=12	PR	30,1	68,9
chr11	30625198	test5	C	<DUP:TANDEM>	.	PASS	IMPRECISE;SVTYPE=DUP;SVLEN=575363;END=31200561;CIPOS=-27,27;CIEND=-69,70;SOMATIC;SOMATICSCORE=39	PR	14,0	38,10
"""

def test_replace_svid():
    vcf = viola.read_vcf(StringIO(HEADER + body))
    vcf.set_value_for_info_by_id("svlen", "test1", 0, -3)
    vcf.set_value_for_info_by_id("svlen", "test4_1", 0, 1100)
    vcf.set_value_for_info_by_id("svtype", "test4_1", 0, 'INV3')
    vcf.set_value_for_info_by_id("imprecise", "test2", 0, True)
    vcf.set_value_for_info_by_id("imprecise", "test4_1", 0, False)
    vcf.set_value_for_info_by_id("imprecise", "test1", 0, False)
    test_df = vcf.get_table('svlen')
    test_df2 = vcf.get_table('svtype')
    test_df3 = vcf.get_table('imprecise')
    expected_df = pd.read_csv(StringIO("""id;value_idx;svlen
test1;0;-3
test2;0;69766915
test4_1;0;1100
test3;0;61679
test5;0;575363"""), sep=';')
    expected_df2 = pd.read_csv(StringIO("""id;value_idx;svtype
test1;0;DEL
test2;0;INV
test4_1;0;INV3
test4_2;0;BND
test3;0;INV
test5;0;DUP"""), sep=';')
    expected_df3 = pd.read_csv(StringIO("""id;value_idx;imprecise
test2;0;True
test4_2;0;True
test3;0;True
test5;0;True"""), sep=';')
    pd.testing.assert_frame_equal(test_df.sort_values(by='id').reset_index(drop=True), expected_df.sort_values(by='id').reset_index(drop=True))
    pd.testing.assert_frame_equal(test_df2.sort_values(by='id').reset_index(drop=True), expected_df2.sort_values(by='id').reset_index(drop=True))
    pd.testing.assert_frame_equal(test_df3.sort_values(by='id').reset_index(drop=True), expected_df3.sort_values(by='id').reset_index(drop=True))
