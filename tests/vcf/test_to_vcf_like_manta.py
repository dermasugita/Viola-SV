import viola
import sys, os
import pandas as pd
from io import StringIO
HERE = os.path.abspath(os.path.dirname(__file__))
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
class TestToVcf:
    #manta_path = os.path.join(HERE, 'data/manta1.inv.vcf')
    #result = viola.read_vcf(manta_path)
    body = """chr1	82550461	test1	G	<DEL>	.	MinSomaticScore	END=82554225;SVTYPE=DEL;SVLEN=-3764;IMPRECISE;CIPOS=-51,52;CIEND=-51,52;SOMATIC;SOMATICSCORE=10	PR:SR	21,0:10,0	43,4:15,3
chr1	22814216	test2	T	<INV>	.	MinSomaticScore;MaxMQ0Frac	END=92581131;SVTYPE=INV;SVLEN=69766915;CIPOS=-51,51;CIEND=-89,90;SOMATIC;SOMATICSCORE=11;INV5	PR	24,0	35,5
chr8	69735694	test4_1	A	A[chr11:330018803[	.	MinSomaticScore	IMPRECISE;SVTYPE=BND;CIPOS=-50,50;MATEID=test4_2;BND_DEPTH=35;MATE_BND_DEPTH=49;SOMATIC;SOMATICSCORE=12	PR	30,1	68,9
chr11	46689527	test3	C	<INV>	.	MinSomaticScore	END=46751206;SVTYPE=INV;SVLEN=61679;IMPRECISE;CIPOS=-64,64;CIEND=-65,65;SOMATIC;SOMATICSCORE=13;INV3	PR	17,1	55,19
chr11	30018803	test4_2	G	]chr8:69735694]G	.	MinSomaticScore	IMPRECISE;SVTYPE=BND;CIPOS=-100,100;MATEID=test4_1;BND_DEPTH=49;MATE_BND_DEPTH=35;SOMATIC;SOMATICSCORE=12	PR	30,1	68,9
chr11	30625198	test5	C	<DUP:TANDEM>	.	PASS	END=31200561;SVTYPE=DUP;SVLEN=575363;IMPRECISE;CIPOS=-27,27;CIEND=-69,70;SOMATIC;SOMATICSCORE=39	PR	14,0	38,10
"""
    expected_out = """chr1	82550461	test1	G	<DEL>	.	MinSomaticScore	IMPRECISE;SVTYPE=DEL;SVLEN=-3764;END=82554225;CIPOS=-51,52;CIEND=-51,52;SOMATIC;SOMATICSCORE=10	PR:SR	21,0:10,0	43,4:15,3
chr1	22814216	test2	T	<INV>	.	MinSomaticScore;MaxMQ0Frac	SVTYPE=INV;SVLEN=69766915;END=92581131;CIPOS=-51,51;CIEND=-89,90;SOMATIC;SOMATICSCORE=11;INV5	PR	24,0	35,5
chr8	69735694	test4_1	A	A[chr11:330018803[	.	MinSomaticScore	IMPRECISE;SVTYPE=BND;CIPOS=-50,50;CIEND=-100,100;MATEID=test4_2;BND_DEPTH=35;MATE_BND_DEPTH=49;SOMATIC;SOMATICSCORE=12	PR	30,1	68,9
chr11	46689527	test3	C	<INV>	.	MinSomaticScore	IMPRECISE;SVTYPE=INV;SVLEN=61679;END=46751206;CIPOS=-64,64;CIEND=-65,65;SOMATIC;SOMATICSCORE=13;INV3	PR	17,1	55,19
chr11	30018803	test4_2	G	]chr8:69735694]G	.	MinSomaticScore	IMPRECISE;SVTYPE=BND;CIPOS=-100,100;CIEND=-50,50;MATEID=test4_1;BND_DEPTH=49;MATE_BND_DEPTH=35;SOMATIC;SOMATICSCORE=12	PR	30,1	68,9
chr11	30625198	test5	C	<DUP:TANDEM>	.	PASS	IMPRECISE;SVTYPE=DUP;SVLEN=575363;END=31200561;CIPOS=-27,27;CIEND=-69,70;SOMATIC;SOMATICSCORE=39	PR	14,0	38,10
"""
    b = StringIO(HEADER + body)
    result = viola.read_vcf(b, patient_name='patient1')
    def test_to_vcf_like(self):
        df_vcf = self.result.to_vcf_like()
        df_expected = pd.read_csv(StringIO(self.expected_out), index_col=False, names=df_vcf.columns, sep='\t')
        print(df_vcf)
        pd.testing.assert_frame_equal(df_vcf, df_expected)
        