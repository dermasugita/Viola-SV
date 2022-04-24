import pytest
import viola
import sys, os
import numpy as np
import pandas as pd
from io import StringIO
delly_HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20201013
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##FILTER=<ID=LowQual,Description="Poor quality and insufficient number of PEs and SRs.">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for POS2 coordinate in case of an inter-chromosomal translocation">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="Genomic position for CHR2 in case of an inter-chromosomal translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description="Median mapping quality of split-reads">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CE,Number=1,Type=Float,Description="Consensus sequence entropy">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Predicted microhomology length using a max. edit distance of 2">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10-scaled genotype likelihoods for RR,RA,AA genotypes">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Per-sample genotype filter">
##FORMAT=<ID=RC,Number=1,Type=Integer,Description="Raw high-quality read counts or base counts for the SV">
##FORMAT=<ID=RCL,Number=1,Type=Integer,Description="Raw high-quality read counts or base counts for the left control region">
##FORMAT=<ID=RCR,Number=1,Type=Integer,Description="Raw high-quality read counts or base counts for the right control region">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Read-depth based copy-number estimate for autosomal sites">
##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">
##FORMAT=<ID=RR,Number=1,Type=Integer,Description="# high-quality reference junction reads">
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="# high-quality variant junction reads">
##reference=/data/share/iGenomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/genome.fa
##contig=<ID=chr10,length=130694993>
##contig=<ID=chr11,length=122082543>
##INFO=<ID=RDRATIO,Number=1,Type=Float,Description="Read-depth ratio of tumor vs. normal.">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic structural variant.">
##bcftools_viewVersion=1.8+htslib-1.8
##bcftools_viewCommand=view mouse1_filtered.bcf; Date=Wed Oct 14 15:37:56 2020
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mouse01_T	mouse01_N
"""

delly_body = """chr10	3485568	DEL00001	C	<DEL>	180	PASS	PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.8.5;END=3485617;PE=0;MAPQ=0;CT=3to5;CIPOS=-38,38;CIEND=-38,38;SRMAPQ=60;INSLEN=0;HOMLEN=38;SR=3;SRQ=1;CONSENSUS=TAGGCATGGTTCTGGAGCAGTAACTAAGAACTTATATTTGATCCATAAGCAGGAGGCAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGAGAGAGAGAGAGAGAGAGAGAGCCTAACTCAGAATAGCTTGGGTTTTTGAAAC;CE=1.84411;RDRATIO=0.832798;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-12.9844,0,-41.0803:130:PASS:1132:1866:1406:1:0:0:15:6	0/0:0,-2.70665,-29.1974:27:PASS:416:859:557:2:0:0:9:0
chr10	8825909	INV00001	A	<INV>	17	LowQual	IMPRECISE;SVTYPE=INV;SVMETHOD=EMBL.DELLYv0.8.5;END=32816608;PE=2;MAPQ=12;CT=5to5;CIPOS=-209,209;CIEND=-209,209;RDRATIO=0.987391;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-2.47951,0,-66.449:25:PASS:1411949:5924972:2958545:3:13:7:0:0	0/0:0,-2.10716,-40:21:PASS:617135:2644311:1308820:3:7:0:0:0
chr10	11350703	DUP00001	G	<DUP>	29	LowQual	IMPRECISE;SVTYPE=DUP;SVMETHOD=EMBL.DELLYv0.8.5;END=122603091;PE=2;MAPQ=18;CT=5to3;CIPOS=-119,119;CIEND=-119,119;RDRATIO=0.999221;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/0:0,-0.108547,-39.5411:4:LowQual:2041981:27256025:1977420:14:8:2:0:0	0/0:0,-2.68066,-40.1714:27:PASS:899570:12088575:881723:14:9:0:0:0
chr10	61016659	INV00002	G	<INV>	10	LowQual	IMPRECISE;SVTYPE=INV;SVMETHOD=EMBL.DELLYv0.8.5;END=111424580;PE=2;MAPQ=5;CT=3to3;CIPOS=-203,203;CIEND=-203,203;RDRATIO=1.00121;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/0:0,-1.01236,-51.7412:10:LowQual:6123243:12379748:4722845:2:9:3:0:0	0/0:0,-0.903089,-18:10:LowQual:2709579:5484294:2101093:2:3:0:0:0
chr11	15249914	BND00001	A	A]chr10:72297333]	193	PASS	IMPRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv0.8.5;END=15249915;CHR2=chr10;POS2=72297333;PE=5;MAPQ=37;CT=3to3;CIPOS=-378,378;CIEND=-378,378;RDRATIO=1;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-14.872,0,-35.0861:149:PASS:32721:68453:35732:2:10:5:0:0	0/0:0,-1.20195,-17.3978:12:LowQual:14993:27908:12915:2:4:0:0:0
"""
delly_body_buf = StringIO(delly_HEADER + delly_body)

class TestReadVcfLumpy:
    vcf_delly = viola.read_vcf(delly_body_buf, variant_caller='delly', patient_name='delly')
    def test_read_vcf_delly(self):
        vcf = self.vcf_delly

        # assert positions table equal
        df_pos_expected = pd.read_csv(StringIO(
            """id,chrom1,pos1,chrom2,pos2,strand1,strand2,ref,alt,qual,svtype
DEL00001,chr10,3485568,chr10,3485618,+,-,C,<DEL>,180,DEL
INV00001,chr10,8825910,chr10,32816609,-,-,A,<INV>,17,INV
DUP00001,chr10,11350704,chr10,122603091,-,+,G,<DUP>,29,DUP
INV00002,chr10,61016659,chr10,111424580,+,+,G,<INV>,10,INV
BND00001,chr11,15249914,chr10,72297333,+,+,A,A]chr10:72297333],193,BND
"""
        ))
        pd.testing.assert_frame_equal(vcf.positions, df_pos_expected)

    