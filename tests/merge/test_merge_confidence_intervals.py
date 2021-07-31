import viola
import sys, os
import pandas as pd
from io import StringIO
HERE = os.path.abspath(os.path.dirname(__file__))

DATA1 = """##fileformat=VCFv4.1
##variantcaller=manta
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mouse01_N	mouse01_T
chr1	100000	M1	N	<DEL>	.	MinSomaticScore	IMPRECISE;SVTYPE=DEL;SVLEN=-100000;END=200000;CIPOS=-51,52;CIEND=-51,52;SOMATIC;SOMATICSCORE=10	PR:SR	21,0:10,0	43,4:15,3
chr1	200030	MD1	N	<DEL>	.	MinSomaticScore	IMPRECISE;SVTYPE=DEL;SVLEN=-100000;END=299920;CIPOS=-18,18;CIEND=-51,52;SOMATIC;SOMATICSCORE=10	PR:SR	21,0:10,0	43,4:15,3
chr1	200000	MDL1	N	<DUP:TANDEM>	.	MinSomaticScore	IMPRECISE;SVTYPE=DUP;SVLEN=100000;END=300000;CIPOS=-30,30;CIEND=-31,31;SOMATIC;SOMATICSCORE=10	PR:SR	21,0:10,0	43,4:15,3
"""
DATA2 = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=20201013
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=BND,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##FILTER=<ID=LowQual,Description="Poor quality and insufficient number of PEs and SRs.">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for POS2 coordinate in case of an inter-chromosomal translocation">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="Genomic position for CHR2 in case of an inter-chromosomal translocation">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">
##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">
##INFO=<ID=CT,Number=1,Type=String,Description="Paired-end signature induced connection type">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="PE confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="PE confidence interval around END">
##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description="Median mapping quality of split-reads">
##INFO=<ID=INSLEN,Number=1,Type=Integer,Description="Predicted length of the insertion">
##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description="Predicted microhomology length using a max. edit distance of 2">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=SRQ,Number=1,Type=Float,Description="Split-read consensus alignment quality">
##INFO=<ID=CONSENSUS,Number=1,Type=String,Description="Split-read consensus sequence">
##INFO=<ID=CE,Number=1,Type=Float,Description="Consensus sequence entropy">
##INFO=<ID=RDRATIO,Number=1,Type=Float,Description="Read-depth ratio of tumor vs. normal.">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic structural variant.">
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
##bcftools_viewVersion=1.8+htslib-1.8
##bcftools_viewCommand=view mouse1_filtered.bcf; Date=Wed Oct 14 15:37:56 2020
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mouse01_T	mouse01_N
chr1	100030	D1	N	<DEL>	180	PASS	PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.8.5;END=600000;PE=0;MAPQ=0;CT=3to5;CIPOS=-38,38;CIEND=-38,38;SRMAPQ=60;INSLEN=0;HOMLEN=38;SR=3;SRQ=1;CONSENSUS=TAGGCATGGTTCTGGAGCAGTAACTAAGAACTTATATTTGATCCATAAGCAGGAGGCAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGAGAGAGAGAGAGAGAGAGAGAGCCTAACTCAGAATAGCTTGGGTTTTTGAAAC;CE=1.84411;RDRATIO=0.832798;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-12.9844,0,-41.0803:130:PASS:1132:1866:1406:1:0:0:15:6	0/0:0,-2.70665,-29.1974:27:PASS:416:859:557:2:0:0:9:0
chr1	200030	MD1	N	<DEL>	180	PASS	PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.8.5;END=299920;PE=0;MAPQ=0;CT=3to5;CIPOS=-18,18;CIEND=-88,88;SRMAPQ=60;INSLEN=0;HOMLEN=38;SR=3;SRQ=1;CONSENSUS=TAGGCATGGTTCTGGAGCAGTAACTAAGAACTTATATTTGATCCATAAGCAGGAGGCAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGAGAGAGAGAGAGAGAGAGAGAGCCTAACTCAGAATAGCTTGGGTTTTTGAAAC;CE=1.84411;RDRATIO=0.832798;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/1:-12.9844,0,-41.0803:130:PASS:1132:1866:1406:1:0:0:15:6	0/0:0,-2.70665,-29.1974:27:PASS:416:859:557:2:0:0:9:0
chr1	200020	MDL1	N	<DUP>	29	LowQual	IMPRECISE;SVTYPE=DUP;SVMETHOD=EMBL.DELLYv0.8.5;END=300080;PE=2;MAPQ=18;CT=5to3;CIPOS=-31,31;CIEND=-61,61;RDRATIO=0.999221;SOMATIC	GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV	0/0:0,-0.108547,-39.5411:4:LowQual:2041981:27256025:1977420:14:8:2:0:0	0/0:0,-2.68066,-40.1714:27:PASS:899570:12088575:881723:14:9:0:0:0
"""

def test_merge():
    manta = viola.read_vcf(os.path.join(HERE, 'data/test.merge.manta.vcf'), variant_caller='manta')
    delly = viola.read_vcf(os.path.join(HERE, 'data/test.merge.delly.vcf'), variant_caller='delly')
    lumpy = viola.read_vcf(os.path.join(HERE, 'data/test.merge.lumpy.vcf'), variant_caller='lumpy')
    gridss = viola.read_vcf(os.path.join(HERE, 'data/test.merge.gridss.vcf'), variant_caller='gridss')

    result = manta._generate_distance_matrix_by_confidence_intervals(viola.TmpVcfForMerge([manta, delly, gridss, lumpy], ['manta', 'delly', 'gridss', 'lumpy']))
    result = viola.merge([manta, delly, lumpy, gridss], mode='confidence_intervals')
    assert result.get_ids() == {'manta_M1', 'manta_MD1', 'manta_ML1', 'manta_MG1', 'manta_MDL1', 'manta_MDG1', 'manta_MLG1', 'manta_MDLG1o', 
    'delly_D1', 'delly_DL1', 'delly_DG1', 'delly_DLG1', 'lumpy_L1', 'lumpy_LG1', 'gridss_G1o'}

def test_merge2():
    manta = viola.read_vcf(StringIO(DATA1), variant_caller='manta')
    delly = viola.read_vcf(StringIO(DATA2), variant_caller='delly')
    result = viola.merge([manta, delly], mode='confidence_intervals', integration=True)
    assert result.sv_count == 4
    assert result.get_ids() == {'manta_M1', 'manta_MD1', 'manta_MDL1', 'delly_D1'}

def test_merge3():
    manta = viola.read_vcf(StringIO(DATA1), variant_caller='manta')
    delly = viola.read_vcf(StringIO(DATA2), variant_caller='delly')
    result = viola.merge([delly, manta], mode='confidence_intervals', integration=True)
    assert result.sv_count == 4
    assert result.get_ids() == {'manta_M1', 'delly_MD1', 'delly_MDL1', 'delly_D1'}