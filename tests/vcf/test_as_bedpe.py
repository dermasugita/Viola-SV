import viola
import pandas as pd
from io import StringIO

data = """##fileformat=VCFv4.1
##contig=<ID=chr1,length=195471971>
##contig=<ID=chr2,length=182113224>
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
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
chr1	11	test1	N	<DEL>	.	MinSomaticScore	IMPRECISE;SVTYPE=DEL;SVLEN=-9;END=20;CIPOS=-51,52;CIEND=-51,52	PR:SR	21,0:10,0	43,4:15,3
chr2	10	test2	N	<DUP:TANDEM>	.	PASS	IMPRECISE;SVTYPE=DUP;SVLEN=10;END=20;CIPOS=-27,27;CIEND=-69,70	PR	14,0	38,10
chr2	10	test3	N	<INV>	.	MinSomaticScore;MaxMQ0Frac	SVTYPE=INV;SVLEN=10;END=20;CIPOS=-51,51;CIEND=-89,90;INV5	PR	24,0	35,5
chr1	10	test4_1	N	N[chr2:10[	.	MinSomaticScore	IMPRECISE;SVTYPE=BND;CIPOS=-100,100;MATEID=test4_2	PR	30,1	68,9
chr2	10	test4_2	N	]chr1:10]N	.	MinSomaticScore	IMPRECISE;SVTYPE=BND;CIPOS=-100,100;MATEID=test4_1	PR	30,1	68,9
"""

def test_as_bedpe():
    vcf = viola.read_vcf(StringIO(data), variant_caller='manta', patient_name='patient1')
    bedpe = vcf.as_bedpe()
    df_positions_expected = vcf.get_table('positions')
    df_positions = bedpe.get_table('positions')
    pd.testing.assert_frame_equal(df_positions_expected, df_positions)
    ls_infonames = vcf._ls_infokeys
    for info in ls_infonames:
        df_expected = vcf.get_table(info)
        df = bedpe.get_table(info)
        pd.testing.assert_frame_equal(df_expected, df)

    