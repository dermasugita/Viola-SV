
import viola
import sys, os
import numpy as np
import pandas as pd
from io import StringIO
lumpy_HEADER = """##fileformat=VCFv4.2
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

lumpy_body_inv = """chr10	39984191	43	N	<INV>	.	.	SVTYPE=INV;STRANDS=++:4,--:2;SVLEN=176;END=39984367;CIPOS=-8,11;CIEND=-17,31;CIPOS95=-3,6;CIEND95=-1,15;IMPRECISE;SU=6;PE=6;SR=0	GT:SU:PE:SR	./.:3:3:0	./.:3:3:0"""
lumpy_body_all = """chr10	39984191	43	N	<INV>	.	.	SVTYPE=INV;STRANDS=++:4,--:2;SVLEN=176;END=39984367;CIPOS=-8,11;CIEND=-17,31;CIPOS95=-3,6;CIEND95=-1,15;IMPRECISE;SU=6;PE=6;SR=0	GT:SU:PE:SR	./.:3:3:0	./.:3:3:0
chr10	121422696	788	N	<DEL>	.	.	SVTYPE=DEL;STRANDS=+-:8;SVLEN=-21;END=121422717;CIPOS=-9,8;CIEND=-9,8;CIPOS95=0,0;CIEND95=0,0;SU=8;PE=0;SR=8	GT:SU:PE:SR	./.:8:0:8	./.:0:0:0
chr14	45812810	7294_1	N	[chr14:65332314[N	.	.	SVTYPE=BND;STRANDS=--:10;EVENT=7294;MATEID=7294_2;CIPOS=-7,6;CIEND=-7,6;CIPOS95=0,0;CIEND95=0,0;SU=10;PE=0;SR=10	GT:SU:PE:SR	./.:10:0:10	./.:0:0:0
chr14	65332314	7294_2	N	[chr14:45812810[N	.	.	SVTYPE=BND;STRANDS=--:10;SECONDARY;EVENT=7294;MATEID=7294_1;CIPOS=-7,6;CIEND=-7,6;CIPOS95=0,0;CIEND95=0,0;SU=10;PE=0;SR=10	GT:SU:PE:SR	./.:10:0:10	./.:0:0:0
chr15	60882465	8801	N	<DUP>	.	.	SVTYPE=DUP;STRANDS=-+:41;SVLEN=3099777;END=63982242;CIPOS=-3,1;CIEND=-3,2;CIPOS95=0,0;CIEND95=0,0;SU=41;PE=19;SR=22	GT:SU:PE:SR	./.:41:19:22	./.:0:0:0
"""

lumpy_body_inv_buf = StringIO(lumpy_HEADER + lumpy_body_inv)
lumpy_body_all_buf1 = StringIO(lumpy_HEADER + lumpy_body_all)
lumpy_body_all_buf2 = StringIO(lumpy_HEADER + lumpy_body_all)

class TestReadVcfLumpy:
    vcf_lumpy_inv = viola.read_vcf2(lumpy_body_inv_buf, variant_caller='lumpy', patient_name='patient1')
    vcf_lumpy_all1 = viola.read_vcf(lumpy_body_all_buf1, variant_caller='lumpy', patient_name='patient1')
    vcf_lumpy_all2 = viola.read_vcf2(lumpy_body_all_buf2, variant_caller='lumpy', patient_name='patient1')
    def test_compatibility(self):
        viola.testing.assert_vcf_equal(self.vcf_lumpy_all1, self.vcf_lumpy_all2)

    def test_read_vcf_lumpy_inv(self):
        vcf = self.vcf_lumpy_inv

        # assert positions table equal
        df_pos_expected = pd.read_csv(StringIO(
            """id,chrom1,pos1,chrom2,pos2,strand1,strand2,qual,ref,alt,svtype
43_1,chr10,39984191,chr10,39984367,+,+,None,N,<INV>,INV
43_2,chr10,39984192,chr10,39984367,-,-,None,N,<INV>,INV"""
        ))
        df_pos_expected['qual'] = None
        pd.testing.assert_frame_equal(vcf.positions, df_pos_expected)

        # assert strands table equal
        df_info_strands_expected = pd.read_csv(StringIO(
            """id,value_idx,strands
43_1,0,++:4
43_2,0,--:2"""
        ))
        pd.testing.assert_frame_equal(vcf.strands, df_info_strands_expected)
    
        # assert su table equal
        df_info_su_expected = pd.read_csv(StringIO(
            """id,value_idx,su
43_1,0,4
43_2,0,2"""
        ))
        pd.testing.assert_frame_equal(vcf.su, df_info_su_expected)

        # assert suorg table equal
        df_info_suorg_expected = pd.read_csv(StringIO(
            """id,value_idx,suorg
43_1,0,6
43_2,0,6"""
        ))
        pd.testing.assert_frame_equal(vcf.suorg, df_info_suorg_expected)