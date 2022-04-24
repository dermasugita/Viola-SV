import viola
import sys, os
import numpy as np
import pandas as pd
from io import StringIO

HEADER = """
##fileformat=VCFv4.2
##gridssVersion=2.10.2-gridss
##contig=<ID=chr15,length=104043685>
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CIRPOS,Number=2,Type=Integer,Description="Confidence interval around remote breakend POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=VF,Number=1,Type=Integer,Description="Count of fragments supporting the variant breakpoint allele and not the reference allele.">
##FORMAT=<ID=ASQ,Number=1,Type=Float,Description="Pro-rata quality score contribution of assemblies supporting breakpoint">
##FORMAT=<ID=ASRP,Number=1,Type=Integer,Description="Count of read pairs incorporated into any breakpoint assembly">
##FORMAT=<ID=ASSR,Number=1,Type=Integer,Description="Count of split, soft clipped or indel-containing reads incorporated into any breakpoint assemblies">
##FORMAT=<ID=BANRP,Number=1,Type=Integer,Description="Count of read pairs at this breakend assembled into a contig that does not support the breakpoint.">
##FORMAT=<ID=BANRPQ,Number=1,Type=Float,Description="Quality score of read pairs at this breakend assembled into a contig that does not support the breakpoint.">
##FORMAT=<ID=BANSR,Number=1,Type=Integer,Description="Count of split reads at this breakend assembled into a contig that does not support the breakpoint.">
##FORMAT=<ID=BANSRQ,Number=1,Type=Float,Description="Quality score of split reads at this breakend assembled into a contig that does not support the breakpoint.">
##FORMAT=<ID=BAQ,Number=1,Type=Float,Description="Pro-rata quality score contribution of assemblies supporting just local breakend">
##FORMAT=<ID=BASRP,Number=1,Type=Integer,Description="Count of read pairs incorporated into any breakend assembly">
##FORMAT=<ID=BASSR,Number=1,Type=Integer,Description="Count of split, soft clipped or indel-containing reads incorporated into any breakend assemblies">
##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Quality score of breakend evidence after evidence reallocation">
##FORMAT=<ID=BSC,Number=1,Type=Integer,Description="Count of soft clips supporting just local breakend per category">
##FORMAT=<ID=BSCQ,Number=1,Type=Float,Description="Quality score of soft clips supporting just local breakend per category">
##FORMAT=<ID=BUM,Number=1,Type=Integer,Description="Count of read pairs (with one read unmapped) supporting just local breakend per category">
##FORMAT=<ID=BUMQ,Number=1,Type=Float,Description="Quality score of read pairs (with one read unmapped) supporting just local breakend per category">
##FORMAT=<ID=BVF,Number=1,Type=Integer,Description="Count of fragments providing breakend for the variant allele.">
##FORMAT=<ID=CASQ,Number=1,Type=Float,Description="Pro-rata quality score of complex compound breakpoint assemblies supporting breakpoint from elsewhere">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=IC,Number=1,Type=Integer,Description="Count of read indels supporting breakpoint per category">
##FORMAT=<ID=IQ,Number=1,Type=Float,Description="Quality score of read indels supporting breakpoint per category">
##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Quality score of breakend evidence after evidence reallocation">
##FORMAT=<ID=RASQ,Number=1,Type=Float,Description="Pro-rata quality score contribution of assemblies supporting breakpoint from remote breakend">
##FORMAT=<ID=REF,Number=1,Type=Integer,Description="Count of reads mapping across this breakend">
##FORMAT=<ID=REFPAIR,Number=1,Type=Integer,Description="Count of reference read pairs spanning this breakend supporting the reference allele">
##FORMAT=<ID=RF,Number=1,Type=Integer,Description="Reference fragments. Count of fragments supporting the reference allele and not the variant allele.">
##FORMAT=<ID=RP,Number=1,Type=Integer,Description="Count of read pairs supporting breakpoint per category">
##FORMAT=<ID=RPQ,Number=1,Type=Float,Description="Quality score of read pairs supporting breakpoint per category">
##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Count of split reads supporting breakpoint per category">
##FORMAT=<ID=SRQ,Number=1,Type=Float,Description="Quality score of split reads supporting breakpoint per category">
##FORMAT=<ID=VF,Number=1,Type=Integer,Description="Count of fragments supporting the variant breakpoint allele and not the reference allele.">
##FILTER=<ID=ASSEMBLY_BIAS,Description="Mismatch between number of directly supporting reads and reads supporting via assembly.">
##FILTER=<ID=ASSEMBLY_ONLY,Description="Variant is supported only by assembly evidence.">
##FILTER=<ID=ASSEMBLY_TOO_FEW_READ,Description="Not enough reads contribute to this assembly as specified by 'assembly.minReads'">
##FILTER=<ID=ASSEMBLY_TOO_SHORT,Description="This assembly is shorter than a read length">
##FILTER=<ID=INSUFFICIENT_SUPPORT,Description="Does not reach the required threshold quality for calling as specified by 'variantcalling.minScore'">
##FILTER=<ID=LOW_QUAL,Description="Low quality call as specified by 'variantcalling.lowQuality'">
##FILTER=<ID=NO_ASSEMBLY,Description="No assembly supporting this variant could be found.">
##FILTER=<ID=NO_RP,Description="Variant does not have any direct read pair support.">
##FILTER=<ID=NO_SR,Description="Variant does not have any direct split read support.">
##FILTER=<ID=REF,Description="Breakpoint corresponds to reference allele">
##FILTER=<ID=SINGLE_ASSEMBLY,Description="Only one side of the breakpoint could be assembled.">
##FILTER=<ID=SINGLE_SUPPORT,Description="Supported by fewer than 'variantcalling.minReads' fragments">
##FILTER=<ID=SMALL_EVENT,Description="Event size is smaller than the minimum reportable size specified by 'variantcalling.minSize'">
##ALT=<ID=INV,Description="Inversion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mouse01_N	mouse01_T
"""
body = """chr15	33050573	single_test	T	.T	1008.31	PASS	CIPOS=-1,2;CIRPOS=-2,1;EVENT=gridss75ff_2435;HOMLEN=3;HOMSEQ=TCA;SVTYPE=BND;VF=10	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:23:0:0:0.00:0:0.00:0	.:330.10:10:23:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:1008.31:348.12:36:2:5:92.98:11:237.12:10
chr15	63184301	gridss68ff_1141o	G	G]chr15:64186262]	3801.43	PASS	EVENT=gridss68ff_1141;MATEID=gridss68ff_1141h;SVTYPE=BND;VF=48	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0.0:0:0:0:0.0:0:0.0:0.0:0:0:0.0:0:0.0:0:0.0:0:0.0:0:0.0:0.0:0.0:24:1:0:0.0:0:0.0:0	.:1323.15:38:88:0:0.0:0:0.0:0.0:0:0:174.1:9:174.1:0:0.0:0:0.0:0:0.0:3801.43:1286.81:67:0:19:353.32:38:838.16:48
chr15	64186262	gridss68ff_1141h	T	T]chr15:63184301]	3801.43	PASS	EVENT=gridss68ff_1141;MATEID=gridss68ff_1141o;SVTYPE=BND;VF=48	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0.0:0:0:0:0.0:0:0.0:0.0:0:0:0.0:0:0.0:0:0.0:0:0.0:0:0.0:0.0:0.0:51:1:0:0.0:0:0.0:0	.:1286.81:38:88:0:0.0:0:0.0:0.0:0:0:95.33:5:95.33:0:0.0:0:0.0:0:0.0:3801.43:1323.15:95:0:19:353.32:38:838.16:48
"""

def test_single():
    vcf = viola.read_vcf(StringIO(HEADER + body), variant_caller='gridss', patient_name='patient1')
    vcf_single = vcf.idx['single_test']
    print(vcf_single.positions)