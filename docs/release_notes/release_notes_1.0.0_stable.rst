# Viola 1.0.0-stable Release Notes

## New Features

### Class conversion from Vcf to Bedpe
We implemented Vcf to Bedpe converter Vcf.as_bedpe as a method of Vcf class.

.. code-block:: python

   import viola
   vcf = viola.read_vcf('https://raw.githubusercontent.com/dermasugita/ViolaDocs/main/docs/html/_static/tutorial.manta.vcf')
   bedpe = vcf.as_bedpe()
   bedpe

   ######## STDOUT ######## 
   # INFO=imprecise,svtype,svlen,end,cipos,ciend,cigar,mateid,event,homlen,homseq,svinslen,svinsseq,left_svinsseq,right_svinsseq,contig,bnd_depth,mate_bnd_depth,somatic,somaticscore,junction_somaticscore,inv3,inv5
   # Documentation of Bedpe object ==> https://dermasugita.github.io/ViolaDocs/docs/html/reference/bedpe.html
   #         id              be1              be2 strand  qual svtype
   # 0    test1    chr1:82550461    chr1:82554226     +-  None    DEL
   # 1    test2    chr1:22814217    chr1:92581132     --  None    INV
   # 2    test3    chr1:60567906    chr1:60675941     +-  None    DEL
   # 3    test4    chr1:69583190    chr1:69590948     +-  None    DEL
   # 4    test5  chr11:104534877  chr11:104536574     +-  None    DEL
   # 5  test6_1  chr11:111134697   chr17:26470495     +-  None    BND
   # 6  test6_2   chr17:26470495  chr11:111134697     -+  None    BND

Note that this process is lossy because only the "positions" table and INFO tables are inherited by Bedpe class.

## Python 3.9.0 support




## Fixed Features

- Fixed feature in viola.Vcf.breadend2breakpoint not calculating SVLEN (`GH72`_).

.. _GH72: https://github.com/dermasugita/Viola-SV/issues/72

## Backwards Incompatible Changes
