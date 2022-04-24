.. _release_notes_1_0_0_stable:

====================================
Viola 1.0.0-stable Release Notes
====================================

---------------
New Features
---------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Class conversion from Vcf to Bedpe
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We implemented Vcf to Bedpe converter :doc:`Vcf.as_bedpe<../reference/api/viola.Vcf.as_bedpe>` as a method of :doc:`Vcf<../reference/api/viola.Vcf>` class.


   >>> vcf = viola.read_vcf('https://raw.githubusercontent.com/dermasugita/ViolaDocs/main/docs/html/_static/tutorial.manta.vcf', patient_name='patient1')
   >>> bedpe = vcf.as_bedpe()
   >>> bedpe
   INFO=imprecise,svtype,svlen,end,cipos,ciend,cigar,mateid,event,homlen,homseq,svinslen,svinsseq,left_svinsseq,right_svinsseq,contig,bnd_depth,mate_bnd_depth,somatic,somaticscore,junction_somaticscore,inv3,inv5
   Documentation of Bedpe object ==> https://dermasugita.github.io/ViolaDocs/docs/html/reference/bedpe.html
           id              be1              be2 strand  qual svtype
   0    test1    chr1:82550461    chr1:82554226     +-  None    DEL
   1    test2    chr1:22814217    chr1:92581132     --  None    INV
   2    test3    chr1:60567906    chr1:60675941     +-  None    DEL
   3    test4    chr1:69583190    chr1:69590948     +-  None    DEL
   4    test5  chr11:104534877  chr11:104536574     +-  None    DEL
   5  test6_1  chr11:111134697   chr17:26470495     +-  None    BND
   6  test6_2   chr17:26470495  chr11:111134697     -+  None    BND

Note that this process is lossy because only the "positions" table and INFO tables are inherited by :doc:`Bedpe<../reference/api/viola.Bedpe>` class.

The similar method :doc:`MultiVcf.as_bedpe_multi<../reference/api/viola.MultiVcf.as_bedpe_multi>` has been also added.


   >>> multi_vcf = viola.read_vcf_multi('/path/to/the/multivcf_dir/', variant_caller='manta')
   >>> multi_bedpe = multi_vcf.as_bedpe_multi()

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Concatenation of MultiBedpe objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Combining data from multiple research projects is effective option for bioinformatics analysis, especially when the size of the data at hand is small.  
Consider the case where you want to mix VCF data from your own research with BEDPE files from PCAWG project.

.. code-block::

   File tree for this section.
   .
   ├── pcawg/
   ├── your_project/
   └── this_script.py
   
First, read the pcawg-BEDPE and your_project-VCF using :doc:`read_bedpe_multi<../reference/api/viola.read_bedpe_multi>` and :doc:`read_vcf_multi<../reference/api/viola.read_vcf_multi>` respectively.

   >>> multi_vcf_your_project = viola.read_vcf_multi('./your_project/', exclude_empty_cases=True)
   >>> multi_bedpe_pcawg = viola.read_bedpe_multi('./pcawg/', svtype_col_name='svclass', exclude_empty_cases=True)

In this update, we've implemented :doc:`viola.concat<../reference/api/viola.concat>` function that concatenate more than one MultiBedpe or Bedpe object into a single MultiBedpe object.  

Since concatenating Vcf-relating objects is not available at this time due to technical reason, here we are going to convert :doc:`MultiVcf<../reference/api/viola.MultiVcf>` object into :doc:`MultiBedpe<../reference/api/viola.MultiBedpe>` object and then concatenate it with pcawg :doc:`MultiBedpe<../reference/api/viola.MultiBedpe>` object.

   >>> multi_bedpe_your_project = multi_vcf_your_project.as_bedpe_multi()
   >>> concat_bedpe = viola.concat([multi_bedpe_your_project, multi_bedpe_pcawg])

Now you've got ``conacate_bedpe`` as the :doc:`MultiBedpe<../reference/api/viola.MultiBedpe>` instance containing all the patients in your project and pcawg!

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A "patient_name" attribute of Bedpe/Vcf classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bedpe/Vcf classes now retain "patient_name" as an attribute.
You can pass the patient name to the read_bedpe and read_vcf functions.


   >>> vcf = read_vcf('/path/to/the/vcf', variant_caller='gridss', patient_name='patient1')
   >>> bedpe = read_bedpe('/path/to/the/bedpe', patient_name='case1')
   >>> bedpe.patient_name
   'patient1'

.. warning::

   Passing NoneType to the patient_name argument is deprecated though the default value of this argument is None.

---------------------
Python 3.9.0 support
---------------------

Python 3.9.0 is now officially supported!

---------------
Fixed Features
---------------

- Fixed feature in :doc:`viola.Vcf.breadend2breakpoint<../reference/api/viola.Vcf.breakend2breakpoint>` not calculating SVLEN (`GH72`_).
- Fixed feature in :doc:`viola.MultiVcf.view<../reference/api/viola.MultiVcf.view>` and :doc:`viola.MultiVcf.view<../reference/api/viola.MultiVcf.view>` showing wrong documentation url.

.. _GH72: https://github.com/dermasugita/Viola-SV/issues/72

---------------
Deprecations
---------------

- Passing NoneType to the 'patient_name' argument of :doc:`viola.read_bedpe<../reference/api/viola.read_bedpe>` or :doc:`viola.read_vcf<../reference/api/viola.read_vcf>` is deprecated.