.. _usecases:

==================================
Use Cases for Bioinformatics
==================================
Here we exemplify the practical use cases of Viola-SV for bioinformatics analysis.

-------------------------------------------------------------------------------
Remove blacklist regions from Vcf object using a blacklist BED file
-------------------------------------------------------------------------------
Some SV callers, such as Gridss, have the ability to exclude blacklist regions by default, while others, such as Manta and Lumpy, do not.
If you are a fan of Manta and want to exclude blacklist regions as well, you can use Viola to do the following.


~~~~~~~~~~~~~~~~~~~~~~~~
Preparations
~~~~~~~~~~~~~~~~~~~~~~~~
You need following two files.

- a VCF file from Manta
- a BED file of blacklist region

In this article, we employ "`DAC Exclusion List Regions`_" provided by ENCODE project for the blacklist BED file of mm10.

For Manta VCF, we use ``tutorial.manta.vcf`` that is also used in :ref`Quick Start<quickstart>` article.

.. _DAC Exclusion List Regions: https://www.encodeproject.org/annotations/ENCSR636HFF/

First, download the blacklist BED file.

.. code:: bash

    wget "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"

Next, import viola and read the Manta VCF file and blacklist BED file as Viola objects.

    >>> import viola
    >>> manta_url = 'https://raw.githubusercontent.com/dermasugita/Viola-SV/master/docs/_static/tutorial.manta.vcf'
    >>> manta_vcf = viola.read_vcf(manta_url, variant_caller='manta', patient_name='filter_test')
    >>> blacklist_bed = viola.read_bed('ENCFF547MET.bed.gz')


