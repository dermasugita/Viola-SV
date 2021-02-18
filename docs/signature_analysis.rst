.. _signature_analysis:

.. meta::
   :robots: noindex

.. meta::
   :robots: nofollow

============================
Signature Analysis Tutorial
============================
This is a tutorial for signature analysis with PySgt.

--------------
Requirements
--------------
In this tutorial, we use BEDPE files provided by `ICGC Data Portal`_.

.. _ICGC Data Portal: https://dcc.icgc.org/releases/PCAWG/consensus_sv

Before proceeding to the next step, download these files and decompress them as follows.

1. Create new directory named "resources" (Of cource you can name it otherwise) and run ``cd resources``.
2. Download ``final_consensus_sv_bedpe_passonly.icgc.public.tgz`` and ``final_consensus_sv_bedpe_passonly.tcga.public.tgz`` to the "resources" directory.
3. Run ``tar zxvf final_consensus_sv_bedpe_passonly.icgc.public.tgz`` and ``tar zxvf final_consensus_sv_bedpe_passonly.icgc.public.tgz``
4. Now you've got the ``icgc/`` directory and the ``tcga/`` directory.
5. Create new directory named "pcawg" in the "resources" directory (This is also customizable).
6. Run ``mv icgc/open/*.gz pcawg`` and ``mv tcga/open/*.gz pcawg``
7. Run ``cd pcawg``, then run ``gunzip *.gz``
8. Now you've got decompressed BEDPE files in the "pcawg" directory.



----------------
Generate Matrix
----------------
A matrix generation step.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Read all BEDPE files into MultiBedpe class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, read all BEDPE files in the 'pcawg' directory into MultiBedpe class (See :ref:`MultiBedpe<multi_bedpe>`).

.. ipython:: python
   :okexcept:

   import sgt

.. ipython:: python
   :okexcept:

   # bedpe_all = sgt.read_bedpe_multi('./resources/pcawg', svtype_col_name='svclass')
   # bedpe_all

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2. Define your own SV type classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is one of the most unique feature of PySgt.

In this tutorial, we are going to classify SV by:

* Common fragile site
* Classify the SVs into DEL, DUP, INV, and TRA, and then subclassify them according to the size of their microhomology.

.. ipython:: python
   :okexcept:

   # fasta = sgt.read_fasta("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz")
   pass



---------------------
Matrix Factorization
---------------------
A NMF step.

------------------
Effect Size Study
------------------
A effect size study step.