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

   bedpe_all = sgt.read_bedpe_multi('./resources/pcawg', svtype_col_name='svclass')

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2. Define your own SV type classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is one of the most unique feature of PySgt.

In this tutorial, we are going to classify SV by:

* Common fragile site
* Classify the SVs into DEL, DUP, INV, and TRA, and then subclassify them according to the size of their microhomology.

First, let's load the BED file of the fragile sites.

.. ipython:: python
   :okexcept:

   fragile_bed = sgt.read_bed('fragile_site.hg19.bed')

Bed class object has been created with this code.

Then, you can make fragile site annotation on the MultiBedpe class.

.. ipython:: python
   :okexcept:

   bedpe_all.annotate_bed(fragile_bed, 'fragile')
   print(bedpe_all.table_list)

As you see, the new table named 'fragile_left' and 'fragile_right' have been added.

Now, let's move on microhomology length inference.

This process requires hg19 fasta file. To read it, run code shown below to get Fasta class object.

.. ipython:: python
   :okexcept:

   fasta = sgt.read_fasta("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz")

:doc:`get_microhomology<reference/api/sgt.core.db.Bedpe.get_microhomology>` method of Bedpe class is quite useful to append INFO table of microhomology.

.. ipython:: python
   :okexcept:

   bedpe_all.get_microhomology(fasta)
   print(bedpe_all.table_list)

The 'HOMLEN' and 'HOMSEQ' tables have been successfully added!

Now, let's define custom SV types and add new INFO table named 'manual_sv_type'.

To do this, you can use :doc:`classify_manual_svtype<reference/api/sgt.core.cohort.MultiBedpe.classify_manual_svtype>` method.

This function requires at least two arguments:

* The first argument (ls_conditions) is a list of functions that defines custom SV classes.
* The second argument (ls_names) is a list of names of SV classes that correspond to 'ls_conditions'.

The following is the specification of the function that is the element of the argument 'ls_conditions'.

For each custom SV class:

* Takes Bedpe, MultiBedpe, and Vcf objects as arguments
* Returns SV id that satisfies the requirements of the respective custom SV class.

Here are the examples.

.. ipython:: python
   :okexcept:

   def fragile_del(x):
       set1 = set(x.filter('fragileleft == True').ids)
       set2 = set(x.filter('fragileright == True').ids)
       set3 = set(x.filter(['svtype == DEL']).ids)
       return (set1 | set2) & set3
   
   def fragile_dup(x):
       set1 = set(x.filter('fragileleft == True').ids)
       set2 = set(x.filter('fragileright == True').ids)
       set3 = set(x.filter(['svtype == DUP']).ids)
       return (set1 | set2) & set3
   
.. ipython:: python
   :okexcept:
   
   def del_micro1(x):
       ids = x.filter(['svtype == DEL', 'homlen <= 1']).ids
       return set(ids)
   
   def del_micro2(x):
       ids = x.filter(['svtype == DEL']).ids
       return set(ids)
   
   def dup_micro1(x):
       ids = x.filter(['svtype == DUP', 'homlen <= 1']).ids
       return set(ids)
   
   def dup_micro2(x):
       ids = x.filter(['svtype == DUP']).ids
       return set(ids)
   
   def inv_micro1(x):
       ids1 = x.filter(['svtype == t2tINV', 'svtype == h2hINV'], query_logic='or').ids
       ids2 = x.filter(['homlen <= 1']).ids
       return set(ids1) & set(ids2)
   
   def inv_micro2(x):
       ids1 = x.filter(['svtype == t2tINV', 'svtype == h2hINV'], query_logic='or').ids
       return set(ids1)

   def tra_micro1(x):
       ids = x.filter(['svtype == TRA', 'homlen <= 1']).ids
       return set(ids)
   
   def tra_micro2(x):
       ids = x.filter(['svtype == TRA']).ids
       return set(ids)
   


Finally, we can obtain a matrix with custom SV classes as columns and the counts for each patient as rows with following operation.

.. ipython:: python
   :okexcept:

   ls_conditions = [fragile_del, fragile_dup, del_micro1, del_micro2, dup_micro1, dup_micro2, inv_micro1, inv_micro2, tra_micro1, tra_micro2]
   ls_names = ['fragile_del', 'fragile_dup', 'del_micor<=1', 'del_micro>1', 'dup_micor<=1', 'dup_micro>1', 'inv_micor<=1', 'inv_micro>1', 'tra_micor<=1', 'tra_micro>1']
   df = bedpe_all.classify_manual_svtype(ls_conditions, ls_names)
   print(df)


---------------------
Matrix Factorization
---------------------
A NMF step.

------------------
Effect Size Study
------------------
A effect size study step.