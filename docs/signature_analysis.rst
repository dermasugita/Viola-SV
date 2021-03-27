.. _signature_analysis:

.. meta::
   :robots: noindex

.. meta::
   :robots: nofollow

============================
Signature Analysis Tutorial
============================
This is a tutorial for SV signature analysis with Viola.

--------------
Requirements
--------------
If you are new to the Viola package, please refer to the :ref:`Quick Start<quickstart>` first. :ref:`Quick Start<quickstart>` explains basic usage of :doc:`Vcf<reference/api/viola.Vcf>` class, which shares a lot of methods with :doc:`MultiBedpe<reference/api/viola.MultiBedpe>` class.

In this tutorial, we use BEDPE files provided by `ICGC Data Portal`_.
Before proceeding to the next step, download these files and decompress them as follows.

.. _ICGC Data Portal: https://dcc.icgc.org/releases/PCAWG/consensus_sv

~~~~~~~~~~~~~~~~~~~~~~~~
PCAWG data preparation
~~~~~~~~~~~~~~~~~~~~~~~~

1. Create new directory named "resources" (Of cource you can name it otherwise) and run ``cd resources``.
2. Download ``final_consensus_sv_bedpe_passonly.icgc.public.tgz`` and ``final_consensus_sv_bedpe_passonly.tcga.public.tgz`` to the "resources" directory.
3. Run ``tar zxvf final_consensus_sv_bedpe_passonly.icgc.public.tgz`` and ``tar zxvf final_consensus_sv_bedpe_passonly.icgc.public.tgz``
4. Now you've got the ``icgc/`` directory and the ``tcga/`` directory.
5. Create new directory named "pcawg" in the "resources" directory (This is also customizable).
6. Run ``mv icgc/open/*.gz pcawg`` and ``mv tcga/open/*.gz pcawg``
7. Run ``cd pcawg``, then run ``gunzip *.gz``
8. Now you've got decompressed BEDPE files in the "pcawg" directory.

.. warning::
   Due to the very large amount of data in PCAWG, the code described below will also take some time. If you want to try it in a shorter time, reduce the number of files in the directory.

----------------
Generate Matrix
----------------

Here, we explain how to generate SV feature matrix with custom SV class.
In this tutorial, we are going to classify SV by:

* Classify the SVs into DEL, DUP, INV, and TRA, and then subclassify them according to their size.
* Common fragile sites
* Replication timing

This can be achieved with only eight lines of code.
The actual code is shown below.

.. code-block:: python

   import viola
   pcawg_bedpe=viola.read_bedpe_multi('./resources/pcawg/')
   bed_fragile = viola.read_bed('./resources/annotation/fragile_site.hg19.bed')
   bedgraph_timing = viola.read_bed('./resources/annotation/replication_timing.bedgraph')
   pcawg_bedpe.annotate_bed(bed=bed_fragile, annotation='fragile', how='flag')
   pcawg_bedpe.annotate_bed(bed=bedgraph_timing, annotation='timing', how='value')
   pcawg_bedpe.calculate_info('(${timingleft} + ${timingright}) / 2', 'timing')
   feature_matrix = pcawg_bedpe.classify_manual_svtype(definitions='./resources/definitions/sv_class_definition.txt', return_data_frame=True)

Now let's look at what each line means!

~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Importing Viola
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import viola

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Reading BEDPE files into viola.MultiBedpe object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   pcawg_bedpe=viola.read_bedpe_multi('./resources/pcawg/')

This code reads all BEDPE files in the ``./resources/pcawg`` directory into MultiBedpe class (See :doc:`MultiBedpe<reference/api/viola.MultiBedpe>`).
Because the BEDPE files created by PCAWG have the ``svclass`` columns, we passed it to the ``svtype_col_name`` argument.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2. Reading BED/BEDGRAPH files for annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   bed_fragile = viola.read_bed('./resources/annotation/fragile_site.hg19.bedgraph')
   bed_timing = viola.read_bed('./resources/annotation/replication_timing.bedgraph')

Reading BED and BEDGRAPH files required for custom SV classification. At the moment we do not make a clear distinction between BED files and BEDGRAPH files. This is because only the first four columns of these files are used for annotation purposes in the first place.

``fragile_site.hg19.bed`` is a BED file specifying the known common fragile site (CFS) regions.
``replication_timing.bedgraph`` is a BEDGRAPH file which records the replication timing for each genome coordinate divided into bins.

These files were built according to the `PCAWG paper`_.

.. _PCAWG paper: https://www.nature.com/articles/s41586-019-1913-9#Sec20

~~~~~~~~~~~~~~~~~~~~~
3. Annotating SV
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   pcawg_bedpe.annotate_bed(bed=bed_fragile, annotation='fragile', how='flag')
   pcawg_bedpe.annotate_bed(bed=bedgraph_timing, annotation='timing', how='value')

In this step, we annotate ``pcawg_bedpe`` with the Bed object we've just loaded. 
After annotation, new INFO – 'fragileleft', 'fragileright', 'timingleft', and 'timingright' – will be added.
Because two breakends form a single SV, 'left' and 'right' suffix are added.
When ``how='flag'``, annotate True/False according wether each breakend is in the range in the Bed (4th column of the Bed is ignored).
When ``how='value'``, annotate the value of 4th column of Bed if the breakends hit.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
4. Get Average values of Replication Timing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python

   pcawg_bedpe.calculate_info('(${timingleft} + ${timingright}) / 2', 'timing')

To get representative values of replication timing for each SV breakpoints, we decided to take mean values of two breakends.
This code adds new INFO named 'timing' by calculating mean values of 'timingleft' and 'timingright'.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
5. Classify SV and Generate Feature Matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   feature_matrix = pcawg_bedpe.classify_manual_svtype(definitions='./resources/definitions/sv_class_definition.txt', return_data_frame=True)

Finally, we classified SV according to its type, size, fragile site, and replication timing. Classification criteria are written in ``sv_class_definition.txt``. Syntax of this file is explained below.

If ``return_data_frame=True``, counts of each custom SV class for each patients are returned as pandas.DataFrame.

Now we successfully obtained feature matrix with custom SV classification!

~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Definition File Syntax
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

   name 'At fragile site DEL'
   0 fragileleft == True
   1 fragileright == True
   2 svtype == DEL
   logic (0 | 1) & 2

   name 'At fragile site DUP'
   0 fragileleft == True
   1 fragileright == True
   2 svtype == DUP
   logic (0 | 1) & 2

   name '<50 kb early DEL'
   0 svlen > -50000
   1 timing > 66.65
   2 svtype == DEL
   logic 0 & 1 & 2

   name '<50 kb mid DEL'
   0 svlen > -50000
   1 timing > 33.35
   2 svtype == DEL
   logic 0 & 1 & 2

   name '<50 kb late DEL'
   0 svlen > -50000
   1 svtype == DEL
   logic 0 & 1

This is an example of definition file of custom SV classification.

Each SV class is defined by a syntax like the following:

.. code-block::

   name '<SV class name>'
   0 <condition>
   1 <condition>
   2 <condition>
   ...
   logic <set operation>

- The syntax of <condition> is the same as query passed to Vcf.filter method (See :ref:`Quick Start<quickstart>`).
- The numbers written in the left of each <condition> can be omitted.
- Use numbers correspoinding to each <condition> for the <set operation>

.. note::

   The order of the SV class definition is very important. The ``classify_manual_svtype`` method reads the definition file in order from the top, so that the SV class definitions written higher up in the file take precedence. Thus, in above example, Deletions that both satisfy 'At fragile site DEL' and '<50 kb early DEL' criteria, they are classified as 'At fragile site DEL', not '<50 kb early DEL'.
