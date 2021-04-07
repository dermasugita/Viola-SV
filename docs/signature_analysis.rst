.. _signature_analysis:


============================
Signature Analysis Tutorial
============================
This is a tutorial for SV signature analysis with Viola.

--------------
Requirements
--------------
If you are new to the Viola package, please refer to the :ref:`Quick Start<quickstart>` first. :ref:`Quick Start<quickstart>` explains basic usage of :doc:`Vcf<reference/api/viola.Vcf>` class, which shares a lot of methods with :doc:`MultiBedpe<reference/api/viola.MultiBedpe>` class.

In this tutorial, we use BEDPE files provided by `ICGC Data Portal`_.
Detailed instruction of data preparation is described below.

.. _ICGC Data Portal: https://dcc.icgc.org/releases/PCAWG/consensus_sv

~~~~~~~~~~~~~~~~~~~~~~~~
PCAWG data preparation
~~~~~~~~~~~~~~~~~~~~~~~~

1. Create new directory named "resources" (Of cource you can name it otherwise) in your working directory and run ``cd resources``.
2. Download ``final_consensus_sv_bedpe_passonly.icgc.public.tgz`` and ``final_consensus_sv_bedpe_passonly.tcga.public.tgz`` to the "resources" directory.
3. Run ``tar zxvf final_consensus_sv_bedpe_passonly.icgc.public.tgz`` and ``tar zxvf final_consensus_sv_bedpe_passonly.icgc.public.tgz``
4. Now you've got the ``icgc/`` directory and the ``tcga/`` directory.
5. Create new directory named "pcawg" in the "resources" directory (This is also customizable).
6. Run ``mv icgc/open/*.gz pcawg`` and ``mv tcga/open/*.gz pcawg``
7. Run ``cd pcawg``, then run ``gunzip *.gz``
8. Now you've got decompressed BEDPE files in the "pcawg" directory.

~~~~~~~~~~~~~~~~~~~~~~~~
Other files preparation
~~~~~~~~~~~~~~~~~~~~~~~~

1. Go back to the "resources" directory.
2. Visit `here`_ and prepare same directories and files (You've already prepared pcawg files by above instruction).
3. Decompress ``replication_timing.bedgraph.gz`` by running ``gunzip replication_timing.bedpe.gz``
4. Go back to your working directory and create new Jupyter Notebook file (You can also use the `ipynb file`_ uploaded on GitHub.).

.. note::
   Example `ipynb file`_ is available on GitHub!

.. _here: https://github.com/dermasugita/Viola-SV/tree/master/examples/demo_sig/resources
.. _ipynb file: https://github.com/dermasugita/Viola-SV/tree/master/examples/demo_sig

.. warning::
   Due to the very large amount of data in PCAWG, the code described below will also take long time (around 50 min on Google Colab). If you want to try it in a shorter time, reduce the number of files in the directory.

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
   pcawg_bedpe=viola.read_bedpe_multi('./resources/pcawg/', exclude_empty_cases=True)
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
2. Reading BEDPE files into viola.MultiBedpe object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   pcawg_bedpe=viola.read_bedpe_multi('./resources/pcawg/', exclude_empty_cases=True)

This code reads all BEDPE files in the ``./resources/pcawg`` directory into MultiBedpe class (See :doc:`MultiBedpe<reference/api/viola.MultiBedpe>`).
Because the BEDPE files created by PCAWG have the ``svclass`` columns, we passed it to the ``svtype_col_name`` argument
Some BEDPE files did not have any SV records. This time we will exclude these by setting ``exclude_empty_cases=True``.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3. Reading BED/BEDGRAPH files for annotation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   bed_fragile = viola.read_bed('./resources/annotation/fragile_site.hg19.bed')
   bed_timing = viola.read_bed('./resources/annotation/replication_timing.bedgraph')

Reading BED and BEDGRAPH files required for custom SV classification. At the moment we do not make a clear distinction between BED files and BEDGRAPH files. This is because only the first four columns of these files are used for annotation purposes in the first place.

``fragile_site.hg19.bed`` is a BED file specifying the known common fragile site (CFS) regions.
``replication_timing.bedgraph`` is a BEDGRAPH file which records the replication timing for each genome coordinate divided into bins.

These files were built according to the `PCAWG paper`_.

.. _PCAWG paper: https://www.nature.com/articles/s41586-019-1913-9#Sec20

~~~~~~~~~~~~~~~~~~~~~
4. Annotating SV
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
5. Get Average values of Replication Timing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python

   pcawg_bedpe.calculate_info('(${timingleft} + ${timingright}) / 2', 'timing')

To get representative values of replication timing for each SV breakpoints, we decided to take mean values of two breakends.
This code adds new INFO named 'timing' by calculating mean values of 'timingleft' and 'timingright'.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
6. Classify SV and Generate Feature Matrix
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


---------------------------
Signature Extraction
---------------------------

Viola offers the function, :doc:`viola.SV_signature_extractor<reference/api/viola.SV_signature_extractor>` which performs non-negative matrix factorization (NMF) and cluster stability evaluation at the same time.

Before diving into the details of :doc:`viola.SV_signature_extractor<reference/api/viola.SV_signature_extractor>`, let's look at an example usage.

.. code-block:: python

   result_silhouette, result_metrics, exposure_matrix, signature_matrix = viola.SV_signature_extractor(
   feature_matrix, n_iter=10, name='testRun', n_components=2, init='nndsvda', solver='mu', beta_loss='kullback-leibler', max_iter=10000, random_state=1
   )

   ######## STDOUT ######## 
   # testRun: finished NMF
   # testRun: finished kmeans clustering
   # testRun: finished all steps
   # Silhouette Score: 0.9994384609504113, kullback-leibler: 98899.5489109343

   # ==================

:doc:`viola.SV_signature_extractor<reference/api/viola.SV_signature_extractor>` outputs four returns. We will explain them one by one.

* result_silhouette:
   This is the score that indicates the stability (reproducibility) of the SV signature. Detailed explanations are given below.
* result_metrics: 
   An error between the original matrix and the product of the factored matrice.
* exposure_matrix: 
   An exposure matrix with (n_samples × n_signatures).
* signature_matrix: 
   A result SV signature matrix with (n_signatures × n_features). Here, ``n_features`` means "the number of custom SV class"