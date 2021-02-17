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

---------------------
Matrix Factorization
---------------------
A NMF step.

------------------
Effect Size Study
------------------
A effect size study step.