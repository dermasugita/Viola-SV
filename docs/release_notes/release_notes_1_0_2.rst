
.. _release_notes_1_0_2:

====================================
Viola 1.0.2 Release Notes
====================================

---------------
Fixed Features
---------------

- Fixed feature in :doc:`viola.Bedpe.filter<../reference/api/viola.Bedpe.filter>` incorrectly raising exception if even one of the condition is NoneType.
- Fixed feature in :doc:`viola.Vcf.filter<../reference/api/viola.Vcf.filter>`  incorrectly raising exception if even one of the condition is NoneType.
- Fixed feature in :doc:`viola.read_bedpe<../reference/api/viola.read_bedpe>` not generating cipos and ciend table on empty BEDPE.