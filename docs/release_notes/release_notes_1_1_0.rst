.. _release_notes_1_1_0:

====================================
Viola 1.1.0 Release Notes
====================================
---------------
Fixed Features
---------------

- Changed data type of length columns of the ``contigs_meta`` table in :doc:`viola.Vcf<../reference/api/viola.Vcf>` class from 'object' into 'int' (except Lumpy).
- Fixed feature in :doc:`viola.Bedpe.filter_by_id<../reference/api/viola.Bedpe.filter_by_id>` not inherit patient name.
- Fixed feature in :doc:`viola.Vcf.filter_by_id<../reference/api/viola.Vcf.filter_by_id>` not inherit patient name.
- Fixed feature in :doc:`viola.Vcf.copy<../reference/api/viola.Vcf.copy>` not copy patient name.
- Fixed feature in :doc:`viola.Vcf.drop_by_id<../reference/api/viola.Vcf.drop_by_id>` not inherit patient name.
- Fixed feature in :doc:`viola.Bedpe.filter_by_id<../reference/api/viola.Bedpe.filter_by_id>` wrongly raising exception when filtering boolean INFO.