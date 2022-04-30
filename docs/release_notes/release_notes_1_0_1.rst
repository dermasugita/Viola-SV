.. _release_notes_1_0_1:

====================================
Viola 1.0.1 Release Notes
====================================
---------------
Fixed Features
---------------

- Fixed feature in :doc:`viola.read_vcf<../reference/api/viola.read_vcf>` where patient_name were not reflected in the Vcf class.
- Fixed feature in :doc:`viola.read_bedpe<../reference/api/viola.read_bedpe>` where the confidence intervals were shifted by 1 to the upstream side (5' end).