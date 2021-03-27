.. _filtering:

======================
Filtering Systems 
======================

| Filtering is an important step in bioinformatics.
| Filtering functions in packaged systems need to be both complete enough to handle even the most complex filtering criteria, as well as available enough to allow simple filtering with simple coding.
| This page explains how PySgt achieves these desired functions.

---------------------------------------------
Completeness: ID based filtering system
---------------------------------------------
In the main PySgt object, including Vcf, Bedpe, and MultiBedpe, filtering-associated functions are eventually passed to the :doc:`filter_by_id<../reference/api/viola.Vcf.filter_by_id>` method.
| :doc:`filter_by_id<../reference/api/viola.Vcf.filter_by_id>` method receives a list of ID (SV ID or global ID) and return filtered object.
| Any table in the object that has a column with an SV ID or global ID will also be filtered.