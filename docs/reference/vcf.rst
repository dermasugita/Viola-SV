.. _vcf:

Vcf
=======

.. currentmodule:: viola

Constructor
-----------
.. autosummary::
   :toctree: api/
   
   Vcf

Attributes
----------
.. autosummary::
   :toctree: api/

   Vcf.sv_count
   Vcf.table_list
   Vcf.ids


Annotation
--------------------
.. autosummary::
   :toctree: api/

   Vcf.annotate_bed
   Vcf.get_microhomology

Conversion
-----------
.. autosummary::
   :toctree: api/

   Vcf.copy
   Vcf.to_bedpe_like
   Vcf.as_bedpe
   Vcf.to_vcf_like
   Vcf.breakend2breakpoint

Export
--------
.. autosummary::
   :toctree: api/

   Vcf.to_bedpe
   Vcf.to_vcf

Filtering
---------
.. autosummary::
   :toctree: api/

   Vcf.filter
   Vcf.filter_by_id

Managing Tables
--------------------

.. autosummary::
   :toctree: api/

   Vcf.add_info_table
   Vcf.append_infos
   Vcf.append_formats
   Vcf.append_filters
   Vcf.get_table
   Vcf.remove_info_table
   Vcf.replace_table
   Vcf.set_value_for_info_by_id

Merging
----------

.. autosummary::
   :toctree: api/

   Vcf.merge
   Vcf.integrate

Representation
--------------
.. autosummary::
   :toctree: api/

   Vcf.view

Signature_analysis
--------------------

.. autosummary::
   :toctree: api/

   Vcf.classify_manual_svtype
   Vcf.get_feature_count_as_series

SV ID/SV records
---------------------

.. autosummary::
   :toctree: api/

   Vcf.replace_svid
   Vcf.drop_by_id