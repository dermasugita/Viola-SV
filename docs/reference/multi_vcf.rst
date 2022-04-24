.. _multi_vcf:


===========
MultiVcf
===========


.. currentmodule:: viola

Constructor
-----------
.. autosummary::
   :toctree: api/
   
   MultiVcf


Attributes
----------

.. autosummary::
   :toctree: api/

   MultiVcf.sv_count
   MultiVcf.table_list
   MultiVcf.ids

Accessing the Tables
--------------------

.. autosummary::
   :toctree: api/

   MultiVcf.get_table
   MultiVcf.to_bedpe_like

Filtering
---------

.. autosummary::
   :toctree: api/

   MultiVcf.filter
   MultiVcf.filter_by_id

Adding Informations
--------------------

.. autosummary::
   :toctree: api/

   MultiVcf.add_info_table

Annotation
--------------------

.. autosummary::
   :toctree: api/

   MultiVcf.annotate_bed
   MultiVcf.get_microhomology

Conversion
---------------

.. autosummary::
   :toctree: api/

   MultiVcf.as_bedpe_multi

Representation
--------------
.. autosummary::
   :toctree: api/

   MultiVcf.view

Signature_analysis
--------------------

.. autosummary::
   :toctree: api/

   MultiVcf.classify_manual_svtype
   MultiVcf.get_feature_count_as_data_frame