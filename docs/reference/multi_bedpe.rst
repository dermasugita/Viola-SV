.. _multi_bedpe:


===========
MultiBedpe
===========


.. currentmodule:: viola

Constructor
-----------
.. autosummary::
   :toctree: api/
   
   MultiBedpe


Attributes
----------

.. autosummary::
   :toctree: api/

   MultiBedpe.sv_count
   MultiBedpe.table_list
   MultiBedpe.ids

Accessing the Tables
--------------------

.. autosummary::
   :toctree: api/

   MultiBedpe.get_table
   MultiBedpe.to_bedpe_like

Filtering
---------

.. autosummary::
   :toctree: api/

   MultiBedpe.filter
   MultiBedpe.filter_by_id

Adding Informations
--------------------

.. autosummary::
   :toctree: api/

   MultiBedpe.add_info_table

Annotation
--------------------

.. autosummary::
   :toctree: api/

   MultiBedpe.annotate_bed
   MultiBedpe.get_microhomology

Representation
--------------
.. autosummary::
   :toctree: api/

   MultiBedpe.view

Signature_analysis
--------------------

.. autosummary::
   :toctree: api/

   MultiBedpe.classify_manual_svtype
   MultiBedpe.get_feature_count_as_data_frame