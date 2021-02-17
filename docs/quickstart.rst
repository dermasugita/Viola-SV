.. _quickstart:

.. meta::
   :robots: noindex

.. meta::
   :robots: nofollow

===========
Quick Start
===========
This is a short example how to manage VCF files.

After you have read this page, you can

* Import PySgt package to your script.
* Read VCF file to create PySgt core object.
* Convert the VCF file into a BEDPE file, with any features you want such as INFO fields, and FORMAT fields.
* Filter the SVs by any features, including genomic positions, INFO fields, and FORMAT fields.

See Also:

*  :ref:`Signature Analysis Tutorial<signature_analysis>`
*  :ref:`API Reference<api_reference>`



---------------------------
Import PySgt to your script
---------------------------
How to import PySgt package is as follow:

.. ipython:: python
   :okexcept:

   import sgt

-----------------
Read the vcf file
-----------------
In order to use PySgt's features, you must run :doc:`read_vcf<reference/api/sgt.io.parser.read_vcf>` to create an :ref:`Vcf<core>` object first.

.. ipython:: python
   :okexcept:

   url = 'https://dermasugita.github.io/PySgtDocs/docs/html/_static/tutorial.vcf'
   sgt_object = sgt.read_vcf(url) # filepath, url, and file-like object are acceptable.

Now you're ready to perform a number of functions that the PySgt package has.

To learn detail data structure of Vcf object, see :ref:`PySgt's Data Structure<data_structure>`.

--------------------------------------
Viewing the contents of Vcf object
--------------------------------------
You may want to view the content of Vcf objects after creating them.
To do so, you have three options.

~~~~~~~~~~~~~~~~~~~~~~~~~~~
1) Simply print the object.
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
   :okexcept:

   print(sgt_object)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2) Generate BEDPE like pandas.DataFrame.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
   
   sgt_bedpe_like = sgt_object.to_bedpe_like()
   print(sgt_bedpe_like)

The way to add INFO, FILTER, and FORMAT to this bedpe-like DataFrame is explained :ref:`here<bedpe_generation>`.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3) Use the :doc:`get_table()<reference/api/sgt.core.db.Vcf.get_table>` method to get individual tables composing the Vcf object.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. ipython:: python

   sgt_object.get_table('positions')

The names of all the tables in ``sgt_object`` are stored in the :doc:`table_list<reference/api/sgt.core.db.Vcf.table_list>` attribute:

.. ipython:: python
   
   sgt_object.table_list

You can get any table you want.

.. ipython:: python

   sgt_object.get_table('formats_meta') # get header information of FORMAT field

------------------------
Export as VCF/BEDPE file
------------------------

Under Programming

---------------------
Filter Vcf object
---------------------

Filtering vcf file is an essential step of bioinformatics study.
Vcf object provides an intuitive way to filter SV in almost any item.

You have two options for filtering. 

~~~~~~~~~~~~~~~~~~~~~~
1) Filter with queries
~~~~~~~~~~~~~~~~~~~~~~
The query system is somewhat specific to this package, but still easy to understand since they are very simple.

First, let's look at a couple of examples.

.. ipython:: python
   
   query1 = 'svtype == DEL'
   sgt_object_deletion = sgt_object.filter(query1)
   df_out = sgt_object_deletion.to_bedpe_like(custom_infonames=['svtype', 'svlen'])
   print(df_out)

Query can be a list.

.. ipython:: python
   
   query2_1 = 'svlen < -4000'
   query2_2 = 'svlen > -10000'
   sgt_object_filter_len = sgt_object.filter([query2_1, query2_2], query_logic='and')
   df_out = sgt_object_filter_len.to_bedpe_like(custom_infonames=['svtype', 'svlen'])
   print(df_out)
