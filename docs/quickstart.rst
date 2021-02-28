.. _quickstart:

.. meta::
   :robots: noindex

.. meta::
   :robots: nofollow

===========
Quick Start
===========
This is a short example how to analyse vcf files.

After you have read this page, you can

* Import PySgt package to your script.
* Read vcf file to create PySgt core object.
* Convert the vcf file into a bedpe file, with any features you want such as INFO fields, and FORMAT fields.
* Filter the SVs by any features, including genomic positions, INFO fields, and FORMAT fields.

If you want to learn how to analyse SV signature with this package, see `Signature Analysis Tutorial<signature_analysis>`.

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
In order to use PySgt's features, you must first run :doc:`read_vcf<reference/api/sgt.io.parser.read_vcf>` to create an :ref:`Vcf<vcf>` object.

.. ipython:: python
   :okexcept:

   url = 'https://dermasugita.github.io/PySgtDocs/docs/html/_static/tutorial.vcf'
   sgt_object = sgt.read_vcf(url) # filepath, url, and file-like object are acceptable.

Now you're ready to perform a number of functions that the PySgt package has.

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

   # The column name 'be' stands for 'breakend'.
   print(sgt_object)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2) Generate bedpe like pandas.DataFrame.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
   
   sgt_bedpe_like = sgt_object.to_bedpe_like()
   print(sgt_bedpe_like)

The way to add INFO, FILTER, and FORMAT to this bedpe-like DataFrame is explained :ref:`here<bedpe_generation>` (no link for now).

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

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1) Filter with queries using filter method 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PySgt has a query system that is easy to understand.

First, let's look at a couple of examples.


**a. Filter with SVTYPE of the INFO field.**
``syntax: "<INFO name> [<value indexer>] <operator> <value>"``

- <value indexer> is optional. 
- The <value indexer> is a 0-origin indexer which allows you to specify which of the comma-separated INFOs, such as CIPOS, should be filtered.
- The following syntax can also be used for other INFO.

.. ipython:: python
   
   # filter with svtype.
   query1_1 = 'svtype == DEL'
   sgt_object.filter(query1_1)

**b. Filter with genomic coordinates.**
``syntax: "<'be1'|'be2'> <chromosome>[:[<start position>]-[<end position>]]"``

- 'be' stands for 'breakend'.
- If you skip <start position> with the minus sign kept, you can get all SV record younger than <end position>, and vice versa if you skip <end position>.
- Note that <start position>-<end position> specifies genomic coordinates with left-closed, right-open interval, that is, [<start position>, <end position>).

.. ipython:: python
   :okexcept:

   # filter with genomic coordinates.
   query1_2 = 'be1 chr11'
   query1_3 = 'be2 chr1:69583189-'
   sgt_object.filter(query1_2)
   sgt_object.filter(query1_3)

**c. Filter with FORMAT table**
``syntax: "<sample name> <FORMAT name> [<FORMAT indexer>] <operator> <value>``

- FORMAT indexer is optional. It is not required when the FORMAT isn't separated by commas.
- FORMAT indexer is 0-origin. Default value is 0.

.. ipython:: python
   :okexcept:

   query1_4 = 'sample1_T PR 1 > 5'
   sgt_object1_4 = sgt_object.filter(query1_4)
   result1_4 = sgt_object1_4.to_bedpe_like(add_formats=True)
   print(result1_4)

**d. Query can be a list**

.. ipython:: python
   
   query2_1 = 'svlen < -4000'
   query2_2 = 'svlen > -10000'
   sgt_object2 = sgt_object.filter([query2_1, query2_2], query_logic='and')
   result2 = sgt_object2.to_bedpe_like(custom_infonames=['svtype', 'svlen'])
   print(result2)
