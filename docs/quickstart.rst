.. _quickstart:

===========
Quick Start
===========
This is a short example how to manage VCF files.

After you have read this page, you can

* Import Viola package to your script.
* Read VCF files to create viola.Vcf object.
* Convert VCF files into BEDPE files, with any features you wish, such as INFO fields and FORMAT fields.
* Filter SV records by any features, including genomic positions, INFO fields, FILTER fields, and FORMAT fields.

See Also:

*  :ref:`Signature Analysis Tutorial<signature_analysis>`
*  :ref:`API Reference<api_reference>`


In this tutorial, we will only deal with the viola.Vcf class, but many of the features available in this class are also available in the Bedpe, MultiBedpe and MultiVcf classes.

Individual tutorials for classes not covered here are currently under development.

---------------------------
Import Viola to your script
---------------------------

.. ipython:: python
   :okexcept:

   import viola

-------------------------------------
Create Vcf object from VCF file.
-------------------------------------
In order to use Viola's features, you must first run :doc:`read_vcf<reference/api/viola.read_vcf>` to create an :ref:`viola.Vcf<vcf>` object.

.. ipython:: python

   url = 'https://raw.githubusercontent.com/dermasugita/ViolaDocs/main/docs/html/_static/tutorial.manta.vcf'
   vcf = viola.read_vcf(url, variant_caller='manta') # filepath, url, and file-like object are acceptable.

Now you're ready to perform a number of functions that the Viola package has.

To learn detail data structure of viola.Vcf object, see :ref:`Viola's Data Structure<data_structure>`.

-------------------------------------------
Viewing the contents of viola.Vcf object
-------------------------------------------
You may want to view the content of Vcf objects after creating them.
To do so, you have three options.

~~~~~~~~~~~~~~~~~~~~~~~~~~~
1) Simply print the object.
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
   :okexcept:

   # The column name 'be' stands for 'breakend'.
   print(vcf)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2) Generate BEDPE like pandas.DataFrame.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
   :okexcept:
   
   vcf_bedpe_like = vcf.to_bedpe_like()
   print(vcf_bedpe_like)

The way to add INFO, FILTER, and FORMAT to this bedpe-like DataFrame is explained :ref:`here<bedpe_generation>`.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
3) Use the :doc:`get_table()<reference/api/viola.Vcf.get_table>` method to get individual tables composing the Vcf object.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. ipython:: python

   # the returned value is pd.DataFrame
   vcf.get_table('positions')

The names of all the tables in ``vcf`` are stored in the :doc:`table_list<reference/api/viola.Vcf.table_list>` attribute:

.. ipython:: python
   
   vcf.table_list

You can get any table you want.

.. ipython:: python

   vcf.get_table('formats_meta') # get header information of FORMAT field

------------------------
Export as VCF/BEDPE file
------------------------

You can export VCF/BEDPE files by ``to_vcf``/``to_bedpe`` method.

.. code-block:: python

   vcf.to_vcf('/path/to/the/output.vcf')
   vcf.to_bedpe('/path/to/the/output.bedpe')

---------------------
Filter Vcf object
---------------------

Filtering VCF file is an essential step of bioinformatics study.
viola.Vcf object provides an intuitive way to filter SV in almost any item.

You have two options for filtering. 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1) Filter with queries using filter method 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Viola has a query system that is easy to understand.

First, let's look at a couple of examples.


**a. Filter with SVTYPE of the INFO field.**
``syntax: "<INFO name> [<value index>] <operator> <value>"``

- Please do not omit spaces.
- <value indexer> is optional. 
- The <value indexer> is a 0-origin indexer which allows you to specify which of the comma-separated INFOs, such as CIPOS, should be filtered.
- The following syntax can also be used for other INFO.

.. ipython:: python
   
   # filter with svtype.
   vcf.filter('svtype == DEL')

.. ipython:: python

   # example of filtering using <value index>
   # The code below means that if the right value of CIPOS (e.g. the value "20" of CIPOS=-10,20) is 
   # lower than 60, the SV record will be output.
   using_value_idx = vcf.filter('cipos 1 < 60').to_vcf_like()
   print(using_value_idx)
   print(using_value_idx['info'].values)

.. note::
   **What is <value index>?**

   Some INFO and FORMAT entries are separated by commas to store more than one value.
   For example, the ``CIPOS`` and ``CIEND`` entries of INFO field always store two values:

   .. code::

      CIPOS=-1,1;CIEND=0,2

   ``<value index>`` is assigned to such data with 0-origin manner as shown below:

   .. image:: ./_static/quickstart/value_index.png
      :width: 300
   
   | 
   As another example, the values in the FORMAT field of Manta VCF are separated by commas.

   .. image:: ./_static/quickstart/value_index2.png
      :width: 350
   

**b. Filter with genomic coordinates.**
``syntax: "<'be1'|'be2'> <chromosome>[:[<start position>]-[<end position>]]"``

- 'be' stands for 'breakend'.
- If you skip <start position> with the minus sign kept, you can get all SV record younger than <end position>, and vice versa if you skip <end position>.
- Note that <start position>-<end position> specifies genomic coordinates with left-closed, right-open interval, that is, [<start position>, <end position>).

.. ipython:: python
   :okexcept:

   # filter with genomic coordinates.
   vcf.filter('be1 chr11')
   vcf.filter('be1 !chr1')
   vcf.filter('be2 chr1:69583189-')

**c. Filter with FORMAT table**
``syntax: "<sample name> <FORMAT name> [<FORMAT index>] <operator> <value>``

- FORMAT indexer is optional. It is not required when the FORMAT isn't separated by commas.
- FORMAT indexer is 0-origin. Default value is 0.

.. ipython:: python
   :okexcept:

   # The meaning of this code is that if the right value of the PR of sample1_T
   # in the FORMAT field (e.g. the value "3" in PR:SR 6,7:8,9  1,3:5,2) is greater than 5,
   # the SV records will be returned.
   vcf.filter('sample1_T PR 1 > 5').to_bedpe_like(add_formats=True)

**c. Filter with FILTER table**
``syntax: "[!]<FILTER name>"``

- When "!" mark is prepended, the SV records excluding <FILTER name> are returned.

.. ipython:: python
   :okexcept:

   vcf.filter('PASS').to_bedpe_like(add_filters=True)
   vcf.filter('!PASS').to_bedpe_like(add_filters=True)

**d. Query can be a list**

.. ipython:: python
   
   query2_1 = 'svlen < -4000'
   query2_2 = 'svlen > -10000'
   vcf2 = vcf.filter([query2_1, query2_2], query_logic='and')
   result2 = vcf2.to_bedpe_like(custom_infonames=['svtype', 'svlen'])
   print(result2)

You can perform set operations by passing expressions to query_logic.

.. ipython:: python
   
   query3_1 = 'svtype == DEL'
   query3_2 = 'svtype == BND'
   query3_3 = 'somaticscore > 20'

   # (query3_1 or query3_2) and (query3_3)
   vcf3 = vcf.filter([query3_1, query3_2, query3_3], query_logic='(0 | 1) & 2')
   print(vcf3.to_bedpe_like(custom_infonames=['svtype', 'somaticscore']))


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2) Filter with SV ID using filter_by_id method 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Vcf.filter`` is very useful, but in some situation, you may have to filter with much more complex criteria.
In such cases we recommend to use :doc:`Vcf.filter_by_id<reference/api/viola.Vcf.filter_by_id>` method.

Suppose you obtained the list of SV id, such as ``['test2', 'test4']``, as a result of quite complex criterion.
In this case, your code should be:

.. ipython:: python

   vcf.filter_by_id(['test2', 'test4'])

