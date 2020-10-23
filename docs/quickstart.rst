Quick Start
===========
This is a short examples how to analyse vcf files.

When you read this page, you can

* Import PySgt package to your script.
* Read vcf file to create PySgt core object.
* Convert the vcf file into a bedpe file, with any features you want such as INFO fields, and FORMAT fields.
* Filter the SVs by any features, including genomic positions, INFO fields, and FORMAT fields.

Import PySgt to your script
---------------------------
How to import PySgt package is as follow:

.. ipython:: python
   :okexcept:

   import sgt

Read the vcf file
-----------------
In order to use PySgt's features, you must first import a vcf file to create an SgtCore object.

.. ipython:: python
   :okexcept:

   sgt_object = sgt.read_vcf('tutorial.vcf')

Now you're ready to perform a number of functions that the PySgt package has.

Viewing the contents of SgtCore object
--------------------------------------
You may want to view the content of SgtCore objects after creating them.
To do so, you have three options.

1) Simply print the object.

.. ipython:: python
   :okexcept:

   print(sgt_object)

2) Generate bedpe like pandas.DataFrame.

.. ipython:: python
   
   sgt_bedpe_like = sgt_object.to_bedpe_like()
   print(sgt_bedpe_like)

If you want to add "SVLEN" and "CIPOS" fields, which recorded in original vcf file, run as follows:

.. ipython:: python
   :okexcept:

   sgt_bedpe_like_with_info = sgt_object.to_bedpe_like(custom_infonames=['svlen', 'cipos'])
   print(sgt_bedpe_like_with_info)