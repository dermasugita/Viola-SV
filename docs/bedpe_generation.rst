.. _bedpe_generation:

.. meta::
   :robots: noindex, nofollow

.. meta::
   :robots: nofollow

Flexible BEDPE generation
=========================

The first step is to create an SgtCore object.

.. ipython:: python
   :okexcept:

   import sgt
   sgt_object = sgt.read_vcf('tutorial.vcf')

Then use ``to_bedpe_like`` method to generate bedpe-formatted pandas DataFrame.
The most simple code is like:

.. ipython:: python
   :okexcept:

   sgt_bedpe_like = sgt_object.to_bedpe_like()
   print(sgt_bedpe_like)

If you want to add "SVLEN" and "CIPOS" fields, run as follows:

.. ipython:: python
   :okexcept:

   sgt_bedpe_like_with_info = sgt_object.to_bedpe_like(custom_infonames=['svlen', 'cipos'])
   print(sgt_bedpe_like_with_info)

To add FORMAT, set ``add_formats = True``:

.. ipython:: python
   :okexcept:

   sgt_bedpe_like_with_format = sgt_object.to_bedpe_like(add_formats=True)
   print(sgt_bedpe_like_with_format)

Do the same for adding FILTER:

.. ipython:: python
   :okexcept:

   sgt_bedpe_like_with_filter = sgt_object.to_bedpe_like(add_filters=True)
   print(sgt_bedpe_like_with_filter)