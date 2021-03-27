.. _bedpe_generation:

.. meta::
   :robots: noindex
   
.. meta::
   :robots: nofollow

Flexible BEDPE generation
=========================

The first step is to create an viola.Vcf object.

.. ipython:: python
   :okexcept:

   import viola
   vcf = viola.read_vcf('https://raw.githubusercontent.com/dermasugita/ViolaDocs/main/docs/html/_static/tutorial.manta.vcf')

Then use ``to_bedpe_like`` method to generate bedpe-formatted pandas DataFrame.

.. ipython:: python
   :okexcept:

   vcf_bedpe_like = vcf.to_bedpe_like()
   print(vcf_bedpe_like)

If you want to add "SVLEN" and "CIPOS" fields, run as follows:

.. ipython:: python
   :okexcept:

   vcf_bedpe_like_with_info = vcf.to_bedpe_like(custom_infonames=['svlen', 'cipos'])
   print(vcf_bedpe_like_with_info)

To add FORMAT, set ``add_formats = True``:

.. ipython:: python
   :okexcept:

   vcf_bedpe_like_with_format = vcf.to_bedpe_like(add_formats=True)
   print(vcf_bedpe_like_with_format)

Do the same thing for adding FILTER:

.. ipython:: python
   :okexcept:

   vcf_bedpe_like_with_filter = vcf.to_bedpe_like(add_filters=True)
   print(vcf_bedpe_like_with_filter)