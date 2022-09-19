************************
Welcom to Viola Package!
************************

.. _here: https://viola-sv.readthedocs.io

.. image:: https://viola-sv.readthedocs.io/en/latest/_static/Viola-logo/JPG/wide.jpg
   

|PyPI|_ |PyVersion|_ |License|_
|DOI|_ |DOC|_ |COV|_
|Issue|_ |Downloads|_

.. |PyPI| image:: https://img.shields.io/pypi/v/viola-sv
.. _PyPI: https://pypi.org/project/Viola-SV/

.. |PyVersion| image:: https://img.shields.io/pypi/pyversions/viola-sv
.. _PyVersion: https://pypi.org/project/Viola-SV/

.. |License| image:: https://img.shields.io/pypi/l/viola-sv
.. _License: https://pypi.org/project/Viola-SV/

.. |DOI| image:: https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtab662-9cf
.. _DOI: https://academic.oup.com/bioinformatics/article/38/2/540/6371863

.. |DOC| image:: https://readthedocs.org/projects/viola-sv/badge/?version=latest
.. _DOC: https://viola-sv.readthedocs.io

.. |COV| image:: https://codecov.io/gh/dermasugita/Viola-SV/branch/master/graph/badge.svg?token=G7TI1S6FY2 
.. _COV: https://codecov.io/gh/dermasugita/Viola-SV

.. |Issue| image:: https://img.shields.io/badge/issue%20tracking-github-brightgreen
.. _Issue: https://github.com/dermasugita/Viola-SV/issues

.. |Downloads| image:: https://pepy.tech/badge/viola-sv/month
.. _Downloads: https://pepy.tech/project/viola-sv

Overview
==============

Viola is a flexible and powerful python package designed specifically for analysis of genomic structural variant (SV) signatures.
We provide following tools for SV signature analysis:

* Custom SV classification tool
* Feature matrix generator 
* SV signature extractor (NMF) with stability evaluation system.

In addition to these, Viola offers a number of other useful utilities, including:

* VCF filter that accepts genomic coordinates and INFO/FORMAT columns.
* VCF merging tool with user-defined thresholds.
* Breakends-to-breakpoint converter with SVTYPE inference.
* Microhomology inference.
* VCF-to-BEDPE conversion.
* Command line tools for light users.

Currently, Viola supports four popular SV callers as input VCF files:

* Manta
* Delly
* Lumpy
* Gridss

In the future, we plan to support more SV callers!  
  
Installation
=========================

The package can be installed with ``pip``:

.. code:: console

   $ pip install viola-sv

To import Viola in your script, simply run below:

.. code:: python
   
   import viola

Command line tools:

.. code:: console

   $ viola <command> [something]


Prerequisites
==============

Python version 3.6 or newer.

Recommended Environment
=======================

* OS
   * Linux
   * Mac
* Script Manager
   * Jupyter Notebook/Jupyter Lab
   * Google Colaboratory

Docker support is coming soon!

How to Learn Viola
===================

As a first step we recommend the `Quick Start`_ page where you can see the basic behaviour of Viola.
A good second step is to read the `Signature Analysis Tutorial`_ page.

For command line interfaces, see `CLI tutorial`_.

If you want to get a deeper understanding of how Viola objects are structured, how filtering works and how merging works, please refer to the `User Guide`_.

Manuals for individual classes/methods/functions are available in the `API Reference`_.

For further discussion, please join our `google group`_! We're waiting for your comments, ideas, and requests.
The questions about how to use Viola-SV are also welcome!  

.. _google group: https://groups.google.com/g/viola-users

Documentation
=============

- `Documentation Home`_
- `Quick Start`_
- `Signature Analysis Tutorial`_
- `API Reference`_

.. _Documentation Home: https://viola-sv.readthedocs.io/
.. _Quick Start: https://viola-sv.readthedocs.io/en/1.0.x-stable-doc/quickstart.html
.. _Signature Analysis Tutorial: https://viola-sv.readthedocs.io/en/1.0.x-stable-doc/signature_analysis.html
.. _API Reference: https://viola-sv.readthedocs.io/en/1.0.x-stable-doc/reference/index.html
.. _User Guide: https://viola-sv.readthedocs.io/en/1.0.x-stable-doc/userguide/index.html
.. _CLI tutorial: https://viola-sv.readthedocs.io/en/1.0.x-stable-doc/userguide/cli.html

Publications
=============

Sugita, I., Matsuyama, S., Dobashi, H., Komura, D. & Ishikawa, S. Viola: a structural variant signature extractor with user-defined classifications. Bioinformatics (2021) doi:10.1093/bioinformatics/btab662.
