************************
Welcom to Viola Package!
************************

| Current version: 1.0.0.dev5
| Official Document --> `here`_.

.. _here: https://dermasugita.github.io/ViolaDocs/docs/html/index.html

.. image:: https://dermasugita.github.io/ViolaDocs/docs/html/_static/Viola-logo/JPG/wide.jpg


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

For command line interfaces, see CLI tutorial.

If you want to get a deeper understanding of how Viola objects are structured, how filtering works and how merging works, please refer to the `User Guide`_.

Manuals for individual classes/methods/functions are available in the `API Reference`_.

Documentation
=============

- `Quick Start`_
- `Signature Analysis Tutorial`_
- `API Reference`_

.. _Quick Start: https://dermasugita.github.io/ViolaDocs/docs/html/quickstart.html
.. _Signature Analysis Tutorial: https://dermasugita.github.io/ViolaDocs/docs/html/signature_analysis.html
.. _API Reference: https://dermasugita.github.io/ViolaDocs/docs/html/reference/index.html
.. _User Guide: https://dermasugita.github.io/ViolaDocs/docs/html/userguide/index.html