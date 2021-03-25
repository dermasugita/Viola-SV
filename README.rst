************************
Welcom to Viola Package!
************************

.. image:: https://dermasugita.github.io/PySgtDocs/docs/html/_static/Viola-logo/JPG/wide.jpg


Viola is a flexible and powerful python package designed specifically for analysis of genomic structural variant (SV) signatures.
We provide following tools for SV signature analysis:

* Custom SV classification tool
* Feature matrix generator 
* SV signature extractor (NMF) with stability evaluation system.
* Signature attributer for given sample

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

To import PySgt in your script, simply run below:

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

As a first step we recommend the :ref:`Quick Start<quickstart>` page where you can see the basic behaviour of Viola.
A good second step is to read the :ref:`Signature Analysis Tutorial<signature_analysis>` page.

For command line interfaces, see CLI tutorial.

If you want to get a deeper understanding of how Viola objects are structured, how filtering works and how merging works, please refer to the :ref:`User Guide<user_guide>`.

Manuals for individual classes/methods/functions are available in the :ref:`API Reference<api_reference>`.

Documentation
=============

- :ref:`Quick Start<quickstart>`
- :ref:`API Reference<api_reference>`

