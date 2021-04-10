.. _merge:

==============
VCF Merging
==============

Here, we describe how to merge VCF files with Viola.

First of all, let's load the four VCF files from the different SV callers.

.. ipython:: python

    manta = viola.read_vcf('https://raw.githubusercontent.com/dermasugita/Viola-SV/master/examples/demo_merge/test.merge.manta.vcf', variant_caller='manta')
    delly = viola.read_vcf('https://raw.githubusercontent.com/dermasugita/Viola-SV/master/examples/demo_merge/test.merge.delly.vcf', variant_caller='delly')
    lumpy = viola.read_vcf('https://raw.githubusercontent.com/dermasugita/Viola-SV/master/examples/demo_merge/test.merge.lumpy.vcf', variant_caller='lumpy')
    gridss = viola.read_vcf('https://raw.githubusercontent.com/dermasugita/Viola-SV/master/examples/demo_merge/test.merge.gridss.vcf', variant_caller='gridss')

These are synthetic files that we have created to briefly test the functionality of :doc:`viola.merge<../reference/api/viola.merge>`, and they share some SVs.

You can see which SV records in a file are sharing SVs with other files by looking at their SV ID.

.. ipython:: python

    print(manta)

In the above case, for example, SV ID "ML1" means "the SV record is intended to be shared by Manta and Lumpy".

.. ipython:: python

    print(delly)
    print(lumpy)
    print(gridss)

-------------------
Default Behaviour
-------------------

Now, let's use :doc:`viola.merge<../reference/api/viola.merge>` function.

.. ipython:: python

    merged_vcf = viola.merge([manta, delly, lumpy, gridss], threshold=100)
    print(merged_vcf)

You can see the SV records that are within 100bp of each other have been successfully merged.

By default, for each merged event, only one SV record is kept; regarding the "ML1", the SV record from Manta is selected and that from Lumpy is discarded.

You can control which SV callers are given priority for selection by changing the order of the list of input Vcf objects.

Here is the example.

.. ipython:: python

    merged_vcf = viola.merge([delly, lumpy, manta, gridss], threshold=100)
    print(merged_vcf)

Look again at "ML1" and you will see that the Lumpy-derived SV record has been chosen. In this case, Lumpy is the 2nd priority and Manta is the 3rd.

------------------------------
INFO relating to VCF Merging
------------------------------

:doc:`viola.merge<../reference/api/viola.merge>` adds three INFOs:

* mergedid: ID given to merged SV records.
* originalid: Original ID of each SV record; the name of the caller is prepended after the merge to avoid conflict between the names of the SV IDs of different callers. 
* caller: The name of the caller which identified the SV record.
* supportingid: IDs of SV records supporting the merged SV record.
* supportingcaller: SV callers supporting the variant.
* supportingidcount: Number of SV records supporting the merged SV record.
* supportingcallercount: Count of SV callers supporting the variant.

As examples, 'mergedid' and 'supportingcallercount' are shown below.

.. ipython:: python

    print(merged_vcf.view(custom_infonames=['mergedid', 'supportingcallercount']))

Now, in order to increase the accuracy of SV detection, let's get SVs detected by more than one SV caller.
To do this, just select the ones with 'supportingcallercount' greater than or equal to 2.

.. ipython:: python

    filtered_vcf = merged_vcf.filter('supportingcallercount >= 2')
    print(filtered_vcf.view(custom_infonames=['mergedid', 'supportingcallercount']))

-----------------------
Export VCF file
-----------------------

Finally, let's export the Vcf object filtered by 'supportingcallercount' as a VCF file.

.. code-block:: python

    filtered_vcf.to_vcf('/path/to/the/output.vcf')