.. _cli:

=======================
CLI Tutorial
=======================

Viola provides command line tools for several features.

Currently, following functions are available:

* VCF to BEDPE conversion
* Feature matrix generation for SV signature analysis
* SV signature extraction

----------------
Preparation
----------------

Several test VCF files are used in this tutorial.

Before going to the next step, download these files in your working directory.

.. code-block:: console

    $ curl -O 'https://raw.githubusercontent.com/dermasugita/ViolaDocs/main/docs/html/_static/tutorial.manta.vcf'
    $ mkdir -p signature_analysis/vcf
    $ curl https://raw.githubusercontent.com/dermasugita/ViolaDocs/main/docs/html/_static/signature_analysis/definitions.txt > signature_analysis/definitions.txt
    $ for i in $(seq 1 3); do curl 'https://raw.githubusercontent.com/dermasugita/ViolaDocs/main/docs/html/_static/signature_analysis/vcf/manta${i}.vcf' > signature_analysis/vcf/manta${i}.vcf; done

-------------------
Overview
-------------------

In your python environment (pyenv, virtualenv, pipenv, etc.), run ``viola``.

.. code-block:: console

    $ viola

    Usage: viola [OPTIONS] COMMAND [ARGS]...

    Options:
    --help  Show this message and exit.

    Commands:
    extract-signature
    generate-feature-matrix  Generate feature matrix from VCF or BEDPE files.
    vcf2bedpe                Convert a VCF file into a BEDPE file.

Actually, ``viola`` command has no option except ``--help``, so the de facto syntax is ``viola COMMAND [ARGS]``.

------------------
VCF to BEDPE
------------------

You can convert VCF files into a BEDPE files with ``viola vcf2bedpe``.

.. code-block:: console

    $ viola vcf2bedpe -h

    Usage: viola vcf2bedpe [OPTIONS] [VCF]

    Convert a VCF file into a BEDPE file.

    A VCF argument is the path to the input VCF file.

    Options:
    --version        Show the version and exit.
    --caller TEXT    The name of SV caller by which the input VCF was generated.
                    [manta, delly, lumpy, gridss] could be acceptable (default,
                    manta).

    -i, --info TEXT  The names of INFO fields to return. To specify multiple
                    INFO, separate them by commas. ex. --info SVTYPE,SVLEN,END

    -f, --filter     If specified, FILTER field of the VCF files is included in
                    output BEDPE.

    -m, --format     If specified, FORMAT field of the VCF files is included in
                    output BEDPE.

    -h, --help       Show this message and exit.


Now let's apply the example VCF file you got (See :ref:`Preparation<userguide/cli:Preparation>`) to the ``vcf2bedpe`` command.

.. code-block:: console

    $ viola vcf2bedpe --caller manta tutorial.manta.vcf

    chrom1     start1       end1 chrom2     start2       end2     name score strand1 strand2
    chr1   82550460   82550461   chr1   82554225   82554226    test1  None       +       -
    chr1   22814216   22814217   chr1   92581131   92581132    test2  None       -       -
    chr1   60567905   60567906   chr1   60675940   60675941    test3  None       +       -
    chr1   69583189   69583190   chr1   69590947   69590948    test4  None       +       -
    chr11  104534876  104534877  chr11  104536573  104536574    test5  None       +       -
    chr11  111134696  111134697  chr17   26470494   26470495  test6_1  None       +       -
    chr17   26470494   26470495  chr11  111134696  111134697  test6_2  None       -       +

The result will be output to the stdout by default.

You can add other VCF features, including FILTER, INFO, and FORMAT.

.. code-block:: console

    $ viola vcf2bedpe --caller manta --filter tutorial.manta.vcf

    chrom1     start1       end1 chrom2     start2       end2     name score strand1 strand2  MinSomaticScore   PASS
    chr1   82550460   82550461   chr1   82554225   82554226    test1  None       +       -             True  False
    chr1   22814216   22814217   chr1   92581131   92581132    test2  None       -       -             True  False
    chr1   60567905   60567906   chr1   60675940   60675941    test3  None       +       -             True  False
    chr1   69583189   69583190   chr1   69590947   69590948    test4  None       +       -            False   True
    chr11  104534876  104534877  chr11  104536573  104536574    test5  None       +       -            False   True
    chr11  111134696  111134697  chr17   26470494   26470495  test6_1  None       +       -             True  False
    chr17   26470494   26470495  chr11  111134696  111134697  test6_2  None       -       +             True  False

    $ viola vcf2bedpe --caller manta --info SVTYPE,SVLEN tutorial.manta.vcf

    chrom1     start1       end1 chrom2     start2       end2     name score strand1 strand2 svtype_0   svlen_0
    chr1   82550460   82550461   chr1   82554225   82554226    test1  None       +       -      DEL     -3764
    chr1   22814216   22814217   chr1   92581131   92581132    test2  None       -       -      INV  69766915
    chr1   60567905   60567906   chr1   60675940   60675941    test3  None       +       -      DEL   -108034
    chr1   69583189   69583190   chr1   69590947   69590948    test4  None       +       -      DEL     -7757
    chr11  104534876  104534877  chr11  104536573  104536574    test5  None       +       -      DEL     -1696
    chr11  111134696  111134697  chr17   26470494   26470495  test6_1  None       +       -      BND         0
    chr17   26470494   26470495  chr11  111134696  111134697  test6_2  None       -       +      BND         0

    $ viola vcf2bedpe --caller manta --format tutorial.manta.vcf

    chrom1     start1       end1 chrom2     start2       end2     name score strand1 strand2  sample1_N_PR_0  sample1_N_PR_1  sample1_N_SR_0  sample1_N_SR_1  sample1_T_PR_0  sample1_T_PR_1  sample1_T_SR_0  sample1_T_SR_1
    chr1   82550460   82550461   chr1   82554225   82554226    test1  None       +       -            21.0             0.0            10.0             0.0            43.0             4.0            15.0             3.0
    chr1   22814216   22814217   chr1   92581131   92581132    test2  None       -       -            24.0             0.0             NaN             NaN            35.0             5.0             NaN             NaN
    chr1   60567905   60567906   chr1   60675940   60675941    test3  None       +       -            23.0             0.0             NaN             NaN            44.0             6.0             NaN             NaN
    chr1   69583189   69583190   chr1   69590947   69590948    test4  None       +       -            21.0             0.0             NaN             NaN            20.0            12.0             NaN             NaN
    chr11  104534876  104534877  chr11  104536573  104536574    test5  None       +       -            22.0             0.0             NaN             NaN            57.0            14.0             NaN             NaN
    chr11  111134696  111134697  chr17   26470494   26470495  test6_1  None       +       -            12.0             0.0             NaN             NaN            45.0             5.0             NaN             NaN
    chr17   26470494   26470495  chr11  111134696  111134697  test6_2  None       -       +            12.0             0.0             NaN             NaN            45.0             5.0             NaN             NaN


-----------------------------------------
Feature Matrix Generation
-----------------------------------------

Simple feature matrix can be generated by ``viola generate-feature-matrix``.

.. code-block:: console

    $ viola generate-feature-matrix -h
    Usage: viola generate-feature-matrix [OPTIONS] OUTPUT

    Generate feature matrix from VCF or BEDPE files.

    Options:
    --version                       Show the version and exit.
    --input-dir TEXT                The directory of input files. When
                                    specified, the --files argument is disabled.

    --input-files TEXT              The input files separeted by comma. When
                                    specified, the --dir argument is disabled.

    --input-files-id TEXT           The sample ID of input files separeted by
                                    comma.

    --format [vcf|bedpe]            File format. vcf or bedpe.
    --caller [manta|delly|lumpy|gridss]
                                    The name of SV caller by which the input VCF
                                    was generated. This option can be specified
                                    when --format=vcf.

    --svtype-col-name TEXT          Name of the column of BEDPE files that
                                    indicate SV type. If not specified, SV type
                                    will be infered. This option can be
                                    specified when --format=bedpe

    --as-breakpoint                 Convert SVTYPE=BND records into breakpoint-
                                    wise SV records and infer its SVTYPE. This
                                    option is used when --format=vcf

    --definitions TEXT              Path to the definition file of custom SV
                                    class.

    -h, --help                      Show this message and exit.

To run this command, custom SV definition file is required. You've may already downloaded ``definitions.txt`` in :ref:`Preparation<userguide/cli:Preparation>`.

This is the content of ``definitions.txt``. Detailed description for the syntax is explained :ref:`here<signature_analysis:Definition File Syntax>`.

.. code-block::

    name "smaller DEL"
    0 SVTYPE == DEL
    1 SVLEN > -100
    logic 0 & 1

    name "larger DEL"
    0 SVTYPE == DEL
    logic 0

    name "smaller DUP"
    0 SVTYPE == DUP
    1 SVLEN < 100
    logic 0 & 1

    name "larger DUP"
    0 SVTYPE == DUP
    logic 0

    name "smaller INV"
    0 SVTYPE == INV
    1 SVLEN < 100
    logic 0 & 1

    name "larger INV"
    0 SVTYPE == INV
    logic 0

    name "translocation"
    0 SVTYPE == TRA
    logic 0

Example command for feature matrix generation.

.. code-block:: console

    $ viola generate-feature-matrix --input-dir signature_analysis/vcf --format vcf --caller manta --as-breakpoint --definitions signature_analysis/definitions.txt output.tsv

    ######### output.tsv ##########
    patients        smaller DEL     larger DEL      smaller DUP     larger DUP      smaller INV     larger INV      translocation   others
    manta1  2       3       1       0       2       1       2       0
    manta3  2       1       0       0       0       1       0       0
    manta2  2       2       1       2       2       1       2       0

Now you have a simple feature matrix.

