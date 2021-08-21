.. _data_structure:

==============================================
Viola's Data Structure – Why tidy data?
==============================================

Viola breaks down the SV information into multiple tidy tables to enable flexible data proccessing.
These tables follow the principles of tidy data, i.e., each SV record is a row, each variable is a column, and each type of observational unit is a table (See Also: `Wickham, 2014`_). 
Consequently, storage of multiple values in one element is avoided, in contrast to the INFO and FORMAT columns of a VCF file.

The most important benefit of this data structure is the extendibility for future functions.
Software always requires updates to keep up with trends and developments. By breaking down the SV information into multiple tidy tables, even complex requirements can be implemented with little effort.

These points are rather the developer side benefits than the user side ones. However, maintaining data structure that is easy to develop is an important part of keeping software healthy, and is ultimately passed on to the user. 

Here, we explain the specific data structure of Bedpe and Vcf class.

.. _`Wickham, 2014`: https://www.jstatsoft.org/article/view/v059i10

-------------------------------
Tables composing Vcf object
-------------------------------

| In the following description tables, the primary keys are in bold.
| If the table has a foreign key (FK), the dependencies are noted below the table.


~~~~~~~~~~~~~~~
positions table
~~~~~~~~~~~~~~~

Name of the table: "positions"

==========  =========================================
id          Identifiers of SV records                
==========  =========================================
chrom1(FK)  Chromosome of the first breakend         
pos1        Position of the first breakend           
chrom2(FK)  Chromosome of the second breakend        
pos2        Position of the second breakend          
strand1     Strand of the first breakend             
strand2     Strand of the second breakend            
ref         Reference nucleotides in the chrom1:pos1 
alt(FK)     Details about the alteration             
svtype      SV type of the event                     
==========  =========================================

| Forein keys:
|    chrom1: "contigs_meta" table
|    chrom2: "contigs_meta" table


~~~~~~~~~~~~~~
filters table
~~~~~~~~~~~~~~

Name of the table: "filters"

======= ==========================
id(FK)  Identifiers of SV records
------- --------------------------
filter  Filters assigned to the id
======= ==========================
======= ==========================

| Forein keys:
|    id: "positions" table


~~~~~~~~~~~~~~
info tables
~~~~~~~~~~~~~~

Vcf object contains the same number of info tables as the number of INFO (e.g. SVLEN, CIPOS, IMPRECISE, ...) in the original vcf file.
Here, we show the 'svlen' table as an example of info tables.

Name of the table: "svlen"

========= =============================================================================================
id(FK)    Identifiers of SV records
--------- ---------------------------------------------------------------------------------------------
value_idx 0 origin index (this column is necessary since some INFOs have multiple values such as CIPOS)
========= =============================================================================================
svlen     Length of SV
========= =============================================================================================

| Forein keys:
|    id: "positions" table

~~~~~~~~~~~~~~
formats table
~~~~~~~~~~~~~~

Name of the table: "formats"

============ ==========================
id(FK)       Identifiers of SV records
------------ --------------------------
sample(FK)   Sample name
------------ --------------------------
format(FK)   Name of the FORMAT
------------ --------------------------
value_idx    0 origin index 
============ ==========================
value        Value of the FORMAT
============ ==========================

| Foreign keys:
|    id: "positions" table
|    sample: "samples_meta" table
|    formats: "formats_meta" table

~~~~~~~~~~~~~~~~~~~
contigs_meta table
~~~~~~~~~~~~~~~~~~~

This table corresponds to the CONTIG section in header of vcf files.

Name of the table: "contigs_meta"

========= ==============================================
id        Name of a contig (in most case, chromosome)
========= ==============================================
length    Length of the contig
========= ==============================================

| Note:
| In this table, "id" doesn't mean the "SV ID" but means the "identifier of the contigs".

~~~~~~~~~~~~~~~~~~~
alts_meta table
~~~~~~~~~~~~~~~~~~~

Name of the table: "alts_meta"

=========== ==============================================
id          Name of an ALT 
=========== ==============================================
description Description of the ALT
=========== ==============================================

~~~~~~~~~~~~~~~~~~~
filters_meta table
~~~~~~~~~~~~~~~~~~~

Name of the table: "filters_meta"

=========== ==============================================
id          Name of a FILTER 
=========== ==============================================
description Description of the FILTER
=========== ==============================================

~~~~~~~~~~~~~~~~~~~
infos_meta table
~~~~~~~~~~~~~~~~~~~

Name of the table: "infos_meta"

=========== ==============================================
id          Name of a INFO
=========== ==============================================
number      Number of values included with the INFO
type        Type of the value of the INFO
description Description of the INFO
source      Source of the annotation (e.g. "dbsnp")
version     Version of the source (e.g. "138")
=========== ==============================================

~~~~~~~~~~~~~~~~~~~
formats_meta table
~~~~~~~~~~~~~~~~~~~

Name of the table: "formats_meta"

=========== ==============================================
id          Name of a FORMAT
=========== ==============================================
number      Number of values included with the FORMAT
type        Type of the value of the FORMTA
description Description of the FORMAT
source      Source of the annotation (e.g. "dbsnp")
version     Version of the source (e.g. "138")
=========== ==============================================

~~~~~~~~~~~~~~~~~~~
samples_meta table
~~~~~~~~~~~~~~~~~~~

This table has only a single column.

Name of the table: "samples_meta"

=========== ==============================================
id          Name of a sample
=========== ==============================================
=========== ==============================================