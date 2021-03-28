.. _data_structure:

======================
Viola's Data Structure
======================

| Viola employs relational-database-like structure to enable flexible data proccessing.
| Here, we explain the specific data structure of Bedpe and Vcf class.

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