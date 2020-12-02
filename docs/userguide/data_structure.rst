.. _data_structure:

======================
PySgt's Data Structure
======================

PySgt employs relational database like structure to enable flexible data proccessing.

-------------------------------
Tables composing Vcf object
-------------------------------

~~~~~~~~~~~~~~~
positions table
~~~~~~~~~~~~~~~

+---------+------------------------------------------+
| id      | Identifiers of SV records                |
+=========+==========================================+
| chrom1  | Chromosome of the first breakend         |
+---------+------------------------------------------+
| pos1    | Position of the first breakend           |
+---------+------------------------------------------+
| chrom2  | Chromosome of the second breakend        |
+---------+------------------------------------------+
| pos2    | Position of the second breakend          |
+---------+------------------------------------------+
| strand1 | Strand of the first breakend             |
+---------+------------------------------------------+
| strand2 | Strand of the second breakend            |
+---------+------------------------------------------+
| ref     | Reference nucleotides in the chrom1:pos1 |
+---------+------------------------------------------+
| alt     | Details about the alteration             |
+---------+------------------------------------------+
| svtype  | SV type of the event                     |
+---------+------------------------------------------+


~~~~~~~~~~~~~~
filters table
~~~~~~~~~~~~~~

======= ==========================
id      Identifiers of SV records
======= ==========================
filters Filters assigned to the id
======= ==========================