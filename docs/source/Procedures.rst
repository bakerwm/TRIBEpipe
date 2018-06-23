.. _procedures:


Procedures
===========


.. note::

  filtering snps by:

  1. not within the first or last 6 bases of read





Here will explain the procedures of TRIBEpipe:

- Trimming, using cutadapt to trim adapter and low quality bases 
- Mapping, using STAR mappping clean reads to reference genome  
- Edits_parsing, extract editing events from BAM file
- Filtering, intersect with gDNA and wtRNA samples


Trimming
--------

Trim 3' adapter ``-a`` and low-quality bases ``-q 20``, cut 7-nt at 5' end (NSR library).

::

    $ cutadapt -a ATCTCGTATGCCGTCTTCTGCTTG \
               -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
               -a GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG \
               -a GCGGTTCAGCAGGAATGCCGAG \
               -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
               --cut=7 -m 19 -q 20 \
               --overlap=1 \
               --error-rate=0.1 \
               --times=4 \
               --cores=16 \
               --trim-n \
               --max-n=0.1 \
               tribe_rep1.fastq


Mapping
-------

Mapping reads to reference genome with STAR.

PCR duplicates were not removed for editing analysis.

::

    $ STAR --runMode alignReads  \
           --genomeDir star_index \
           --readFilesIn tribe_rep1.fastq \
           --readFilesCommand cat \
           --outFileNamePrefix out. \
           --runThreadN 16 \
           --limitOutSAMoneReadBytes 1000000 \
           --genomeLoad LoadAndKeep \
           --limitBAMsortRAM 10000000000 \
           --outSAMtype BAM SortedByCoordinate \
           --outReadsUnmapped Fastx \
           --outFilterMismatchNoverLmax 0.05 \
           --seedSearchStartLmax 20


Extract editing events
----------------------

Extract the nucleotide frequencies at each position on chromosome with ``samtools mpileup``

:: 

    # forward strand
    $ samtools view -f 16 -bhS in.bam | \
        samtools mpileup -f dm6.fa - > out_fwd.mpileup

    # reverse strand
    $ samtools view -F 16 -bhS in.bam | \
        samtools mpileup -f dm6.fa - > out_rev.mpileup

    # The output format of mpileup
    $ head out_fwd.mpileup
    chr2L   10598   A       1       ^~,     B
    chr2L   10599   A       1       ,       E
    chr2L   10600   T       1       ,       E
    chr2L   10601   C       1       ,       E
    chr2L   10602   T       1       ,       E
    chr2L   10603   T       1       ,       E
    chr2L   10604   T       1       ,       E
    chr2L   10605   G       1       ,       E
    chr2L   10606   A       1       ,       E
    chr2L   10607   A       1       ,       E

    column        content
    1            chromosome
    2            1-based coordinate
    3            reference base
    4            the number of reads covering the site
    5            read bases
    6            base qualities


Find more about ``mpileup`` at http://samtools.sourceforge.net/pileup.shtml and https://davetang.org/muse/2015/08/26/samtools-mpileup/ 

Converting mpileup to base frequency at each position

using the Python script: ``pileup2acgt`` from ``sequenza-utils`` (source_).

.. _source: https://bitbucket.org/sequenza_tools/sequenza-utils

The output looks like this:

::

    chr     n_base  ref_base        read.depth      A       C       G       T       strand
    chr2L   10598   A       1       1       0       0       0       0:0:0:0
    chr2L   10599   A       1       1       0       0       0       0:0:0:0
    chr2L   10600   T       1       0       0       0       1       0:0:0:0
    chr2L   10601   C       1       0       1       0       0       0:0:0:0

**TRIBE** 

- Forward strand, ref_base = A, G% >= 10%, read.depth >= 20  

- Reverse strand, ref_base = T, C% >= 10%, read.depth >= 20  

**gDNA**

- Forward strand, ref_base = A, A% >= 80%, G% == 0%

- Reverse strand, ref_base = T, T% >= 80%, C% == 0%

**wtRNA**

- Forward strand, ref_base = A, G% >= 10%, read.depth >= 10  

- Reverse strand, ref_base = T, C% >= 10%, read.depth >= 10  


Filtering
----------


Final results = (TRIBE intersect gDNA) exclude wtRNA


The criteria to define editing events:

- >= 20 reads in each replicate 

- in gDNA, A >= 80%, G = 0  

- A >= 10% in mRNA (editing)

Append the gene name to the editing record.


.. warning::

    HyperTRIBE is an improved version of TRIBE. The following are criteria to define editing events in TRIBE.

    Overall, A > 80% and G = 0 in gDNA, G > 0% in RNA

    + S2 cell

    20 reads and 10% editing

    + In neurons 

    A lower threshold (10 reads, 10% editing) was used to define endogenous editign events.

    All endogenous editing events detected were excluded from downstream analysis of TRIBE-expressing neurons.



``TRIBE`` RNA editing (A to I (G)) events were defined by the following rules:

In **TRIBE** samples:  

- read depth >= 20  
- editing percentage >= 10%  

In **gDNA** sample:  

- editing percentage = 0% (A to G)  
- A percentage >= 80%  

In **wtRNA** sample:  

- read depth >= 10  
- editing percentage >= 10%  

**TRIBE** sites = (TRIBE & gDNA) not wtRNA



