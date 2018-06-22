.. _tutorial:


Tutorial
=========

Before we go through the tutorial, I suppose ``TRIBEpipe`` was installed on your machine. If not, go to :ref:`installation` for more details.

``TRIBEpipe`` works like this with default parameters:

::

    TRIBEpipe -i tribe_rep1.fastq tribe_rep2.fastq \
              -gDNA gDNA.fastq \
              -wtRNA wtmRNA.fastq \
              -g dm6 \
              -o results \
              --cut 7 \
              --threads 8


Preparing genome index
-----------------------

RNA-seq reads mapped to reference genome using ``STAR``, user should provide the **STAR_index**, **Reference FASTA** and **GTF genome annotation** files.

You can download these files from the FTP site of ensembl_. and saved in your ``$HOME/data/genome``.

for example, download files for ``dm6``:

::

    # create folders for dm6
    $ mdkir -p $HOME/data/genome/dm6
    $ cd $HOME/data/genome/dm6
    $ mkdir annotation_and_repeats bigZips STAR_index

    # download files
    $ wget -O annotation_and_repeats/dm6.ensembl.gtf ftp://ftp.ensembl.org/pub/release-92/gtf/drosophila_melanogaster
    $ wget -O bigZips/dm6.fa ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/dna/

    # create STAR index
    $ STAR --runMode genomeGenerate --genomeDir STAR_index/ --genomeFastaFiles bigZips/dm6.fa --sjdbGTFfile annotation_and_repeats/dm6.ensembl.gtf --runThreadN 8

.. _ensembl: http://asia.ensembl.org/info/data/ftp/index.html


Preparing sequencing data
---------------------------

``TRIBEpipe`` requires **three** types of sequencing datasets:

- gDNA: genomic DNA sequencing 

- wtRNA: wild type RNA-seq data

- tribe: TRIBE RNA-seq data

Currently, only single-end (SE) reads supported.

If you have paired-end sequencing reads, you can choose the read represent the sense strand of RNA to ``TRIBEpipe``. For example, read2 of dUTP strand-specific RNA-seq.


Running demo data
------------------

Download the demo data from the following link: ftp://demo.tar.gz. 

:: 
    # demo data
    $ wget ftp://demo.tar.gz
    $ tar zxvf demo.tar.gz
    $ cd demo

Or, you can download the Hrp48 HyperTRIBE data from GEO with the following accession number:
(~ 7 GB to download in SRA format), see `Xu et al., RNA, 2018`_ for details.

.. _`Xu et al., RNA, 2018`: http://rnajournal.cshlp.org/content/24/2/173.long

::

    # HyperTRIBE data 
    SRR5944748 GSM2746730 Hrp48_HyperTRIBE_s2_rep1
    SRR5944749 GSM2746731 Hrp48_HyperTRIBE_s2_rep2
    SRR6426146 GSM2905815 s2_wt_mRNA
    SRR3177714 GSM2065948 s2_gDNA

The following tutorial using the **demo** data

:: 

    $ TRIBEpipe -i tribe_rep1.fastq tribe_rep2.fastq \
                -gDNA gDNA.fastq \
                -wtRNA wtmRNA.fastq \
                -g dm6 \
                -o results \
                --cut 7 \
                --threads 8

The RNA-seq libraries were constructed using NSR methods (REF [#]_), we will cut the fist 7 nucleutides, ``--cut=7``. and saving output in ``results`` directory.

.. [#] NSR reference

The directory structure of ``results`` should like this:

::

    results_20M/
    ├── gDNA
    │   ├── edits
    │   │   └── gDNA.edits.bedgraph
    │   ├── input_reads
    │   │   ├── gDNA.clean.fq
    │   │   ├── gDNA.clean.txt
    │   │   ├── gDNA.cutadapt.json
    │   │   └── gDNA.cutadapt.log
    │   └── mapping
    │       ├── demo
    │       ├── gDNA
    │       └── gDNA.mapping_stat.csv
    ├── TRIBE
    │   ├── edits
    │   │   ├── tribe_rep1.nodup.edits.bedgraph
    │   │   └── tribe_rep2.nodup.edits.bedgraph
    │   ├── edits_filted
    │   │   ├── tribe_rep1.nodup.edits.bedgraph
    │   │   └── tribe_rep2.nodup.edits.bedgraph
    │   ├── input_reads
    │   │   ├── tribe_rep1.clean.fq
    │   │   ├── tribe_rep1.clean.txt
    │   │   ├── tribe_rep1.cutadapt.json
    │   │   ├── tribe_rep1.cutadapt.log
    │   │   ├── tribe_rep2.clean.fq
    │   │   ├── tribe_rep2.clean.txt
    │   │   ├── tribe_rep2.cutadapt.json
    │   │   └── tribe_rep2.cutadapt.log
    │   └── mapping
    │       ├── demo
    │       ├── tribe_rep1
    │       ├── tribe_rep1.mapping_stat.csv
    │       ├── tribe_rep2
    │       └── tribe_rep2.mapping_stat.csv
    └── wt_RNA
        ├── edits
        │   └── wtRNA.edits.bedgraph
        ├── input_reads
        │   ├── wtRNA.clean.fq
        │   ├── wtRNA.clean.txt
        │   ├── wtRNA.cutadapt.json
        │   └── wtRNA.cutadapt.log
        └── mapping
            ├── demo
            ├── wtRNA
            └── wtRNA.mapping_stat.csv


There are three folders within ``resutls``: ``gDNA``, ``TRIBE`` and ``wt_RNA``.


within each folder, there are three sub-folders:

- ``input_reads`` : save the clean reads and \*.json statistics file 

- ``mapping`` : save the \*.bam files and \*.csv statistics file  

- ``edits``: save the \*.bedgraph file, not filtered editing events  

The finall results were saved in ``results/TRIBE/edits_filted`` in **BedGraph** format.



About results
---------------

Editing events were saved in **BedGraph** format in ``results/TRIBE/edits_filter/``

Before filter, editing events is a 12-column file in BedGraph format

::

    $ head results/TRIBE/edits/tribe_rep1.nodup.edits.bedgraph
    chr2L   73868   73869   15      chr2L_73869_32_15%      A       32      0       0       5       0       27
    chr2L   75790   75791   40      chr2L_75791_32_40%      A       32      17      0       13      0       2
    chr2L   103699  103700  12      chr2L_103700_24_12%     A       24      15      0       3       0       6
    chr2L   103718  103719  30      chr2L_103719_20_30%     A       20      8       0       6       0       6
    chr2L   103720  103721  10      chr2L_103721_20_10%     A       20      12      0       2       0       6
    chr2L   108661  108662  11      chr2L_108662_45_11%     A       45      37      0       5       0       3
    chr2L   108739  108740  12      chr2L_108740_48_12%     A       48      39      0       6       0       3
    chr2L   108740  108741  39      chr2L_108741_46_39%     A       46      25      0       18      0       3
    chr2L   108756  108757  10      chr2L_108757_47_10%     A       47      39      0       5       0       3
    chr2L   108766  108767  31      chr2L_108767_48_31%     A       48      30      0       15      0       3


    Column      Content
    1           chromosome
    2           chromStart
    3           chromEnd
    4           percentage of editing events
    5           name, chr_start_depth_pct%
    6           reference base
    7           read depth
    8           A count
    9           C count
    10          G count
    11          T count
    12          N count


After filter, editing events contains one extra column for the **gene** information


Running your data
-------------------

If your project contains mulitple TRIBE RNA-seq datasets while sharing the same group of **gDNA** and **wtRNA** samples, you could run the command like this:


**Option1.** Run all the TRIBE RNA-seq datasets in one command:

::
    
    $ TRIBEpipe -i tribe_A_rep1.fastq tribe_A_rep2.fastq tribe_B_rep1.fastq ... \
                -gDNA gDNA.fastq \
                -wtRNA wtmRNA.fastq \
                -g dm6 \
                -o results \
                --cut 7 \
                --threads 8


**Option2.** Run different TRIBE datasets separately

Make sure you are using the same folder each time, eg: ``results``

::

    $ TRIBEpipe -i tribe_A_rep1.fastq tribe_A_rep2.fastq \
                -gDNA gDNA.fastq \
                -wtRNA wtmRNA.fastq \
                -g dm6 \
                -o results \
                --cut 7 \
                --threads 8

    $ TRIBEpipe -i tribe_B_rep1.fastq tribe_B_rep2.fastq \
                -gDNA gDNA.fastq \
                -wtRNA wtmRNA.fastq \
                -g dm6 \
                -o results \
                --cut 7 \
                --threads 8

