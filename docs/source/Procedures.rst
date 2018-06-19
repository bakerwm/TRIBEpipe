Procedures
===========


Here will explain the procedures of TRIBEpipe:

- Trimming, using cutadapt to trim adapter and low quality bases 
- Mapping, using STAR mappping clean reads to reference genome  
- Edits_parsing, extract editing events from BAM file
- Filtering, intersect with gDNA and wtRNA samples


Trimming
----------

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
---------

Mapping reads to reference genome with STAR

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
-----------------------

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

	column		content
	1			chromosome
	2			1-based coordinate
	3			reference base
	4			the number of reads covering the site
	5			read bases
	6			base qualities


Find more about ``mpileup`` at http://samtools.sourceforge.net/pileup.shtml and https://davetang.org/muse/2015/08/26/samtools-mpileup/ 
