Installation
=============

Download this repository and enter the main folder, run ``python setup.py install --user`` to install ``TRIBEpipe``.

::

    git clone <repo>
    cd edits_parser
    python setup.py install --user


Requirements
-------------

- Python 3.6.4  
- fastqc >= 0.11.6  
- cutadapt >= 1.16  
- STAR >= 2.5.2a  
- samtools >= 1.4  
- bedtools >= 2.25.0  
- fastx_toolkit >= 0.0.14  
- numpy >= 1.14.0  
- pandas >= 0.22  
- pybedtools >= 0.7.10   
- pysam >= 0.11.1


Download genome data
---------------------

This pipeline require ``FASTA``, ``GTF``, ``STAR_index`` files located in ``$HOME/data/genome`` with the following structure.  

FASTA and GTF file could be downlaoded from ensembl ftp site: http://asia.ensembl.org/info/data/ftp/index.html 

For example, dm6

:: 

    dm6
    ├── annotation_and_repeats/dm6.ensembl.gtf
    ├── bigZips/dm6.fa
    └── STAR_index/


Test
-----

Installation is OK, if you can find the following message when typing ``TRIBEpipe -h`` in your console.

::

    $ TRIBEpipe -h
    usage: edits_parser [-h] -i TRIBE [TRIBE ...] [-gDNA gDNA] [-wtRNA wtRNA] -g
                        Genome -o OUTDIR [--tribe_depth_cutoff tribe_depth]
                        [--tribe_pct_cutoff percentage]
                        [--gDNA_depth_cutoff gDNA_depth]
                        [--gDNA_pct_cutoff percentage]
                        [--wtRNA_pct_cutoff percentage]
                        [--wtRNA_depth_cutoff wtRNA_depth] [-a adapter]
                        [-m min_length] [--cut cut N bases] [--threads Threads]
                        [--genome_data genome_data]

    Parsing editing events

    optional arguments:
      -h, --help            show this help message and exit
      -i TRIBE [TRIBE ...]  TRIBE sample in fastq format, (Single-end reads)
      -gDNA gDNA            genomic DNA sample in fastq format, as control
      -wtRNA wtRNA          wild-type RNA-seq data in fastq format, as control
      -g Genome             genome build of the reference, default: [dm6]
      -o OUTDIR             the directory to save results
      --tribe_depth_cutoff tribe_depth
                            minimum read depth at editing position for tribe
                            samples, default: 20
      --tribe_pct_cutoff percentage
                            minimum editing percentage [1-100]% for tribe samples,
                            default: 10
      --gDNA_depth_cutoff gDNA_depth
                            minimum read depth at editing position for genomic DNA
                            samples, default: 10
      --gDNA_pct_cutoff percentage
                            minimum percentage [1-100]% of reference base,
                            default: 80
      --wtRNA_pct_cutoff percentage
                            minimum editing percentage [1-100]% for wtRNA samples,
                            default: 10
      --wtRNA_depth_cutoff wtRNA_depth
                            minimum read depth at editing position for wtRNA
                            samples, default: 10
      -a adapter            adapter sequences at 3-prime end, default: TruSeq
      -m min_length         minimum length of reads, default: 19
      --cut cut N bases     cut N bases from reads, plus value at left of read,
                            minus value at right of read default: 0
      --threads Threads     number of threads to use for the pipeline, default: 1
      --genome_data genome_data
                            specify the directory contains genome data, eg: fasta,
                            gtf, index default: [$HOME/data/genome/]

    Output:
