#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TRIBEpipe 
parsing TRIBE data and control data


TRIBE sample:
1. count > N (20)
2. freqG > 10%

genomic DNA:
A% > 80%, G% = 0

background mRNA:
count > N (10)
freqG > 10%


forward strand: A -> G conversion
reverse strand: T -> C conversion

pipeline

1 fastqc and mapping reads to reference genome
2 extract edits positions
2.1 samtools view: separate forward strand (-f 16) and reverse strand (-F 16)
2.2 samtools mpileup: filt and extract nucleotide content at each position
2.3 convert mpileup to BED format
2.4 filter positions
3 filt TRIBE records by controls, gDNA, wt_RNA

Compute each sample respectively:
TRIBE - gDNA - wt_mRNA = sites
"""

import os, sys, argparse, tempfile, logging
import shlex, subprocess
import pysam
import logging
from TRIBEpipe import trim
from TRIBEpipe import map
from TRIBEpipe import edits_filter
from TRIBEpipe import edits_parser


logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'edits_parser', description = 'Parsing editing events',
        epilog = 'Output: ')
    parser.add_argument('-i', nargs = '+', required = True, metavar = 'TRIBE', 
        type = argparse.FileType('r'),
        help = 'TRIBE sample in fastq format, (Single-end reads)')
    parser.add_argument('-gDNA', required = True, metavar = 'gDNA', 
        type = argparse.FileType('r'),
        help = 'genomic DNA sample in fastq format, as control')
    parser.add_argument('-wtRNA', required = True, metavar = 'wtRNA', 
        type = argparse.FileType('r'),
        help = 'wild-type RNA-seq data in fastq format, as control')
    parser.add_argument('-g', required = True, metavar = 'Genome',
        type = argparse.FileType('r'),
        help = 'genome build of the reference, default: [dm6]')
    parser.add_argument('-o', required = True, metavar = 'OUTDIR',
        help = 'the directory to save results')
    parser.add_argument('--tribe_depth_cutoff', default = 20, type = int, 
        metavar = 'tribe_depth',
        help = 'minimum read depth at editing position for tribe samples, default: 20')
    parser.add_argument('--tribe_pct_cutoff', default = 10, type = int,
        metavar = 'percentage',
        help = 'minimum editing percentage [1-100]%% for tribe samples, default: 10')
    parser.add_argument('--gDNA_depth_cutoff', default = 10, type = int, 
        metavar = 'gDNA_depth',
        help = 'minimum read depth at editing position for genomic DNA samples, default: 10')
    parser.add_argument('--gDNA_pct_cutoff', default = 80, type = int,
        metavar = 'percentage',
        help = 'minimum percentage [1-100]%% of reference base, default: 80')
    parser.add_argument('--wtRNA_pct_cutoff', default = 10, type = int,
        metavar = 'percentage',
        help = 'minimum editing percentage [1-100]%% for wtRNA samples, default: 10')
    parser.add_argument('--wtRNA_depth_cutoff', default = 10, type = int, 
        metavar = 'wtRNA_depth',
        help = 'minimum read depth at editing position for wtRNA samples, default: 10')
    parser.add_argument('-a', metavar = 'adapter', 
        help = 'adapter sequences at 3-prime end, default: TruSeq')
    parser.add_argument('-m', metavar = 'min_length', type = int, default = 19,
        help = 'minimum length of reads, default: 19')
    parser.add_argument('--cut', metavar = 'cut N bases', type = int, default = 0,
        help = 'cut N bases from reads, plus value at left of read, minus value at right of read\
                default: 0')
    parser.add_argument('--threads', metavar = 'Threads', type = int, default = 1,
        help = 'number of threads to use for the pipeline, default: 1')
    parser.add_argument('--genome_data', metavar = 'genome_data', 
        help = 'specify the directory contains genome data, eg: fasta, gtf, index \
                default: [$HOME/data/genome/]')
    args = parser.parse_args()
    return args



## trimming





















def main():
    args = get_args()
    print("This is the TRIBEpipe")



if __name__ == '__main__':
    main()


## EOF
