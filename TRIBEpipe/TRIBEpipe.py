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
import shlex, subprocess, pathlib
import pysam
import logging
from TRIBEpipe import trim
from TRIBEpipe import map
from TRIBEpipe import edits_parser
from TRIBEpipe import edits_filter


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
    parser.add_argument('-gDNA', required = False, metavar = 'gDNA', 
        type = argparse.FileType('r'),
        help = 'genomic DNA sample in fastq format, as control')
    parser.add_argument('-wtRNA', required = False, metavar = 'wtRNA', 
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
    parser.add_argument('--genome_data', metavar = 'genome_data', default = None,
        help = 'specify the directory contains genome data, eg: fasta, gtf, index \
                default: [$HOME/data/genome/]')
    args = parser.parse_args()
    return args



## trimming




def tribe_edits_parser(fqs, outdir, genome, genome_index, ad3, len_min, cut, 
                       threads, depth_cutoff, pct_cutoff):
    """extract editing events from TRIBE samples"""
    ## Trimming
    tribe_trim_dir = os.path.join(outdir, 'input_reads')
    tribe_clean_fq = trim.trim(fqs = fqs, adapter3 = ad3,
                               out_path = tribe_trim_dir, len_min = len_min, 
                               cut = cut, multi_cores = threads)

    ## Mapping
    tribe_map_dir = os.path.join(outdir, 'mapping')
    tribe_map_bam = []
    for fq in tribe_clean_fq:       
        b = map.map(fq, 'demo', tribe_map_dir, genome, threads, genome_index,
                    map_tools = 'STAR')
        bx = map.pcr_dup_remover(b[0])
        tribe_map_bam.append(bx)

    ## extract edits
    tribe_edits_dir = os.path.join(outdir, 'edits')
    tribe_edits = []
    for bam in tribe_map_bam:
        prefix = os.path.split(os.path.basename(bam))[0]
        bam_edits = os.path.join(tribe_edits_dir, prefix + '.edits.bedgraph')
        edits_parser.edits_parser(bam, genome, bam_edits, 'RNA', 
                                  depth_cutoff, pct_cutoff)
        tribe_edits.append(bam_edits)
    return tribe_edits


def gDNA_edits_parser(fq, outdir, genome, genome_index, ad3, len_min, cut, 
                       threads, depth_cutoff, pct_cutoff):
    """extract editing events from genomic DNA-seq sample"""
    ## Trimming
    gDNA_trim_dir = os.path.join(outdir, 'input_reads')
    gDNA_clean_fq = trim.trim(fqs = [fq], adapter3 = ad3,
                              out_path = gDNA_trim_dir, len_min = len_min, 
                              cut = cut, multi_cores = threads)
    
    ## Mapping
    gDNA_map_dir = os.path.join(outdir, 'mapping')
    gDNA_map_bam = map.map(gDNA_clean_fq[0], 'demo', gDNA_map_dir, genome,
                           threads, genome_index, map_tools = 'STAR')
    # bx = map.pcr_dup_remover(b[0])

    ## extract edits
    gDNA_edits_dir = os.path.join(outdir, 'edits')
    prefix = os.path.split(os.path.basename(gDNA_map_bam[0]))[0]
    bam_edits = os.path.join(gDNA_edits_dir, prefix + '.edits.bedgraph')
    gDNA_edits = edits_parser.edits_parser(gDNA_map_bam[0], genome, bam_edits, 
                                           'DNA', depth_cutoff, pct_cutoff)



def wtRNA_edits_parser(fq, outdir, genome, genome_index, ad3, len_min, cut, 
                       threads, depth_cutoff, pct_cutoff):
    """extract editing events from genomic DNA-seq sample"""
    ## Trimming
    wtRNA_trim_dir = os.path.join(outdir, 'input_reads')
    wtRNA_clean_fq = trim.trim(fqs = [fq], adapter3 = ad3,
                              out_path = wtRNA_trim_dir, len_min = len_min, 
                              cut = cut, multi_cores = threads)
    
    ## Mapping
    wtRNA_map_dir = os.path.join(outdir, 'mapping')
    wtRNA_map_bam = map.map(wtRNA_clean_fq[0], 'demo', wtRNA_map_dir, genome,
                           threads, genome_index, map_tools = 'STAR')
    # bx = map.pcr_dup_remover(b[0])

    ## extract edits
    wtRNA_edits_dir = os.path.join(outdir, 'edits')
    prefix = os.path.split(os.path.basename(wtRNA_map_bam[0]))[0]
    bam_edits = os.path.join(wtRNA_edits_dir, prefix + '.edits.bedgraph')
    wtRNA_edits = edits_parser.edits_parser(wtRNA_map_bam[0], genome, bam_edits, 
                                           'DNA', depth_cutoff, pct_cutoff)



def genome_parser(genome, path = None, aligner = 'STAR'):
    """parsing genome data"""
    if path is None:
        path = os.path.join(str(pathlib.Path.home()), 'data', 'genome')
    genome_fa = os.path.join(path, genome, 'bigZips', genome + '.fa')
    genome_gtf = os.path.join(path, genome, 'annotation_and_repeats', 
                              genome + '.ensembl.gtf')
    genome_index = os.path.join(path, genome, 'STAR_index')
    assert os.path.exists(genome_fa)
    assert os.path.exists(genome_gtf)
    assert os.path.exists(genome_index)
    return [genome_fa, genome_gtf, genome_index]



def main():
    args = get_args()


    ## prepare genome
    genome = args.g
    genome_index, genome_gtf, genome_fa = genome_parser(args.g, args.genome_data)

    ## prepare dir
    subdirs = [os.path.join(args.o, i) for i in ['TRIBE', 'gDNA', 'wt_RNA']]

    tribe_edits = tribe_edits_parser([i.name for i in args.i], subdirs[0], 
                                     genome, genome_index, args.a, args.m, 
                                     args.cut, args.threads, 
                                     args.tribe_depth_cutoff, 
                                     args.tribe_pct_cutoff)

    gDNA_edits = gDNA_edits_parser(args.gDNA.name, subdirs[1], 
                                   genome, genome_index, args.a, args.m, 
                                   args.cut, args.threads, 
                                   args.gDNA_depth_cutoff, 
                                   args.gDNA_pct_cutoff)

    wtRNA_edits = wtRNA_edits_parser(args.wtRNA.name, subdirs[2], 
                                     genome, genome_index, args.a, args.m, 
                                     args.cut, args.threads, 
                                     args.wtRNA_depth_cutoff, 
                                     args.wtRNA_pct_cutoff)

    ## filt
    for i in tribe_edits:
        i_dir = os.path.join(os.path.dirname(os.path.dirname(i)), 'edits_filted')
        i_name = os.path.basename(i)
        i_file = os.path.join(i_dir, i_name)
        edits_filter.edits_filter(i, [gDNA_edits, wtRNA_edits] , genome_gtf, i_file)

    print("This is the TRIBEpipe")


if __name__ == '__main__':
    main()


## EOF
