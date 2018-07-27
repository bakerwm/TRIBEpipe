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
3 filt TRIBE records by controls, gDNA, wtRNA

Compute each sample respectively:
TRIBE - gDNA - wt_mRNA = sites
"""

__author__ = "Ming Wang <wangm08@hotmail.com>"
__copyright__ = "2018 by Ming Wang <wangm08@hotmail.com>"
__license__ = "MIT"
__email__ = "wangm08@hotmail.com"
__version__ = "0.1"


import os, sys, argparse, tempfile, logging
import shlex, subprocess, pathlib
import pysam
import logging
from TRIBEpipe import trim
from TRIBEpipe import map
from TRIBEpipe import edits_parser
from TRIBEpipe import edits_filter
from TRIBEpipe.helper import *


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
    parser.add_argument('-gDNA', nargs = '+', required = False, 
        metavar = 'gDNA', type = argparse.FileType('r'),
        help = 'genomic DNA sample in fastq format, as control')
    parser.add_argument('-wtRNA', nargs = '+', required = False, 
        metavar = 'wtRNA', type = argparse.FileType('r'),
        help = 'wild-type RNA-seq data in fastq format, as control')
    parser.add_argument('-g', required = True, metavar = 'Genome',
        help = 'genome build of the reference, default: [dm6]')
    parser.add_argument('-o', required = True, metavar = 'OUTDIR',
        help = 'the directory to save results')
    parser.add_argument('--trimmed', action='store_true',
        help = 'if specified, do not run trimming, input fastq files are \
        supposed to be clean')
    parser.add_argument('--rm_dup', action='store_true',
        help = 'if specified, remove PCR duplicates')
    parser.add_argument('--tribe_depth_cutoff', default = 20, type = int, 
        metavar = 'tribe_depth',
        help = 'minimum read depth at editing position for tribe samples, \
        default: 20')
    parser.add_argument('--tribe_pct_cutoff', default = 10, type = int,
        metavar = 'percentage',
        help = 'minimum editing percentage [1-100]%% for tribe samples, \
        default: 10')
    parser.add_argument('--gDNA_depth_cutoff', default = 1, type = int, 
        metavar = 'gDNA_depth',
        help = 'minimum read depth at editing position for genomic DNA samples,\
        default: 1')
    parser.add_argument('--gDNA_pct_cutoff', default = 80, type = int,
        metavar = 'percentage',
        help = 'minimum percentage [1-100]%% of reference base, default: 80')
    parser.add_argument('--wtRNA_depth_cutoff', default = 10, type = int, 
        metavar = 'wtRNA_depth',
        help = 'minimum read depth at editing position for wtRNA samples, \
        default: 10')
    parser.add_argument('--wtRNA_pct_cutoff', default = 10, type = int,
        metavar = 'percentage',
        help = 'minimum editing percentage [1-100]%% for wtRNA samples, \
        default: 10')
    parser.add_argument('-a', metavar = 'adapter', default = None,
        help = 'adapter sequences at 3-prime end, default: TruSeq')
    parser.add_argument('-m', metavar = 'min_length', type = int, default = 19,
        help = 'minimum length of reads, default: 19')
    parser.add_argument('--cut', metavar = 'cut N bases', type = int, 
        default = 0,
        help = 'cut N bases from reads, plus value at left of read, minus value \
        at right of read default: 0')
    parser.add_argument('--threads', metavar = 'Threads', type = int, default = 1,
        help = 'number of threads to use for the pipeline, default: 1')
    parser.add_argument('--path_data', 
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args



def tribe_edits_parser(fqs, outdir, genome, ad3, len_min, cut, 
                       threads, depth_cutoff, pct_cutoff, merge = True, 
                       trimmed = False, rm_dup=True, path_data = None, 
                       overwrite = False):
    """extract editing events from TRIBE samples"""
    genome_fa = Genome_info(genome).get_fa()

    ## Trimming
    tribe_trim_dir = os.path.join(outdir, 'input_reads')
    if trimmed:
        logging.info('trimming skipped by --trimmed')
        tribe_clean_fq = fqs
    else:
        tribe_clean_fq = trim.trim(fqs = fqs, adapter3 = ad3,
                                   out_path = tribe_trim_dir, 
                                   len_min = len_min, 
                                   cut = cut, 
                                   multi_cores = threads,
                                   overwrite = False)
def map(fns, smp_name, path_out, genome, spikein=None, multi_cores=1, 
        aligner='bowtie2', path_data=None, overwrite=False):
    ## Mapping
    tribe_map_dir = os.path.join(outdir, 'mapping')
    # tribe_map_bam = []
    merge_name = str_common([os.path.basename(f) for f in tribe_clean_fq])
    merge_name = re.sub(".rep|_rep", "", merge_name)
    if merge is True:
        b = map.map(tribe_clean_fq, merge_name, tribe_map_dir, 
                    genome, multi_cores=threads, aligner = 'STAR', 
                    path_data = path_data, overwrite = overwrite)
        if rm_dup is True:
            tribe_map_bam = [map.pcr_dup_remover(i) for i in b]
        else:
            tribe_map_bam = b
    else:
        for fq in tribe_clean_fq:       
            b = map.map([fq], merge_name, tribe_map_dir, genome, 
                        multi_cores=threads, aligner = 'STAR', 
                        path_data = path_data, overwrite = overwrite)
            if rm_dup is True:
                bx = map.pcr_dup_remover(b[0])
            else:
                bx = b[0]
            tribe_map_bam.append(bx)

    ## extract edits
    tribe_edits_dir = os.path.join(outdir, 'edits')
    tribe_edits = []
    for bam in tribe_map_bam:
        prefix = os.path.splitext(os.path.basename(bam))[0]
        bam_edits = os.path.join(tribe_edits_dir, prefix + '.edits.bedgraph')
        edits_parser.edits_parser(bam, genome_fa, bam_edits, 'RNA', 
                                  depth_cutoff, pct_cutoff)
        tribe_edits.append(bam_edits)

    return tribe_edits



def gDNA_edits_parser(fqs, outdir, genome, ad3, len_min, cut, 
                       threads, depth_cutoff, pct_cutoff, merge = True,
                       trimmed = False, path_data = None, overwrite = False):
    """extract editing events from genomic DNA-seq sample"""
    genome_fa = Genome_info(genome).get_fa()
    ## Trimming
    gDNA_trim_dir = os.path.join(outdir, 'input_reads')
    if trimmed:
        logging.info('trimming skipped by --trimmed')
        gDNA_clean_fq = fqs
    else:
        gDNA_clean_fq = trim.trim(fqs = fqs, 
                                  adapter3 = ad3,
                                  out_path = gDNA_trim_dir, 
                                  len_min = len_min, 
                                  cut = cut, 
                                  multi_cores = threads,
                                  overwrite = overwrite)

    ## Mapping
    gDNA_map_dir = os.path.join(outdir, 'mapping')
    gDNA_map_bam = []
    merge_name = str_common([os.path.basename(f) for f in gDNA_clean_fq])
    merge_name = re.sub(".rep|_rep", "", merge_name)
    if merge is True:
        b = map.map(gDNA_clean_fq, merge_name, gDNA_map_dir, 
                    genome, multi_cores=threads, aligner = 'STAR', 
                    path_data = path_data, overwrite = overwrite)
        gDNA_map_bam.append(b[-1])
    else:
        for fq in gDNA_clean_fq:
            b = map.map([fq], merge_name, gDNA_map_dir, genome, 
                        multi_cores=threads, aligner = 'STAR', 
                        path_data = path_data, overwrite = overwrite)
            gDNA_map_bam.append(b)

    ## extract edits
    gDNA_edits_dir = os.path.join(outdir, 'edits')
    gDNA_edits = []
    for bam in gDNA_map_bam:
        prefix = os.path.splitext(os.path.basename(bam))[0]
        bam_edits = os.path.join(gDNA_edits_dir, prefix + '.edits.bedgraph')
        edits_parser.edits_parser(bam, genome_fa, bam_edits, 'DNA',
                                  depth_cutoff, pct_cutoff)
        gDNA_edits.append(bam_edits)

    return gDNA_edits



def wtRNA_edits_parser(fqs, outdir, genome, ad3, len_min, cut, 
                       threads, depth_cutoff, pct_cutoff, merge = True,
                       trimmed = False, rm_dup=True,
                       path_data = None, overwrite = False):
    """extract editing events from genomic DNA-seq sample"""
    genome_fa = Genome_info(genome).get_fa()
    ## Trimming
    wtRNA_trim_dir = os.path.join(outdir, 'input_reads')
    if trimmed:
        logging.info('trimming skipped by --trimmed')
        wtRNA_clean_fq = fqs
    else:
        wtRNA_clean_fq = trim.trim(fqs = fqs, 
                                   adapter3 = ad3,
                                   out_path = wtRNA_trim_dir, 
                                   len_min = len_min, 
                                   cut = cut, 
                                   multi_cores = threads,
                                   overwrite = overwrite)

    ## Mapping
    wtRNA_map_dir = os.path.join(outdir, 'mapping')
    # wtRNA_map_bam = []
    merge_name = str_common([os.path.basename(f) for f in wtRNA_clean_fq])
    merge_name = re.sub(".rep|_rep", "", merge_name)
    if merge is True:
        b = map.map(wtRNA_clean_fq, merge_name, wtRNA_map_dir, 
                    genome, multi_cores=threads, aligner = 'STAR', 
                    path_data = path_data, overwrite = overwrite)
        if rm_dup is True:
            tribe_map_bam = [map.pcr_dup_remover(i) for i in b]
        else:
            wtRNA_map_bam = b
    else:
        for fq in wtRNA_clean_fq:
            b = map.map([fq], merge_name, wtRNA_map_dir, genome, 
                        multi_cores=threads, aligner = 'STAR', 
                        path_data = path_data, overwrite = overwrite)
            if rm_dup is True:
                bx = map.pcr_dup_remover(b[0])
            else:
                bx = b[0]
            wtRNA_map_bam.append(bx)
    
    ## extract edits
    wtRNA_edits_dir = os.path.join(outdir, 'edits')
    wtRNA_edits = []
    for bam in wtRNA_map_bam:
        prefix = os.path.splitext(os.path.basename(bam))[0]
        bam_edits = os.path.join(wtRNA_edits_dir, prefix + '.edits.bedgraph')
        edits_parser.edits_parser(bam, genome_fa, bam_edits, 'RNA', 
                                  depth_cutoff, pct_cutoff)
        wtRNA_edits.append(bam_edits)

    return wtRNA_edits


def main():
    args = get_args()

    ## prepare genome
    genome = args.g
    genome_gtf = Genome_info(genome).ensembl_gtf()

    ## prepare dir
    subdirs = [os.path.join(args.o, i) for i in ['TRIBE', 'gDNA', 'wtRNA']]
    ## extract edits
    logging.info('step 1. processing TRIBE samples')
    if isinstance(args.i, str):
        tribe_edits = tribe_edits_parser([args.i.name, ], subdirs[0], 
                                         genome, args.a, 
                                         args.m, 
                                         args.cut, args.threads, 
                                         args.tribe_depth_cutoff, 
                                         args.tribe_pct_cutoff,
                                         trimmed = args.trimmed,
                                         rm_dup = args.rm_dup,
                                         path_data = args.path_data,
                                         overwrite = args.overwrite)
    elif isinstance(args.i, list):
        tribe_edits = tribe_edits_parser([i.name for i in args.i], subdirs[0], 
                                         genome, args.a, 
                                         args.m, 
                                         args.cut, args.threads, 
                                         args.tribe_depth_cutoff, 
                                         args.tribe_pct_cutoff,
                                         trimmed = args.trimmed,
                                         rm_dup = args.rm_dup,
                                         path_data = args.path_data,
                                         overwrite = args.overwrite)
    else:
        raise ValueError('illegal tribe samples, -i:')


    logging.info('step 2. processing genomic DNA sample')
    if isinstance(args.gDNA, str):
        gDNA_edits = gDNA_edits_parser(args.gDNA.name, subdirs[1], 
                                       genome, args.a, 
                                       args.m, 
                                       args.cut, args.threads, 
                                       args.gDNA_depth_cutoff, 
                                       args.gDNA_pct_cutoff,
                                       trimmed = args.trimmed,
                                       path_data = args.path_data,
                                       overwrite = args.overwrite)
    elif isinstance(args.gDNA, list):
        gDNA_edits = gDNA_edits_parser([i.name for i in args.gDNA], subdirs[1], 
                                       genome, args.a, 
                                       args.m, 
                                       args.cut, args.threads, 
                                       args.gDNA_depth_cutoff, 
                                       args.gDNA_pct_cutoff,
                                       trimmed = args.trimmed,
                                       path_data = args.path_data,
                                       overwrite = args.overwrite)
    else:
        gDNA_edits = None
        logging.error('gDNA skipped, -gDNA')


    logging.info('step 3. processing wildtype RNA-seq sample')
    if isinstance(args.gDNA, str):
        wtRNA_edits = wtRNA_edits_parser(args.wtRNA.name, subdirs[2], 
                                         genome, args.a, 
                                         args.m, 
                                         args.cut, args.threads, 
                                         args.wtRNA_depth_cutoff, 
                                         args.wtRNA_pct_cutoff,
                                         trimmed = args.trimmed,
                                         rm_dup = args.rm_dup,
                                         path_data = args.path_data,
                                         overwrite = args.overwrite)
    elif isinstance(args.wtRNA, list):
        wtRNA_edits = wtRNA_edits_parser([i.name for i in args.wtRNA], 
                                         subdirs[2], 
                                         genome, args.a, 
                                         args.m, 
                                         args.cut, args.threads, 
                                         args.wtRNA_depth_cutoff, 
                                         args.wtRNA_pct_cutoff,
                                         trimmed = args.trimmed,
                                         rm_dup = args.rm_dup,
                                         path_data = args.path_data,
                                         overwrite = args.overwrite)
    else:
        wtRNA_edits = None
        logging.error('wtRNA skipped, -wtRNA')

    logging.info('step 4. filtering TRIBE editing events')
    final_dir = os.path.join(args.o, 'TRIBE', 'edits_filted')
    for i in tribe_edits:
        i_file = os.path.join(final_dir, os.path.basename(i))
        edits_final = edits_filter.edits_filter(i, gDNA_edits, wtRNA_edits,
                                                genome_gtf, i_file, False)

    logging.info('finish.')


if __name__ == '__main__':
    main()


## EOF
