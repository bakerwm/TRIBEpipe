#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mapping reads to reference genome
1. to spike-in genome
2. to repbase
3. to rRNA, tRNA, MT
4. to genome

Note:
1. filtering unique mapped reads
STAR --outFilterMultimapNmax 1 

"""

__author__ = "Ming Wang <wangm08@hotmail.com>"
__copyright__ = "2018 by Ming Wang <wangm08@hotmail.com>"
__license__ = "MIT"
__email__ = "wangm08@hotmail.com"
__version__ = "0.1"


import os
import sys
import re
import json
import glob
import argparse
import shlex
import subprocess
import pandas as pd
import binascii
import logging
import gzip
import pysam
import pybedtools
from TRIBEpipe.helper import *

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'map',
        description = 'Mapping RNA-seq reads to genome',
        epilog = 'Example: map -i f1.fq f2.fq -n demo -o output')
    parser.add_argument('-i', nargs = '+', required = True, metavar = 'INPUT', 
        type = argparse.FileType('r'),
        help = 'Sequencing reads in FASTQ format, 1-4 files.')
    parser.add_argument('-n', required = True, metavar = 'NAME',
        help = 'Name of the experiment')
    parser.add_argument('-o', required = True, default = None, 
        metavar = 'OUTPUT', help = 'The directory to save results.')
    parser.add_argument('-g', default = 'hg19', metavar = 'GENOME',
        help = 'Reference genome : dm3, dm6, hg19, GRCh38, mm10, GRCm38')
    parser.add_argument('-x', required = True, metavar = 'index',
        help = 'Index of the genome for alignment tools')
    parser.add_argument('-t', required = True, default = 'bowtie', 
        metavar = 'Aligner',
        help = 'Aligner for the mapping, bowtie, bowtie2, STAR')
    parser.add_argument('-p', default = 1, metavar = 'Threads', type = int, 
        help = 'Number of threads to launch, default [1].')
    parser.add_argument('--rmdup', action = 'store_true',
        help = 'remove PCR duplicates using Picard, if specified')
    parser.add_argument('--path_data', 
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args



def bowtie2_se(fn, idx, path_out, para=1, multi_cores=1, overwrite=False):
    """
    Mapping SE reads to idx using Bowtie
    filter uniquely mapped reads by samtools:
    samtools view -bhS -q 10 in.bam > out.bam
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert is_idx(idx, 'bowtie2')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    assert isinstance(para, int)
    ## parameters
    para_v = {1: '--sensitive', 2: '--local'}
    para_bowtie2 = para_v[para] if para in para_v else ''
    fn_type = seq_type(fn)
    if fn_type == 'fasta':
        para_bowtie2 += ' -f'
    elif fn_type == 'fastq':
        para_bowtie2 += ' -q'
    else:
        raise ValueError('unknown type of reads')
    ## prefix
    fn_prefix = file_prefix(fn)[0]
    fn_prefix = re.sub('\.clean|\.nodup|\.cut', '', fn_prefix)
    # fn_prefix = re.sub('_[12]|_R[12]$', '', fn_prefix)
    idx_name = os.path.basename(idx)
    fn_unmap_file = os.path.join(path_out, '%s.not_%s.%s' % (fn_prefix, idx_name, fn_type))
    fn_map_prefix = os.path.join(path_out, fn_prefix)
    fn_map_bam = fn_map_prefix + '.map_%s.bam' % idx_name
    fn_map_bed = fn_map_prefix + '.map_%s.bed' % idx_name
    fn_map_log = fn_map_prefix + '.map_%s.bowtie2.log' % idx_name
    if os.path.exists(fn_map_bam) and os.path.exists(fn_unmap_file) and overwrite is False:
        logging.info('file exists: %s' % fn_map_bam)
    else:
        c1 = 'bowtie2 %s -p %s --mm --no-unal --un %s -x %s %s' % (para_bowtie2,
              multi_cores, fn_unmap_file, idx, fn)
        c2 = 'samtools view -bhS -q 10 -F 0x4 -@ %s -' % multi_cores
        c3 = 'samtools sort -@ %s -o %s -' % (multi_cores, fn_map_bam)
        with open(fn_map_log, 'wt') as fo:
            p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE, stderr=fo)
            p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
            p4 = p3.communicate()
        pysam.index(fn_map_bam)
        pybedtools.BedTool(fn_map_bam).bam_to_bed().saveas(fn_map_bed)
    ## statistics
    d = bowtie2_log_parser(fn_map_log)

    return [fn_map_bam, fn_unmap_file]



def star_se(fn, idx, path_out, para=1, multi_cores=1, overwrite=False):
    """
    mapping single read to one index using STAR
    Input: fastq|a
    Output: bam (sorted), unmapped reads
    #
    filtering unique mapped reads by samtools
    #
    STAR --runMode alignReads \
         --genomeDir /path/to/genome \
         --readFilesIn /path/to/reads \
         --readFilesCommand cat \
         --outFileNamePrefix /name \
         --runThreadN 8 \
         --limitOutSAMoneReadBytes 1000000 \
         --genomeLoad LoadAndRemove \
         --limitBAMsortRAM 10000000000 \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMismatchNoverLMax 0.05 \
         --seedSearchStartLmax 20


    STAR  --runThreadN 8 \
          --outFilterMismatchNoverLmax 0.07 \
          --outFileNamePrefix $prefix"_" \
          --outFilterMatchNmin 16 \
          --outFilterMultimapNmax 1 \
          --genomeDir $star_indices \
          --readFilesIn $input
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert is_idx(idx, 'star')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    fn_type = seq_type(fn)
    freader = 'zcat' if is_gz(fn) else '-'
    ## prefix
    fn_prefix = file_prefix(fn)[0]
    fn_prefix = re.sub('\.clean|\.nodup|\.cut', '', fn_prefix)
    # fn_prefix = re.sub('_[12]|_R[12]$', '', fn_prefix)
    idx_name = os.path.basename(idx)
    fn_unmap_file = os.path.join(path_out, '%s.not_%s.%s' % (fn_prefix, idx_name, fn_type))
    fn_map_prefix = os.path.join(path_out, fn_prefix)
    fn_map_bam = fn_map_prefix + '.map_%s.bam' % idx_name
    fn_map_bed = fn_map_prefix + '.map_%s.bed' % idx_name
    fn_map_log = fn_map_prefix + '.map_%s.star.log' % idx_name
    ## skip exist files
    if os.path.exists(fn_map_bam) and overwrite is False:
        logging.info('file exists: %s' % fn_map_bam)
    else:
        c1 = 'STAR --runMode alignReads \
              --genomeDir %s \
              --readFilesIn %s \
              --readFilesCommand %s \
              --outFileNamePrefix %s \
              --runThreadN %s \
              --outFilterMismatchNoverLmax 0.07 \
              --outFilterMultimapNmax 1 \
              --limitOutSAMoneReadBytes 1000000 \
              --genomeLoad LoadAndRemove \
              --limitBAMsortRAM 10000000000 \
              --outSAMtype BAM SortedByCoordinate \
              --outReadsUnmapped Fastx'  % (idx, fn, freader, fn_map_prefix,
                                            multi_cores)
        p1 = subprocess.run(shlex.split(c1))
        # rename exists file
        os.rename(fn_map_prefix + 'Aligned.sortedByCoord.out.bam', fn_map_bam)
        os.rename(fn_map_prefix + 'Unmapped.out.mate1', fn_unmap_file)
        os.rename(fn_map_prefix + 'Log.final.out', fn_map_log)
        pysam.index(fn_map_bam)
        d = star_log_parser(fn_map_log)

    return [fn_map_bam, fn_unmap_file]



def map_se_batch(fn, idxes, path_out, para=1, multi_cores=1, aligner='STAR', 
                 overwrite=False):
    """
    mapping fastq to multiple indexes
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert isinstance(idxes, list)
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    if aligner.lower() == 'star':
        align_se = star_se
    elif aligner.lower() == 'bowtie2':
        align_se = bowtie2_se
    else:
        raise ValueError('unknown aligner: %s' % aligner)
    # iterate index
    fn_bam_files = []
    fn_input = fn
    for idx in idxes:
        para = 2 if idx is idxes[-1] else para
        fn_bam_idx, fn_unmap_idx = align_se(fn_input, idx, path_out, 
                                            para=para,
                                            multi_cores=multi_cores,
                                            overwrite=overwrite)
        fn_input = fn_unmap_idx
        fn_bam_files.append(fn_bam_idx)
    return fn_bam_files



def pcr_dup_remover(bam_in):
    """
    remove PCR duplicates using Picard
    """
    picard_jar = '/data/biosoft/picard/build/libs/picard.jar'
    if not os.path.exists(picard_jar):
        logging.error('file not found - picard.jar')
        return None
    bam_nodup = os.path.splitext(bam_in)[0] + '.nodup.bam'
    metrics_nodup = os.path.splitext(bam_in)[0] + '.nodup.metrics'
    log_nodup = os.path.splitext(bam_in)[0] + '.nodup.log'
    c1 = 'java -Xmx8g -jar {} MarkDuplicates INPUT={} OUTPUT={} METRICS_FILE={} \
          REMOVE_DUPLICATES=true ASSUME_SORTED=true'.format(picard_jar, 
          bam_in, bam_nodup, metrics_nodup)
    if os.path.exists(bam_in) and not os.path.exists(bam_nodup):
        with open(log_nodup, 'w') as fo:
            p1 = subprocess.run(shlex.split(c1), stdout = fo, stderr = fo)
        pysam.index(bam_nodup)
    return bam_nodup



def align(fns, smp_name, path_out, genome, spikein=None, multi_cores=1, 
          aligner='STAR', path_data=None, overwrite=False):
    """
    mapping reads to multiple indexes, one-by-one
    """
    assert isinstance(fns, list)
    assert isinstance(genome, str)
    assert isinstance(smp_name, str)

    # get indexes
    sg = idx_picker(genome, path_data=path_data, aligner=aligner)
    idxes = [sg]
    if isinstance(spikein, str) and not spikein == genome:
        sp = idx_picker(spikein, path_data=path_data, aligner=aligner) # 
        idxes.append(sp)
    idxes = list(filter(None.__ne__, idxes)) # idxes
    if len(idxes) == 0:
        raise ValueError('genome index not exists: %s' % path_data)

    # mapping se reads
    fn_bam_files = []
    # mapping 
    for fn in fns:
        logging.info('mapping file: %s' % fn)
        fn_prefix = file_prefix(fn)[0]
        fn_prefix = re.sub('\.clean|\.nodup|\.cut', '', fn_prefix)
        # fn_prefix = re.sub('_[12]$|_R[12]$', '', fn_prefix)
        path_out_fn = os.path.join(path_out, fn_prefix)
        b = map_se_batch(fn, idxes, path_out_fn, multi_cores=multi_cores,
                         aligner=aligner, overwrite=overwrite) # list
        fn_bam_files.append(b) # bam files
        rep_map_wrapper(path_out_fn)

    # merge bam files
    path_out_merge = os.path.join(path_out, smp_name)
    merge_bam_files = []
    if len(fn_bam_files) > 1:
        assert is_path(path_out_merge)
        for i in range(len(fn_bam_files[0])): # merge each sub-index
            se_bam_files = [b[i] for  b in fn_bam_files]
            merge_suffix = str_common(se_bam_files, suffix=True)
            merge_suffix = re.sub('^_[12]|_R[12]', '', merge_suffix)
            merge_bam_name = smp_name + merge_suffix
            merge_bam_file = os.path.join(path_out_merge, merge_bam_name)
            merge_bed_file = re.sub('.bam$', '.bed', merge_bam_file)
            if os.path.exists(merge_bam_file) and overwrite is False:
                logging.info('file exists: %s' % merge_bam_name)
            else:
                tmp = bam_merge(se_bam_files, merge_bam_file)
                pybedtools.BedTool(merge_bam_file).bam_to_bed().saveas(merge_bed_file)
            merge_bam_files.append(merge_bam_file)
        merge_map_wrapper(path_out_merge)
        fn_bam_files.append(merge_bam_files)

    # get genome mapping files (the last one)
    genome_bam_files = [f[-1] for f in fn_bam_files]

    # rename genome bam, to a shorter name
    # remove "not_index." suffix
    gbam_files = []
    gbed_files = []
    for i in range(len(genome_bam_files)):
        bam_from = genome_bam_files[i]
        bam_to = os.path.join(os.path.dirname(bam_from), 
                              filename_shorter(bam_from))
        if not os.path.exists(bam_to):
            os.symlink(os.path.basename(bam_from), bam_to)
        if not os.path.exists(bam_to + '.bai'):
            if not os.path.exists(bam_from + '.bai'):
                pysam.index(bam_from)
            os.symlink(os.path.basename(bam_from) + '.bai', 
                       bam_to + '.bai')
        gbam_files.append(bam_to)
        # rename .bed
        bed_from = re.sub('.bam$', '.bed', bam_from)
        bed_to = re.sub('.bam$', '.bed', bam_to)
        if os.path.exists(bed_from) and not os.path.exists(bed_to):
            os.symlink(os.path.basename(bed_from), bed_to)
        gbed_files.append(bed_to)

    return gbam_files # [gbam_files, gbed_files]


def main():
    args = get_args()
    fqs = [f.name for f in args.i]
    smp_name = args.n
    path_out = args.o
    genome = args.g
    multi_cores = args.p
    aligner = args.t.lower()
    rm_dup = args.rmdup
    path_data = args.path_data
    overwrite = args.overwrite
    # mapping
    p = map(fqs, smp_name, path_out, genome, 
        multi_cores=multi_cores, aligner=aligner, 
        path_data=path_data, overwrite=overwrite)
    p_out = p # bam files
    if args.rmdup:
        px = []
        # remove dup
        for b in p:
            x = pcr_dup_remover(b)
            px.append(x)
        p_out = px
    return p_out


if __name__ ==  '__main__':
    p = main()

## EOF
