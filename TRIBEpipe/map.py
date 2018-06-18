#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mapping reads to reference genome
1. to spike-in genome
2. to repbase
3. to rRNA, tRNA, MT
4. to genome

"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-06-13"
__version__ = "0.1"

import os, sys, re, datetime, json, glob
import argparse, shlex, subprocess
import numpy as np
import pandas as pd
import binascii
import pysam
import logging

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
    parser.add_argument('-o', required = True, default = None, metavar = 'OUTPUT',  
        help = 'The directory to save results.')
    parser.add_argument('-g', default = 'hg19', metavar = 'GENOME',
        help = 'Reference genome : dm3, dm6, hg19, GRCh38, mm10, GRCm38')
    parser.add_argument('-x', required = True, metavar = 'index',
        help = 'Index of the genome for alignment tools')
    parser.add_argument('-t', required = True, default = 'bowtie', metavar = 'Aligner',
        help = 'Aligner for the mapping, bowtie, bowtie2, STAR')
    parser.add_argument('-p', default = 1, metavar = 'Threads', type = int, 
        help = 'Number of threads to launch, default [1].')
    parser.add_argument('--rmdup', action = 'store_true',
        help = 'remove PCR duplicates using Picard, if specified')
    args = parser.parse_args()
    return args


## functions
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def seq_type(path, top_n = 1000):
    """
    Check the input file is FASTA or FASTQ format
    support only, squence in one line
    extract the first 1000 lines
    '>' is fasta, '@' is fastq, 'others' is None
    """
    tag = set()
    if not isinstance(path, str):
        logging.info('only string accepted, error.')
        return None
    f_reader = gzip.open if(is_gz_file(path)) else open
    with f_reader(path, 'r') as f:
        for i, line in enumerate(f):
            if i > top_n:
                break
            elif i % 4 == 0:
                tag.add(line[0])
            else:
                continue
    if tag ==  {'@'}:
        return 'fastq'
    elif tag ==  {'>'}:
        return 'fasta'
    else:
        return None


def bam_merge(bam_ins, bam_out):
    """
    merge multiple bam files
    input: list of bam files
    input: out.bam
    """
    # check input files
    bam_flag = []
    for b in bam_ins:
        if not os.path.exists(b) is True:
            bam_flag.append(b)
    if len(bam_flag) > 0:
        sys.exit('BAM files not exists:' + '\n'.join(bam_flag))
    # check output file
    if os.path.exists(bam_out) is True:
        pass
        # sys.exit('BAM exists:' + bam_out)
    else:
        # merge
        pysam.merge('-f', bam_out + '.unsorted.bam', *bam_ins) # overwrite output BAM
        pysam.sort('-o', bam_out, bam_out + '.unsorted.bam')
        pysam.index(bam_out)
        os.remove(bam_out + '.unsorted.bam')
        

def aligner_index_validator(path, aligner = 'bowtie'):
    """validate the alignment index, bowtie, bowtie2, hisat2"""
    # bowtie index
    inspect = aligner + '-inspect'
    cmd = [inspect, '-s', path]
    n = subprocess.run(cmd, check = False, stdout = subprocess.PIPE).stdout
    if len(n) > 0:
        return True
    else:
        return False


def bowtie_map_se(read_in, align_index, out_path = None, *, para = 1, 
                  multi_cores = 1):
    """
    mapping single read to one index using bowtie
    Input: fastq|a
    Output: bam (sorted), unmapped reads
    """
    aligner_index_validator(align_index, 'bowtie')
    freader = 'zcat' if is_gz_file(read_in) else 'cat'
    mypara = '-v 2 -m 1 --strata' if para == 1 else '-v 2 -k 100'# unique mapper
    # # outdir exists
    if out_path is None:
        out_path = os.path.dirname(read_in)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    ## 
    read_type = seq_type(read_in)
    if read_type is 'fasta':
        mypara +=  ' -f'
    elif seq_type(read_in) is 'fastq':
        mypara +=  ' -q'
    else:
        raise Exception('unkown type of read file(s)')
    ## filename
    idx_prefix = os.path.basename(align_index)
    read_prefix = os.path.basename(read_in) # remove *.nodup.clean    
    read_prefix = re.sub('\.f[astq]+(.gz)?', '', read_prefix)
    read_prefix = re.sub('\.clean|\.nodup|\.cut\-?\d+', '', read_prefix) #
    read_unmap_name = '{}.unmapped.{}'.format(read_prefix, idx_prefix, read_type)
    read_unmap = os.path.join(out_path, read_unmap_name)
    bam_out = os.path.join(out_path, read_prefix + '.bam')
    log_out = re.sub('.bam', '.log', bam_out)
    ## skip exist files
    if os.path.exists(bam_out) and os.path.exists(read_unmap):
        logging.info('BAM file exists, mapping skipped...: ' + bam_out)
    else:
        c0 = '{} {}'.format(freader, read_in)
        c1 = 'bowtie {} -p {} --mm --best --sam --no-unal --un {} {} {}'.\
            format(mypara, multi_cores, read_unmap, align_index, read_in)
        c2 = 'samtools view -bhS -F 4 -@ {} -'.format(multi_cores)
        c3 = 'samtools sort -@ {} -'.format(multi_cores)
        c4 = 'samtools index {}'.format(bam_out)
        with open(bam_out, 'w') as fbam, open(log_out, 'w') as fo:
            p0 = subprocess.Popen(shlex.split(c0), stdout = subprocess.PIPE)
            p1 = subprocess.Popen(shlex.split(c1), stdin = p0.stdout, stdout = subprocess.PIPE, stderr = fo)
            p2 = subprocess.Popen(shlex.split(c2), stdin = p1.stdout, stdout = subprocess.PIPE)
            p3 = subprocess.Popen(shlex.split(c3), stdin = p2.stdout, stdout = fbam)
            px = p3.communicate()
            p4 = subprocess.run(shlex.split(c4)) # index
        d = bowtie_log_parser(log_out)

    return [bam_out, read_unmap]


def bowtie2_map_se(read_in, align_index, out_path = None, *, para = 1, 
                  multi_cores = 1):
    """
    mapping single read to one index using bowtie2
    Input: fastq|a
    Output: bam (sorted), unmapped reads
    """
    aligner_index_validator(align_index, 'bowtie2')
    freader = 'zcat' if is_gz_file(read_in) else 'cat'
    mypara = '--very-sensitive' if para == 1 else '--local'# unique mapper
    # # outdir exists
    if out_path is None:
        out_path = os.path.dirname(read_in)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    ## 
    read_type = seq_type(read_in)
    if read_type is 'fasta':
        mypara +=  ' -f'
    elif seq_type(read_in) is 'fastq':
        mypara +=  ' -q'
    else:
        raise Exception('unkown type of read file(s)')
    ## filename
    idx_prefix = os.path.basename(align_index)
    read_prefix = os.path.basename(read_in) # remove *.nodup.clean    
    read_prefix = re.sub('\.f[astq]+(.gz)?', '', read_prefix)
    read_prefix = re.sub('\.clean|\.nodup|\.cut\-?\d+', '', read_prefix) #
    read_unmap_name = '{}.unmapped.{}'.format(read_prefix, idx_prefix, read_type)
    read_unmap = os.path.join(out_path, read_unmap_name)
    bam_out = os.path.join(out_path, read_prefix + '.bam')
    log_out = re.sub('.bam', '.log', bam_out)
    ## skip exist files
    if os.path.exists(bam_out) and os.path.exists(read_unmap):
        logging.info('BAM file exists, mapping skipped...: ' + bam_out)
    else:
        c0 = '{} {}'.format(freader, read_in)
        c1 = 'bowtie {} -p {} --mm --best --sam --no-unal --un {} {} -'.\
            format(mypara, multi_cores, read_unmap, align_index, read_in)
        c2 = 'samtools view -bhS -F 4 -@ {} -'.format(multi_cores)
        c3 = 'samtools sort -@ {} -'.format(multi_cores)
        c4 = 'samtools index {}'.format(bam_out)
        with open(bam_out, 'w') as fbam, open(log_out, 'w') as fo:
            p0 = subprocess.Popen(shlex.split(c0), stdout = subprocess.PIPE)
            p1 = subprocess.Popen(shlex.split(c1), stdin = p0.stdout, stdout = subprocess.PIPE, stderr = fo)
            p2 = subprocess.Popen(shlex.split(c2), stdin = p1.stdout, stdout = subprocess.PIPE)
            p3 = subprocess.Popen(shlex.split(c3), stdin = p2.stdout, stdout = fbam)
            px = p3.communicate()
            p4 = subprocess.run(shlex.split(c4)) # index
        d = bowtie2_log_parser(log_out)

    return [bam_out, read_unmap]


def star_map_se(read_in, align_index, out_path = None, *, para = 1, 
                multi_cores = 1):
    """
    mapping single read to one index using STAR
    Input: fastq|a
    Output: bam (sorted), unmapped reads
    """
    freader = 'zcat' if is_gz_file(read_in) else 'cat'
    mypara = '--very-sensitive' if para == 1 else '--local'# unique mapper
    # # outdir exists
    if out_path is None:
        out_path = os.path.dirname(read_in)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    ## 
    read_type = seq_type(read_in)
    if read_type is 'fasta':
        mypara +=  ' -f'
    elif seq_type(read_in) is 'fastq':
        mypara +=  ' -q'
    else:
        raise Exception('unkown type of read file(s)')
    ## filename
    read_prefix = re.sub('.gz$', '', os.path.basename(read_in))
    read_prefix = os.path.splitext(read_prefix)[0]
    read_prefix = re.sub('\.clean|\.nodup|\.cut\-?\d+', '', read_prefix) #
    read_unmap_name = '{}.unmapped.{}'.format(read_prefix, read_type)
    read_unmap = os.path.join(out_path, read_unmap_name)
    bam_out = os.path.join(out_path, read_prefix + '.bam')
    log_out = re.sub('.bam', '.log', bam_out)
    ## skip exist files
    if os.path.exists(bam_out) and os.path.exists(read_unmap):
        logging.info('BAM file exists, mapping skipped...: ' + bam_out)
    else:
        c0 = 'STAR --genomeDir {} --readFilesIn {} \
              --readFilesCommand {} --outFileNamePrefix {} \
              --runThreadN {} '.format(
                align_index, read_in, 
                freader, os.path.splitext(bam_out)[0],
                multi_cores)
        c1 = c0 + '--runMode alignReads \
             --limitOutSAMoneReadBytes 1000000 \
             --genomeLoad LoadAndKeep \
             --limitBAMsortRAM 10000000000 \
             --outSAMtype BAM SortedByCoordinate \
             --outReadsUnmapped Fastx \
             --outFilterMismatchNoverLmax 0.05 \
             --seedSearchStartLmax 20'
        c2 = 'samtools index {}'.format(bam_out)
        with open(bam_out, 'w') as fbam, open(log_out, 'w') as fo:
            p1 = subprocess.run(shlex.split(c1), stdout = fo)
            ## rename files
            os.rename(os.path.splitext(bam_out)[0] + 'Aligned.sortedByCoord.out.bam', bam_out)
            os.rename(os.path.splitext(bam_out)[0] + 'Unmapped.out.mate1', read_unmap)
            os.rename(os.path.splitext(bam_out)[0] + 'Log.final.out', log_out)
            ## make index
            p2 = subprocess.run(shlex.split(c2))
        d = star_log_parser(log_out)

    return [bam_out, read_unmap]


def pcr_dup_remover(bam_in):
    """
    remove PCR duplicates using Picard
    """
    picard_jar = '/data/biosoft/picard/build/libs/picard.jar'
    if not os.path.exists(picard_jar):
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
    return bam_nodup


##############################
## wrapper results          ##
##############################
def bowtie_log_parser(fn):
    """
    Parsing the log file of bowtie (to stderr)
    fetch Input, mapped, unmapped reads
    save in JSON format
    return dict of all values
    """
    # !!! to-do
    # convert to DataFrame, json, dict
    # unmap = input - reported
    logdict = {}
    _input = []
    _one_hit = []
    _not_hit = []
    _rpt = []
    with open(fn, 'r') as f:
        for line in f.readlines():
            if 'reads processed' in line:
                _num1 = re.sub(',', '', line.split()[-1])
                _input.append(int(_num1))
            elif 'at least one reported' in line:
                _num2 = re.sub(',', '', line.split()[-2])
                _one_hit.append(int(_num2))
            elif 'Reported' in line:
                _num4 = re.sub(',', '', line.split()[1])
                _rpt.append(int(_num4))
            elif 'failed to align' in line:
                _num3 = re.sub(',', '', line.split()[-2])
                _not_hit.append(int(_num3))
            else:
                continue
    logdict['input_reads'] = _input[0] # first one
    # logdict['mapped'] = sum(_one_hit) # sum all
    logdict['mapped'] = sum(_rpt) # sum all
    logdict['unmapped'] = int(_input[-1]) - int(_one_hit[-1]) # -m suppress
    logdict['map_pct'] = '{:.2f}%'.\
        format(int(logdict['mapped']) / int(logdict['input_reads'])*100)
    json_out = os.fn.splitext(fn)[0] + '.json'
    with open(json_out, 'w') as fo:
        json.dump(logdict, fo, indent = 4)
    return logdict


def bowtie2_log_parser(fn):
    """
    Parsing the log file of bowtie (to stderr)
    fetch Input, mapped, unmapped reads
    save in JSON format
    return dict of all values
    """
    # !!! to-do
    # convert to DataFrame, json, dict
    # unmap = input - reported
    logdict = {}
    _input = []
    _one_hit = []
    _multi_hit = []
    _not_hit = []
    with open(fn, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            num = int(line.split()[0])
            if 'reads; of these' in line:
                _input.append(num)
            elif 'aligned exactly 1 time' in line:
                _one_hit.append(num)
            elif 'aligned >1 times' in line:
                _multi_hit.append(num)
            elif 'aligned 0 times' in line:
                _not_hit.append(num)
            else:
                continue
    # group
    logdict['input_reads'] = _input[0] # first one
    logdict['unique'] = _one_hit[0]
    logdict['multiple'] = _multi_hit[0]
    logdict['mapped'] = _input[0] - _not_hit[0]
    logdict['unmapped'] = _not_hit[0]
    logdict['map_pct'] = '{:.2f}%'.\
        format(int(logdict['mapped']) / int(logdict['input_reads'])*100)
    json_out = os.fn.splitext(fn)[0] + '.json'
    with open(json_out, 'w') as fo:
        json.dump(logdict, fo, indent = 4)
    return logdict



def star_log_parser(fn):
    """
    Parsing the log file of star (to stderr)
    fetch Input, mapped, unmapped reads
    save in JSON format
    return dict of all values
    """
    # !!! to-do
    # convert to DataFrame, json, dict
    # unmap = input - reported
    logdict = {}
    _input = []
    _one_hit = []
    _multi_hit = []
    with open(fn, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            # num = line.split()[-1]
            if 'Number of input reads' in line:
                _input.append(int(line.split()[-1]))
            elif 'Uniquely mapped reads number' in line:
                _one_hit.append(int(line.split()[-1]))
            elif 'Number of reads mapped to multiple loci' in line:
                _multi_hit.append(int(line.split()[-1]))
            else:
                continue
    # group
    logdict['input_reads'] = _input[0] # first one
    logdict['unique'] = _one_hit[0]
    logdict['multiple'] = _multi_hit[0]
    logdict['mapped'] = _one_hit[0] + _multi_hit[0]
    logdict['unmapped'] = _input[0] - _one_hit[0] - _multi_hit[0]
    logdict['map_pct'] = '{:.2f}%'.\
        format(int(logdict['unique'] + int(logdict['multiple'])) / int(logdict['input_reads']) * 100)
    json_out = os.path.splitext(fn)[0] + '.json'
    with open(json_out, 'w') as fo:
        json.dump(logdict, fo, indent = 4)
    return logdict



def rep_map_wrapper(fn, save = True):
    """
    wrap all bowtie log files, [only for this script] namespace 
    summarize mapping and RTStops 
    header: name, group, read, RTStop
    input: list of json files
    output: pd.DataFrame
    """
    def _json_wrapper(file):
        """
        parsing only one json file
        output: type, count
        """
        try:
            group = file.split('.')[-2]
            group = re.sub('map_', '', group)
            name = os.path.basename(file).split('.')[0]
            with open(file, 'r') as f:
                da = json.load(f)
            return [name, group, da['input_reads'], da['mapped'], 
                    da['unmapped']]
        except IOError:
            print('json file faild: ' + file)
    ## sort files by len
    json_files = sorted(glob.glob(fn + '/*.json'), key = len)
    rep_prefix = os.path.basename(fn)
    df = pd.DataFrame(columns = ['name', 'group', 'read'])

    #for ff in json_files:
    for n in range(len(json_files)):
        _, g, _, c, _ = _json_wrapper(json_files[n])
        # the first one - genome
        g = 'spikein' if g == 'genome' and n == 0 else g
        df = df.append(pd.DataFrame([[rep_prefix, g, c]], 
            columns = ['name', 'group', 'read']), ignore_index = True)
    n1, g1, _, _, c1 = _json_wrapper(json_files[-1]) # unmap
    df = df.append(pd.DataFrame([[n1, 'unmapped', c1]], 
            columns = ['name', 'group', 'read']), ignore_index = True)
    save_csv = os.path.join(os.path.dirname(fn), rep_prefix + '.mapping_stat.csv')
    if save:
        try:
            df.to_csv(save_csv, ',', header = True, index = False)
        except IOError:
            print('[failed] saving data to file: ' + save_csv)
    return df



def merge_map_wrapper(fn, save = True):
    """
    count BAM files
    Output: pd.DataFrame
    """
    b_files = sorted(glob.glob(fn + '/*.bam'), key = len) # bam files
    b_prefix = os.path.basename(fn)
    df = pd.DataFrame(columns = ['name', 'group', 'read'])
    #for b in b_files:
    for n in range(len(b_files)):
        b_cnt = pysam.AlignmentFile(b_files[n], 'rb').count()
        group = os.path.basename(b_files[n]).split('.')[-2] #
        group = re.sub('^map_', '', group)
        group = 'spikein' if n == 0 and group == 'genome' else group
        df = df.append(pd.DataFrame([[b_prefix, group, b_cnt]],
                       columns = ['name', 'group', 'read']),
                       ignore_index = True)
    save_csv = os.path.join(os.path.dirname(fn), b_prefix + '.mapping_stat.csv')
    if save:
        try:
            df.to_csv(save_csv, ',', header = True, index = False)
        except IOError:
            print('[failed] saving data to file: ' + save_csv)
    return df


def map(fqs, smp_name, out_path, genome, multi_cores, genome_index, 
        map_tools = 'bowtie'):
    """
    mapping reads to multiple indexes, one-by-one
    call RT stops
    merge RT stops 
    filter merged_rt by snoRNA, miRNA, ...
    Input: read1, read2, ... [1 to 4]
    Input: idx1, idx2, ...
    Output: RTStops
    """    
    # idxs = aligner_index_picker(genome, spikein)
    bam_files = []
    logging.info('mapping reads')
    map_dict = {'bowtie': bowtie_map_se,
                'bowtie2': bowtie2_map_se,
                'star': star_map_se}
    mapper = map_dict[map_tools]
    for read in fqs:
        read_prefix = re.sub(r'.f[astq]*(.gz)?', '', os.path.basename(read))
        read_prefix = re.sub('\.clean|\.nodup|\.cut\-?\d+', '', read_prefix) #
        path_sub = os.path.join(out_path, read_prefix)
        p = mapper(read, genome_index, path_sub, multi_cores = multi_cores)
        rep_map_wrapper(path_sub) # wrapper
        bam_files.append(p[0])

    # merge bam files
    path_merge = os.path.join(out_path, smp_name)
    if not os.path.exists(path_merge):
        os.makedirs(path_merge)
    merge_bam_files = [] # map to multiple indexes
    if len(bam_files) > 1: # multiple replicates
        merged_bam = os.path.join(path_merge, smp_name + '.bam')
        if os.path.exists(merged_bam):
            logging.info('BAM file exists, merging skipped...: ' + merged_bam)
        else:
            tmp = bam_merge(bam_files, merged_bam)
        merge_bam_files.append(merged_bam)

        merge_map_wrapper(path_merge) # wrapper
        bam_files.append(merge_bam_files[0])

    return bam_files


def main():
    args = get_args()
    fqs = [f.name for f in args.i]
    smp_name = args.n
    out_path = args.o
    genome = args.g
    genome_index = args.x
    multi_cores = args.p
    map_tools = args.t.lower()
    rm_dup = args.rmdup
    # mapping
    p = map(fqs, smp_name, out_path, genome, multi_cores, genome_index, map_tools = map_tools)
    p_out = p
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
