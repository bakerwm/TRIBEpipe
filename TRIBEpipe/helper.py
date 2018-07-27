

import os
import sys
import re
# import datetime
import json
import glob
import argparse
import shlex
import subprocess
# import numpy as np
import pandas as pd
import binascii
import pysam
import logging
import gzip



## functions
def is_gz(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def seq_type(path, top_n=1000):
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
    f_reader = gzip.open if(is_gz(path)) else open
    with f_reader(path, 'rt') as f:
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


def is_idx(path, aligner='bowtie'):
    """
    check aligner index, bowtie, bowtie2, STAR
    """
    # bowtie/bowtie2
    c = [aligner + '-inspect', '-s', path]
    if aligner.lower() == 'star':
        pg = os.path.join(path, 'Genome')
        flag = True if os.path.exists(pg) else False
    else:
        p = subprocess.run(c, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
        flag = True if len(p) > 0 else False
    return flag



def idx_picker(genome, group='genome', path_data=None, aligner='bowtie'):
    """
    return the path of index
    group: genome, rRNA, tRNA, ...
    aligner: bowtie, bowie2
    #
    default: path ~/data/genome/
    """
    assert isinstance(genome, str)
    assert isinstance(group, str)
    if path_data is None:
        path_data = os.path.join(pathlib.Path.home(), 'data', 'genome')
    idx = os.path.join(path_data, genome, aligner + '_index', group)
    if aligner.lower() == 'star':
        idx = os.path.join(path_data, genome, 'STAR_index')
    if is_idx(idx, aligner):
        return idx
    else:
        return None



##############################
## wrapper results          ##
##############################

def bowtie2_log_parser(path):
    """
    Parsing the log file of bowtie
    fetch Input, unique, multiple, unmapped
    save in JSON format
    return dict of all values
    """
    logdict = {}
    _input = []
    _one_hit = []
    _multi_hit = []
    _not_hit = []
    _rpt = []
    with open(path, 'rt') as fi:
        for line in fi.readlines():
            line = line.strip()
            _num = line.split(' ')[0]
            _num = _num.strip('%')
            _num = float(_num)
            if line.endswith('reads; of these:'):
                _input.append(_num)
            elif ') aligned 0 times' in line:
                _not_hit.append(_num)
            elif ') aligned exactly 1 time' in line:
                _one_hit.append(_num)
            elif ') aligned >1 times' in line:
                _multi_hit.append(_num)
            else:
                continue
    # save to dict
    logdict['input_reads'] = int(_input[0]) # first one
    logdict['mapped'] = int(sum(_one_hit + _multi_hit))
    logdict['unique'] = int(sum(_one_hit))
    logdict['multi'] = int(sum(_multi_hit))
    logdict['unmapped'] = int(_not_hit[-1]) # -m suppress
    logdict['map_pct'] = '{:.2f}%'.\
        format(int(logdict['mapped']) / int(logdict['input_reads'])*100)
    json_out = os.path.splitext(path)[0] + '.json'
    with open(json_out, 'w') as fo:
        json.dump(logdict, fo, indent=4)
    return logdict



def star_log_parser(path):
    logdict = {}
    with open(path) as f:
        for line in f:
            sep = line.strip().split('|')
            if 'Number of input reads' in line:
                logdict['input_reads'] = int(sep[1].strip())
            elif 'Uniquely mapped reads number' in line:
                logdict['unique'] = int(sep[1].strip())
            elif 'Number of reads mapped to multiple loci' in line:
                logdict['multi'] = int(sep[1].strip())
            else:
                pass
    logdict['mapped'] = logdict['unique'] + logdict['multi']
    logdict['unmapped'] = logdict['input_reads'] - logdict['mapped']
    logdict['map_pct'] = '{:.2f}%'.format(logdict['mapped'] / logdict['input_reads'] * 100)
    json_out = os.path.splitext(path)[0] + '.json'
    with open(json_out, 'wt') as fo:
        json.dump(logdict, fo, indent=4)
    return logdict




def rep_map_wrapper(path, save=True):
    """
    wrap all bowtie log files, [only for this script] namespace 
    summarize mapping and RTStops 
    header: name, group, read, RTStop
    input: list of json files
    output: pd.DataFrame
    """
    def _json_wrapper(fn): 
        """
        Only for bowtie map statistics
        parsing only one json file
        output: type, count
        """
        group = fn.split('.')[-3] # group name, *map_genome.bowtie.json
        group = group.split('_')[1] # reference name
        name = os.path.splitext(os.path.basename(fn))[0]
        with open(fn, 'r') as f:
            da = json.load(f)
        df = [name, group, da['input_reads'], da['mapped'], 
              da['unmapped']]
        return df

    # multiple json files
    json_files = sorted(glob.glob(path + '/*.json'), key=len)
    rep_prefix = os.path.basename(path)
    df = pd.DataFrame(columns=['name', 'group', 'read'])
    # 
    for n in range(len(json_files)):
        _, g, _, m1, _ = _json_wrapper(json_files[n])
        if n == 0 and len(json_files) > 1 and g == 'genome':
            g = 'spikein' # the first one - genome
        df = df.append(pd.DataFrame([[rep_prefix, g, m1]], 
            columns = ['name', 'group', 'read']), ignore_index=True)
    _, gx, _, _, un = _json_wrapper(json_files[-1]) # unmap
    df = df.append(pd.DataFrame([[rep_prefix, 'unmapped', un]], 
                   columns=['name', 'group', 'read']), ignore_index=True)
    save_csv = os.path.join(os.path.dirname(path), 
                            rep_prefix + '.mapping_stat.csv')
    if save:
        df.to_csv(save_csv, ',', header=True, index=False)
    return df



def merge_map_wrapper(path, save=True):
    """
    count BAM files
    Output: pd.DataFrame
    """
    bam_files = sorted(glob.glob(path + '/*.bam'), key=len) # bam files
    merge_prefix = os.path.basename(path)
    df = pd.DataFrame(columns=['name', 'group', 'read'])
    bam_files = [f for f in bam_files if not os.path.islink(f)]
    # iterate
    for n in range(len(bam_files)):
        b_cnt = pysam.AlignmentFile(bam_files[n], 'rb').count()
        group = bam_files[n].split('.')[-2] # group name*.map_genome.bam
        group = group.split('_')[1] # reference name
        if n == 0 and len(bam_files) > 1 and group == 'genome':
            group = 'spikein' # the first one - genome
        dfx = pd.DataFrame([[merge_prefix, group, b_cnt]],
                            columns=['name', 'group', 'read'])
        df = df.append(dfx, ignore_index=True)
    save_csv = os.path.join(os.path.dirname(path), 
                            merge_prefix + '.mapping_stat.csv')
    if save:
        df.to_csv(save_csv, ',', header=True, index=False)
    return df



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

