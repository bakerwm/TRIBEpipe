

import os
import sys
import re
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
import pathlib


##
def xopen(fn, mode='r', bgzip=False):
    """
    Read / Write regular and gzip file, also support stdin
    """
    assert isinstance(fn, str)
    if fn == '-':
        return sys.stdin if 'r' in mode else sys.stdout
    if fn.endswith('.gz') and mode.startswith('w') or is_gz(fn):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)



## functions
def is_gz(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def is_path(path, create = True):
    """
    Check path, whether a directory or not
    if not, create it
    """
    assert isinstance(path, str)
    if os.path.exists(path):
        return True
    else:
        if create:
            try:
                os.makedirs(path)
                return True
            except IOError:
                logging.error('failed to create directories: %s' % path)
        else:
            return False


def seq_type(fn, top_n = 1000):
    """
    Check the top 1000 rows of fn
    identify @ for fastq, > for fasta, * unknown
    """
    assert isinstance(fn, str)
    tag = set()
    with xopen(fn, 'rt') as fi:
        for i, line in enumerate(fi):
            if i > top_n:
                break
            elif i % 4 == 0:
                b = line[0] # the first base
                if b.lower() in 'acgtn':
                    continue
                else:
                    tag.add(line[0])
            else:
                continue
    if tag ==  {'@'}:
        return 'fastq'
    elif tag ==  {'>'}:
        return 'fasta'
    else:
        return None


def is_fastq(fn):
    if seq_type(fn) == 'fastq':
        return True
    else:
        return False


def is_fasta(fn):
    if seq_type(fn) == 'fasta':
        return True
    else:
        return False


def str_common(strList, suffix = False):
    # extract longest prefix/suffix from list of strings
    # default: prefix
    # sort strings by len
    def iterStop(exp):
        if exp is False:
            raise StopIteration
        else:
            return True    

    def commonPrefix(s1, s2):
        # prefix
        return ''.join(list(val for i, val in enumerate(s1) 
                       if iterStop(s2[i] is val)))

    def fact(l):
        if len(l) ==  1:
            return l[0]
        else:
            la = l[0:2]
            lb = l[2:]
            s = commonPrefix(la[0], la[1])
            lb.insert(0, s)
            return fact(lb)

    ## empty or single item 
    if len(strList) ==  0:
        return ''
    elif len(strList) ==  1:
        return strList[0]
    else:
        ## save a copy of list
        L2 = sorted(strList, key = len)
        c = fact(L2)
    
    ## suffix, reverse strings
    if suffix is True:
        L2 = [i[::-1] for i in L2]
        c = fact(L2)
        c = c[::-1]

    return c # string 0-index


def file_prefix(fn, with_path = False):
    """
    extract the prefix of a file
    remove extensions
    .gz, .fq.gz
    """
    assert isinstance(fn, str)
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    if px.endswith('gz') or px.endswith('.bz'):
        px = os.path.splitext(p1)[1] + px
        p1 = os.path.splitext(p1)[0]
    if not with_path:
        p1 = os.path.basename(p1)
    return [p1, px]


def rm_suffix1(fn):
    """
    simplify the name of bam files
    from: {name}.not_{}.not_{}.....map_{}
    to: {name}
    """
    if '.' in fn:
        p = os.path.splitext(fn)[0]
        px = os.path.splitext(fn)[1]
        if px.startswith('.not_') or px.startswith('.map_'):
            return rm_suffix1(p)
        else:
            return fn
    else:
        return fn
        

def filename_shorter(fn, with_path=False):
    """
    input: name1.not_spikein.not_mtrRNA.map_genome.bam 
           name2.not_spikein.not_mtrRNA.map_genome.bam
    output: name1.bam
            name2.bam
    """
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    p2 = rm_suffix1(p1)
    if not with_path:
        p2 = os.path.basename(p2)
    return p2 + px


def file_row_counter(fn):
    """
    count the file rows
    count '\n' 
    from @glglgl on stackoverflow, modified
    https://stackoverflow.com/a/9631635/2530783
    """
    def blocks(files, size = 1024 * 1024):
        while True:
            b = files.read(size)
            if not b: break
            yield b
    freader = gzip.open if is_gz(fn) else open
    with freader(fn, 'rt', encoding="utf-8", errors='ignore') as fi:
        return sum(bl.count('\n') for bl in blocks(fi))


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
    bam_files = [f for f in bam_files if not f.endswith('nodup.bam')]
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


class Genome_info():
    """
    including the information of genome
    index, annotation, ...
    """

    def __init__(self, genome, **kwargs):
        assert isinstance(genome, str)
        self.kwargs = kwargs
        self.kwargs['genome'] = genome
        if not 'path_data' in kwargs:
            self.kwargs['path_data'] = os.path.join(pathlib.Path.home(), 
                                                    'data', 'genome')
        

    def get_fa(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        gfa = os.path.join(path_data, genome, 'bigZips', genome + '.fa')
        assert os.path.exists(gfa)
        return gfa


    def get_fasize(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        gsize = os.path.join(path_data, genome, 'bigZips', genome + '.chrom.sizes')
        assert os.path.exists(gsize)
        return gsize


    def bowtie_index(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        return idx_picker(genome, path_data=path_data, aligner='bowtie')


    def bowtie2_index(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        return idx_picker(genome, path_data=path_data, aligner='bowtie2')


    def hisat2_index(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        return idx_picker(genome, path_data=path_data, aligner='hisat2')


    def star_index(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        return idx_picker(genome, path_data=path_data, aligner='star')


    def phylop100(self):
        """
        only support hg19
        """
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        phylop100 = os.path.join(self.kwargs['path_data'],
                            genome, 'phyloP100way', 
                            genome + '.100way.phyloP100way.bw')
        if not os.path.exists(phylop100):
            phylop100 = None
        return phylop100

        
    def gene_bed(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        gbed = os.path.join(path_data, 
                            genome,
                            'annotation_and_repeats', 
                            genome + '.refseq.bed')
        if not os.path.exists(gbed):
            gbed = None
        return gbed


    def gene_rmsk(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        grmsk= os.path.join(path_data, 
                            genome,
                            'annotation_and_repeats', 
                            genome + '.rmsk.bed')
        if not os.path.exists(grmsk):
            grmsk = None
        return grmsk


    def ensembl_gtf(self):
        genome = self.kwargs['genome']
        path_data = self.kwargs['path_data']
        ggtf = os.path.join(path_data, 
                            genome,
                            'annotation_and_repeats', 
                            genome + '.ensembl.gtf')
        if not os.path.exists(ggtf):
            ggtf = None
        return ggtf

