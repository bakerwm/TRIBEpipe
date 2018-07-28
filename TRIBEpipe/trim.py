#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Trimming FASTQ reads, SE an PE
 - 3' adapter
 - low quality bases at both ends
 - for NSR, trim 7 nt at the beginning of reads
 - poly(N) sequences 
 - (optional) barcode at 5' or 3' end  
 - (optional) additional adpaters  

## To-do
 - add subparser / subprogram
 https://docs.python.org/3/library/argparse.html

## ad feature
--cut_before_trim
--cut_after_trim

# NSR: --cut_before_trim 10  --cut_after_trim 10
# RNAseq: --cut_before_trim 0 --cut_after_trim 0


"""

__author__ = "Ming Wang <wangm08@hotmail.com>"
__copyright__ = "2018 by Ming Wang <wangm08@hotmail.com>"
__license__ = "MIT"
__email__ = "wangm08@hotmail.com"
__version__ = "0.1"

import os
import sys
import re
import argparse
import subprocess
import shlex
import logging
from TRIBEpipe.helper import *
logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'trimmer',
        description = 'Trimming reads by cutadapt',
        epilog = 'Example: trimmer -i input.fq')
    parser.add_argument('-i', nargs = '+', required = True, metavar = 'input', 
        type = argparse.FileType('r'),
        help = 'Sequencing reads in FASTQ format, support (*.gz), 1-4 files.')
    parser.add_argument('-a', 
        default = 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG', 
        metavar = 'adapter',
        help = '3 prime adapter sequence, \
        default [AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG].')
    parser.add_argument('-o', default = None, metavar = 'out_path', 
        help = 'The directory to save results.')
    parser.add_argument('--read12', type=int, default=1, metavar='read12',
        help='which one of PE reads, 1=read1, 2=read2, default: 1')
    parser.add_argument('-m', default = 15, metavar = 'len_min', 
        type = int, help = 'Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('-q', default = 20, metavar = 'quality', 
        type = int, help = 'The cutoff of base quality, default [20]')
    parser.add_argument('-p', default = 80, metavar = 'PERCENT', type = int,
        help = 'minimum percent of bases that must have -q quality, default [80]')
    parser.add_argument('--rm_untrim', action='store_false',
        help='if specified, discard reads without adapter')
    parser.add_argument('--rm_dup', action='store_true',
        help='if specified, remove duplicated reads' )
    parser.add_argument('--cut_before_trim', default='0', metavar='cut1', 
        help='cut bases before trimming adapter, Number of bases to cut from each \
        read, plus on 5-prime end, minus on 3-prime end, could be \
        single, or double numbers, eg: 3 or -4 or 3,-4, default [0]')
    parser.add_argument('--cut_after_trim', default='0', metavar='cut2', 
        help='cut bases after trimming adapter, Number of bases to cut from each \
        read, plus on 5-prime end, minus on 3-prime end, , could be \
        single, or double numbers, eg: 3 or -4 or 3,-4, default [0]')
    parser.add_argument('--cut_before_trim', default=0,
        help = 'cut N bases from each read before trimming adapters, plus value at 5 end, \
                minus value at 3 end, ')
    parser.add_argument('--threads', default = 1,  metavar = 'THREAD', type = int,
        help = 'Number of threads to launch, default [1]')
    parser.add_argument('-e', default = 0.1, metavar = 'err_rate', type = float,
        help = 'Maximum allowed error rate, default [0.1]')
    parser.add_argument('-O', default = 1, metavar = 'OVERLAP',
        help = 'Required N bases overlap between reads and adapter, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    args = parser.parse_args()
    return args



def adapter_chopper(s, step=2, window=15):
    """
    chop the adapter by given length
    a series of adapters to trim
    """
    assert isinstance(s, str)
    assert isinstance(step, int)
    assert isinstance(window, int)
    p = []
    if len(s) <= window:
        p.append(s)
    else:
        for i in range(int(len(s) / step)):
            a = i * step
            b = a + window
            #b = b if b < len(s) else len(s)
            if b > len(s):
                continue
            p.append(s[a:b])
    return p


def dup_remover(fn, path_out, *, q=20, p=80):
    """
    Remove duplicate fastq sequences using fastx-collapser (fq to fa)
    convert fa to fastq
    optional:
    trim N-bases
    """
    pkg_dir, _ = os.path.split(goldclip.__file__)
    fa2fq = os.path.join(pkg_dir, 'bin', 'fasta_to_fastq.pl')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    ## q, p, m
    fn_out_name = file_prefix(fn)[0] + '.nodup.fastq'
    fn_out_file = os.path.join(path_out, fn_out_name)
    if not os.path.exists(fn_out_file):
        freader = 'zcat' if is_gz(fn) else 'cat'
        c1 = '{} {}'.format(freader, fn)
        c2 = 'fastq_quality_filter -Q33 -q {} -p {}'.format(q, p)
        c3 = 'fastx_collapser -Q33'
        c4 = 'perl {} -'.format(fa2fq)
        cmd1 = shlex.split(c1)
        cmd2 = shlex.split(c2)
        cmd3 = shlex.split(c3)
        cmd4 = shlex.split(c4)
        p1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdin = p1.stdout, stdout = subprocess.PIPE)
        p3 = subprocess.Popen(cmd3, stdin = p2.stdout, stdout = subprocess.PIPE)
        with open(fn_out_file, 'wt') as fo:
            p4 = subprocess.Popen(cmd4, stdin = p3.stdout, stdout = fo)
            p5 = p4.communicate()
    return fn_out_file


def cutadapt_cut(s, cut_para=True):
    """
    recognize para: cut for cutadapt
    eg: cut=6, cut=-3, cut=6,-3
    """
    if ',' in s:
        n = s.split(',')
        if len(n) > 2:
            raise ValueError('illegal ad_cut: %s' % s)
        else:
            c1, c2 = (int(n[0]), int(n[1]))
            if c1 < 0 or c2 > 0:
                raise ValueError('illegal ad_cut: %s' % s)
        if cut_para:
            c_para = '--cut %s --cut %s' % (c1, c2)
        else:
            c_para = [c1, c2]
    else:
        if cut_para:
            c_para = '--cut %s' % int(s)
        else:
            c_para = [int(s), ]
    return c_para


def ends_trimmer(fn, path_out, cut_after_trim='0', len_min=15):
    """
    trim N-bases at either end of the read
    """
    assert os.path.exists(fn)
    assert isinstance(len_min, int)
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    fn_out_name = file_prefix(fn)[0] + '.cut.fastq'
    fn_out_file = os.path.join(path_out, fn_out_name)
    # trim either ends
    tmp = cutadapt_cut(cut_after_trim, False)
    if len(tmp) == 2:
        trim_5, trim_3 = tmp
    elif len(tmp) == 1:
        trim_5 = 0 if(int(tmp[0]) < 0) else tmp[0]
        trim_3 = tmp[0] if(int(tmp[0]) < 0) else 0
    else:
        trim_5 = trim_3 = 0
    with xopen(fn, 'rt') as fi, open(fn_out_file, 'wt') as fo:
        while True:
            try:
                fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(), 
                                                   next(fi).strip(), 
                                                   next(fi).strip(), 
                                                   next(fi).strip(),]
                if len(fq_seq) < len_min + abs(trim_5) + abs(trim_3): 
                    # skip short reads
                    continue
                fq_seq = fq_seq[trim_5:trim_3] if(trim_3 < 0) else fq_seq[trim_5:]
                fq_qual = fq_qual[trim_5:trim_3] if(trim_3 < 0) else fq_qual[trim_5:]
                fo.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
            except StopIteration:
                break
    return fn_out_file



def se_trimmer(fn, path_out=None, len_min=15, adapter3=None, double_trim=True, 
               qual_min=20, err_rate=0.1, multi_cores=1, rm_untrim=False,
               overlap=4, pcr_r2=None, truseq_p7a=None, truseq_p7b=None, 
               truseq_p7c=None, truseq_univ=None, adapter5=None,
               cut_before_trim=0, overwrite=False):
    """
    using cutadapt to trim SE read
    support one input only
    """
    assert os.path.exists(fn)
    assert isinstance(adapter3, str)
    path_out = os.path.dirname(read_in) if path_out is None else path_out
    assert is_path(path_out)
    fn_out_file = os.path.join(path_out, file_prefix(fn)[0] + '.clean.fastq')
    log_out_file = os.path.join(path_out, file_prefix(fn)[0] + '.cutadapt.log')
    # sliding adapter3
    ads = adapter_chopper(adapter3)
    for i in [pcr_r2, truseq_p7a, truseq_p7b, truseq_p7c, truseq_univ]:
        if i is None:
            continue
        else:
            ads.append(i)
    para_ad = ' '.join(['-a {}'.format(i) for i in ads])
    # for adapter5
    if isinstance(adapter5, str):
        para_ad += ' -g %s' % adapter5
    # for untrim
    if rm_untrim is False:
        fn_untrim_file = os.path.join(path_out, 
                                      file_prefix(fn)[0] + '.untrim.fastq')
        para_adx = ' --untrimmed-output=%s --cores=%s' % (fn_untrim_file, 1)
        # untrim, not support multiple threads
    else:
        para_adx = ' --cores=%s' % multi_cores
    # cut adapter *before* trimming
    para_adx  += ' ' + cutadapt_cut(cut_before_trim)
    ## command line
    c1 = 'cutadapt %s %s -m %s -q %s --overlap=%s --error-rate=%s --times=4 \
          --trim-n --max-n=0.1 %s' % (para_ad, para_adx, len_min, qual_min, 
                                      overlap, err_rate, fn)
    c2 = 'cutadapt %s -m %s -q %s --overlap=%s --error-rate=%s --times=4 \
          --cores=%s --trim-n --max-n=0.1 -' % (para_ad, len_min, 
                                                      qual_min, overlap, 
                                                      err_rate, multi_cores)
    if not os.path.isfile(fn_out_file) or overwrite is True:
        with open(log_out_file, 'w') as fo1, open(log_out_file, 'a') as fo2, \
             open(fn_out_file, 'wt') as fr:
            p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE, 
                                  stderr=fo1)
            p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=fr,
                                  stderr=fo2)
            p3 = p2.communicate()
        tmp1 = cutadapt_log_parser(log_out_file) # processing log
    return fn_out_file



def pe_trimmer():
    """
    using cutadapt to trim PE reads
    """
    pass


def trim(fns, path_out, adapter3=None, len_min=15, qual_min=20,
         err_rate=0.1, overlap=1, multi_cores=1, read12=1,
         rm_untrim=False, pcr_r2=None, truseq_p7a=None, 
         truseq_univ=None, adapter5=None, rm_dup=False, 
         cut_before_trim='0', cut_after_trim='0', overwrite=False):
    """
    processing reads
    """
    # default Truseq adapters
    ad3 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    ad5 = 'GATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    ad5rc = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATC'
    ad_p7a = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    fn_out_files = []
    for fn in fns:
        assert os.path.exists(fn)
        logging.info('trimming reads: %s' % file_prefix(fn)[0])
        if rm_dup is True:
            if not cut_after_trim == '0':
                fn_name = file_prefix(fn)[0] + '.clean.nodup.cut.fastq'
            else:
                fn_name = file_prefix(fn)[0] + '.clean.nodup.fastq'
        elif not cut_after_trim == '0':
            fn_name = file_prefix(fn)[0] + '.clean.cut.fastq'
        else:
            fn_name = file_prefix(fn)[0] + '.clean.fastq'
        fn_out_check = os.path.join(path_out, fn_name)
        # file existence
        if os.path.exists(fn_out_check) and not overwrite:
            fn_out_files.append(fn_out_check)
            logging.info('file exists: %s' % fn_name)
            continue # next
        if read12 == 1:
            # read1 of PE reads
            adapter3 = ad3 if adapter3 is None else adapter3
            truseq_p7a = ad_p7a if truseq_p7a is None else truseq_p7a
            fn_out_file = se_trimmer(fn, path_out, len_min, adapter3, 
                                     qual_min=qual_min, 
                                     err_rate=err_rate, 
                                     multi_cores=multi_cores, 
                                     overlap=overlap, 
                                     truseq_p7a=truseq_p7a,
                                     rm_untrim=rm_untrim,
                                     cut_before_trim=cut_before_trim,
                                     overwrite=overwrite)
        elif read12 == 2:
            # read2 or pE reads
            adapter3 = ad5rc if adapter3 is None else adapter3
            fn_out_file = se_trimmer(fn, path_out, len_min, adapter3,
                                     qual_min=qual_min,
                                     err_rate=err_rate, 
                                     multi_cores=multi_cores, 
                                     overlap=overlap,
                                     rm_untrim=rm_untrim,
                                     cut_before_trim=cut_before_trim,
                                     overwrite=overwrite)
        else:
            logging.error('unknown --read12: %s' % read12)
        ## post-processing
        if rm_dup is True:
            if not cut_after_trim == '0':
                f1 = dup_remover(fn_out_file, path_out, q=qual_min)
                f2 = ends_trimmer(f1, path_out, cut_after_trim=cut_after_trim, len_min=len_min)
                os.remove(fn_out_file)                
                os.remove(f1)
                fn_out_file = f2
            else:
                f3 = dup_remover(fn_out_file, path_out, q=qual_min)
                os.remove(fn_out_file)
                fn_out_file = f3
        elif not cut_after_trim == '0':
            f4 = ends_trimmer(fn_out_file, path_out, cut_after_trim=cut_after_trim, len_min=len_min)
            os.remove(fn_out_file)
            fn_out_file = f4
        else:
            pass # return fn_out_file
        ## post-stat
        fn_out_files.append(fn_out_file)
        fn_n = int(file_row_counter(fn_out_file) / 4)
        fn_n_file = os.path.join(path_out, file_prefix(fn)[0] + '.reads.txt')
        with open(fn_n_file, "wt") as f:
                f.write(str(fn_n) + '\n')
    return fn_out_files



def main():
    """
    Processing CLIP reads
    """
    args = get_args()
    fq_files = [f.name for f in args.i]
    path_out = args.o

    p = trim(fq_files, path_out, adapter3=args.a, len_min=args.m, 
               qual_min=args.q, err_rate=args.e, multi_cores=args.threads,
               rm_untrim=args.rm_untrim, overlap=args.O, 
               rm_dup=args.rm_dup, cut_before_trim=args.cut_before_trim,
               cut_after_trim=args.cut_after_trim, overwrite=args.overwrite)
    # p = trim(fqs = [f.name for f in args.i], adapter3 = args.a, 
    #     out_path = args.o, len_min = args.m, qual_min = args.q, cut = args.cut,
    #     err_rate = args.e, multi_cores = args.threads, overlap = args.O)
    

if __name__ ==  '__main__':
    main()


## EOF
