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
"""

__author__ = "Ming Wang <wangm08@hotmail.com>"
__copyright__ = "2018 by Ming Wang <wangm08@hotmail.com>"
__license__ = "MIT"
__email__ = "wangm08@hotmail.com"
__version__ = "0.1"

import os, sys
import re, datetime, json, itertools
import argparse, subprocess, shlex
import binascii
import logging
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
    parser.add_argument('-m', default = 15, metavar = 'len_min', 
        type = int, help = 'Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('-q', default = 20, metavar = 'quality', 
        type = int, help = 'The cutoff of base quality, default [20]')
    parser.add_argument('-p', default = 80, metavar = 'PERCENT', type = int,
        help = 'minimum percent of bases that must have -q quality, default [80]')
    parser.add_argument('--cut', default = 0, metavar = 'CUT', type = int,
        help = 'cut N bases from each read before trimming, plus: at 5 end, \
                minus at 3 end default [0]')
    parser.add_argument('-e', default = 0.1, metavar = 'err_rate', type = float,
        help = 'Maximum allowed error rate, default [0.1]')
    parser.add_argument('--threads', default = 1,  metavar = 'THREAD', type = int,
        help = 'Number of threads to launch, default [1]')
    parser.add_argument('-O', default = 1, metavar = 'OVERLAP',
        help = 'Required N bases overlap between reads and adapter, default [1]')
    args = parser.parse_args()
    return args


## Functions
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def file_row_counter(file):
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

    f_reader = gzip.open if(is_gz_file(file)) else open
    with f_reader(file, "rt", encoding="utf-8", errors='ignore') as f:
        return (sum(bl.count("\n") for bl in blocks(f)))



def se_trimmer(read_in, adapter3, len_min = 15, *, out_path = None, 
    qual_min = 20, err_rate = 0.2, multi_cores = 1, overlap = 4, 
    cut = 0):
    """
    using cutadapt to trim SE read
    support one input only
    """
    pcr_r2 = 'ATCTCGTATGCCGTCTTCTGCTTG'
    truseq_p7a = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    truseq_p7b = 'GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG'
    truseq_p7c = 'GCGGTTCAGCAGGAATGCCGAG'
    truseq_univ = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    if adapter3 is None:
        ads = [pcr_r2, truseq_p7a, truseq_p7b, truseq_p7c, truseq_univ]
    else:
        ads = [adapter3, pcr_r2, truseq_p7a, truseq_p7b, truseq_p7c, truseq_univ]
    para_ad = ' '.join(['-a {}'.format(i) for i in ads])
    ## cut N
    if not cut == 0:
        para_ad += ' --cut={}'.format(cut)
    ## para checker
    if out_path is None:
        out_path = os.path.dirname(read_in) # default with read_in
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    ## filenames
    read_out_name = re.sub(r'.f[ast]*q(.gz)?', '.clean.fq', os.path.basename(read_in))
    read_out = os.path.join(out_path, read_out_name)
    log_out = read_out.replace('.clean.fq', '.cutadapt.log')
    ## commands
    c1 = 'cutadapt {} -m {} -q {} --overlap={} \
          --error-rate={} --times=4 --cores={} --trim-n --max-n=0.1 {}'\
          .format(para_ad, len_min, qual_min, overlap, err_rate, 
                  multi_cores, read_in)
    if not os.path.isfile(read_out):
        with open(log_out, 'w') as fo, open(read_out, 'w') as fr:
            p1 = subprocess.run(shlex.split(c1), stdout = fr, stderr = fo)
            tmp1 = cutadapt_log_parser(log_out)
    return read_out


def pe_trimmer():
    """
    using cutadapt to trim PE reads
    """
    pass


def cutadapt_log_parser(log):
    """
    Parsing the output of cutadpat
    support multiple run output
    save data in JSON format
    """
    logdict = {}
    _cutadapt_ver = []
    _python_ver = []
    _cmd = []
    _fname = []
    _raw = []
    _clean = []
    outfile = os.path.splitext(log)[0] + ".json"
    with open(log, 'r') as f:
        for line in f.readlines():
            if(len(re.findall('^This is cutadapt', line)) == 1):
                r1 = re.findall(r'(cutadapt|Python) (\d+\.\d+\.?\d+?)', line) # 1st line, tools
                for i in r1: logdict[i[0]] = i[1]
            elif('Command line parameters' in line):
                _cmd.append(re.sub(r'Command line parameters: ', '', line).strip('\n'))
                _fname.append(os.path.basename(line.split()[-1])) # file name
            elif('Total reads processed:' in line):
                 _raw.append(re.sub(r',', '', line.split()[-1])) # first input
            elif('Reads written (passing filters):' in line):
                 _clean.append(re.sub(r',', '', line.split()[-2])) # output
            else:
                 continue
    logdict['filename'] = _fname[0]
    logdict['run_cutadapt_times'] = len(_raw)
    logdict['command_lines'] = _cmd
    logdict['raw'] = _raw[0]
    logdict['clean'] = _clean[-1]
    logdict['clean_pct'] = '{:.2f}%'.format(int(logdict['clean'])/int(logdict['raw'])*100)
    with open(outfile, 'w') as fo:
        json.dump(logdict, fo, indent = 4)
    return outfile



def trim(fqs, adapter3, out_path, len_min = 15, qual_min = 20, cut = 0, 
    multi_cores = 1, err_rate = 0.1, overlap = 1, read12 = 1, overwrite = False):
    """
    processing reads
    trim - nodup
    """
    logging.info('trimming reads')
    fq_out = []
    for fq in fqs:
        q_prefix = re.sub(r'.f[ast]*q(.gz)?', '', os.path.basename(fq))
        q_out = os.path.join(out_path, q_prefix + '.clean')
        if abs(cut) > 0:
            q_out += '.cut' + str(cut) + '.fq'
        else:
            q_out += '.fq'
        if os.path.isfile(q_out) and not overwrite:
            fq_out.append(q_out)
            logging.info('file exists, trimming skipped...' + q_out)
        else:
            q_out = se_trimmer(fq, adapter3, len_min, out_path = out_path,
                               qual_min = qual_min, err_rate = err_rate,
                               multi_cores = multi_cores, overlap = overlap, 
                               cut = cut)
            # count reads
            fq_n = int(file_row_counter(q_out) / 4)
            fq_prefix = re.sub(r'.f[ast]*q(.gz)?', '', os.path.basename(fq))
            fq_n_file = os.path.join(out_path, fq_prefix + ".clean.txt")
            with open(fq_n_file, "w") as f:
                f.write(str(fq_n) + '\n')
            fq_out.append(q_out)
    return fq_out


def main():
    """
    Processing CLIP reads
    """
    args = get_args()
    p = trim(fqs = [f.name for f in args.i], adapter3 = args.a, 
        out_path = args.o, len_min = args.m, qual_min = args.q, cut = args.cut,
        err_rate = args.e, multi_cores = args.threads, overlap = args.O)
    

if __name__ ==  '__main__':
    main()


## EOF
