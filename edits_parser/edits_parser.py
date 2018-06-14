#!/usr/bin/env python
"""
extract edit sites from BAM

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

1. samtools view: separate forward strand (-f 16) and reverse strand (-F 16)
2. samtools mpileup: filt and extract nucleotide content at each position
3. sequenza-utils pileup2acgt: extract ACGT content and frequency
4. functions: further filt editsites


Compute each sample respectively:
TRIBE - gDNA - wt_mRNA = sites


output:
1       chrom
2       chromStart
3       chromEnd
4       percentage_of_edits%
5       name of edits, [chrom]_[chromEnd]_[read_depth]_[percentage_of_edit]%
6       refBase
7       read_depth
8       A count
9       C count
10      G count
11      T count
12      N count       
"""

import os, sys, argparse, tempfile, logging
import shlex, subprocess

import logging
logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'edits_parser', description = 'Parsing editing events',
        epilog = 'Output: ')
    parser.add_argument('-t', default = 'RNA', choices = ['RNA', 'DNA'],
        help = 'SNP type (RNA|DNA) for input BAM file')
    parser.add_argument('-i', required = True, metavar = 'BAM', 
        type = argparse.FileType('r'),
        help = 'BAM file')
    parser.add_argument('-g', required = True, metavar = 'Genome',
        type = argparse.FileType('r'),
        help = 'Reference sequence in FASTA format, indexed')
    parser.add_argument('-o', required = True, metavar = 'OUTPUT',
        help = 'path to save the results')    
    parser.add_argument('--min_depth', default = 2, type = int, 
        metavar = 'depth',
        help = 'minimum read depth at editing position, default: 1')
    parser.add_argument('--min_pct', default = 10, type = int,
        metavar = 'percentage',
        help = 'minimum editing percentage [1-100]%%, default: 10')
    args = parser.parse_args()
    return args



def mpileup2acgt(pileup, quality, depth, reference, qlimit=53,
         noend=False, nostart=False):
    """
    This function was adopted from: sequenza-utils pileup2acgt 
    URL: https://bitbucket.org/sequenza_tools/sequenza-utils
    original code were protected under GPLv3 license.

    Parse the mpileup format and return the occurrence of
    each nucleotides in the given positions.

    pileup format:

    1    chr
    2    1-based coordinate
    3    reference base
    4    depth
    5    base qualities
    6    mapping qualities
    """
    nucleot_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    strand_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    n = 0
    block = {'seq': '', 'length': 0}
    start = False
    del_ins = False
    l_del_ins = ''
    last_base = None
    ins_del_length = 0
    for base in pileup:
        if block['length'] == 0:
            if base == '$':
                if noend:
                    if last_base:
                        nucleot_dict[last_base.upper()] -= 1
                        if last_base.isupper():
                            strand_dict[last_base.upper()] -= 1
                    last_base = None
            elif base == '^':
                start = True
                block['length'] += 1
                block['seq'] = base
            elif base == '+' or base == '-':
                del_ins = True
                block['length'] += 1
                block['seq'] = base
            elif base == '.' or base == ',':
                if ord(quality[n]) >= qlimit:
                    nucleot_dict[reference] += 1
                    if base == '.':
                        strand_dict[reference] += 1
                        last_base = reference
                    else:
                        last_base = reference.lower()
                else:
                    last_base = None
                n += 1
            elif base.upper() in nucleot_dict:
                if ord(quality[n]) >= qlimit:
                    nucleot_dict[base.upper()] += 1
                    if base.isupper():
                        strand_dict[base.upper()] += 1
                    last_base = base
                else:
                    last_base = None
                n += 1
            else:
                n += 1
        else:
            if start:
                block['length'] += 1
                block['seq'] += base
                if block['length'] == 3:
                    if not nostart:
                        if base == '.' or base == ',':
                            if ord(quality[n]) >= qlimit:
                                nucleot_dict[reference] += 1
                                if base == '.':
                                    strand_dict[reference] += 1
                        elif base.upper() in nucleot_dict:
                            if ord(quality[n]) >= qlimit:
                                nucleot_dict[base.upper()] += 1
                                if base.isupper():
                                    strand_dict[base.upper()] += 1
                    block['length'] = 0
                    block['seq'] = ''
                    start = False
                    last_base = None
                    n += 1
            elif del_ins:
                if base.isdigit():
                    l_del_ins += base
                    block['seq'] += base
                    block['length'] += 1
                else:
                    ins_del_length = int(l_del_ins) + 1 + len(l_del_ins)
                    block['seq'] += base
                    block['length'] += 1
                    if block['length'] == ins_del_length:
                        block['length'] = 0
                        block['seq'] = ''
                        l_del_ins = ''
                        # ins_del = False
                        ins_del_length = 0

    nucleot_dict['Z'] = [strand_dict['A'], strand_dict[
        'C'], strand_dict['G'], strand_dict['T']]
    return nucleot_dict




def rna_snp_filter(fs, ncount, npct, strand = '+'):
    """
    filt bases: count>ncount, pct>npct
    fit: samtools pileup output
    chr n_base ref_base read.depth A C G T strand

    example:

    column  header
    1       chr
    2       coord
    3       refbase
    4       total
    5       A
    6       C
    7       G
    8       T
    9       strand
    """
    refbase, numTotal, numA, numC, numG, numT = fs[2:8]
    numN = int(numTotal) - int(numA) - int(numC) - int(numG) - int(numT)
    freqA = int(float(numA) / float(numTotal) * 100)
    freqC = int(float(numC) / float(numTotal) * 100)
    freqG = int(float(numG) / float(numTotal) * 100)
    freqT = int(float(numT) / float(numTotal) * 100)
    freqHit = 0
    if strand == "-":
        freqHit = freqC
        if not refbase == "T" or freqC < npct or int(numTotal) < ncount:
            return None
    else:
        freqHit = freqG
        if not refbase == "A" or freqG < npct or int(numTotal) < ncount:
            return None
    # BED output
    start = int(fs[1]) - 1
    name = fs[0] + '_' + fs[1] + '_{}_{}%'.format(numTotal, freqHit)
    # fout = [fs[0], str(start), fs[1], name, '200', strand, fs[2]] + fs[3:]
    fout = [fs[0], str(start), fs[1], str(freqHit), name, fs[2]] + fs[3:8] + [str(numN)]
    return fout


def dna_snp_filter(fs, ncount, npct, strand = '+'):
    """
    filt bases: freqA > 80%, freqG = 0
    fit: samtools pileup output
    chr n_base ref_base read.depth A C G T strand

    example:

    column  header
    1       chr
    2       coord
    3       refbase
    4       total
    5       A
    6       C
    7       G
    8       T
    9       strand
    """
    refbase, numTotal, numA, numC, numG, numT = fs[2:8]
    numN = int(numTotal) - int(numA) - int(numC) - int(numG) - int(numT)
    freqA = int(float(numA) / float(numTotal) * 100)
    freqC = int(float(numC) / float(numTotal) * 100)
    freqG = int(float(numG) / float(numTotal) * 100)
    freqT = int(float(numT) / float(numTotal) * 100)
    freqHit = 0
    if strand == "-":
        freqHit = freqT
        if not refbase == "T" or freqT < npct or freqC > 0:
            return None
    else:
        freqHit = freqA
        if not refbase == "A" or freqA < npct or freqG > 0:
            return None
    # BED output
    start = int(fs[1]) - 1
    name = fs[0] + '_' + fs[1] + '_{}_{}%'.format(numTotal, freqHit)
    # fout = [fs[0], str(start), fs[1], name, '200', strand, fs[2]] + fs[3:]
    fout = [fs[0], str(start), fs[1], str(freqHit), name, fs[2]] + fs[3:8] + [str(numN)]
    return fout


def snp_parser(fn, gfa, out_file, strand = '+', snp_type = 'rna',
               min_depth = 1, npct = 10):
    """extract nucleotide count table at each position"""
    # parameters
    assert isinstance(min_depth, int)
    out_path = os.path.dirname(out_file)
    if (os.path.exists(out_path) or out_path == ''):
        pass
    else:
        os.makedirs(out_path)
    # faidx indexed file
    gfa_idx = gfa + '.fai'
    if not os.path.exists(gfa_idx):
        cx = 'samtools faidx {}'.format(gfa)
        px = subprocess.run(shlex.split(cx))
    # dna or rna
    if snp_type.lower() == 'dna':
        snp_filter = dna_snp_filter
    elif snp_type.lower() == 'rna':
        snp_filter = rna_snp_filter
    else:
        return None
    # check strand
    if strand == '+':
        c1 = 'samtools view -f 16 -bhS {}'.format(fn)
    elif strand == '-':
        c1 = 'samtools view -F 16 -bhS {}'.format(fn)
    else:
        return None
    # devnull
    fnull = open(os.devnull, 'w')
    c2 = 'samtools mpileup -f {} -'.format(gfa)
    c3 = 'sequenza-utils pileup2acgt -n {} -p -'.format(min_depth)
    p1 = subprocess.Popen(shlex.split(c1), stdout = subprocess.PIPE)
    p2 = subprocess.Popen(shlex.split(c2), stdin = p1.stdout, 
                          stdout = subprocess.PIPE, stderr = fnull)
    p3 = subprocess.Popen(shlex.split(c3), stdin = p2.stdout, 
                           stdout = subprocess.PIPE,
                           universal_newlines = True)
    # filtering
    with open(out_file, 'w') as fo:
        while True:
            line = p3.stdout.readline().strip()
            if not line:
                break
            if 'ref_base' in line:
                continue
            tabs = line.split('\t')
            tabs_new = snp_filter(tabs, min_depth, npct, strand = strand)
            if not tabs_new is None:
                fo.write('\t'.join(tabs_new) + '\n')
    fo.close() # 
    return True


def main():
    logging.info('calling edits')
    args = get_args() #
    assert args.min_depth > 0
    assert args.min_pct >= 0 and args.min_pct <= 100
    # # for DNA SNPs
    # if args.t.lower() == 'dna':
    #     args.min_pct = 80 # predefined
    # forward
    file_fwd = args.o + '.fwd.tmp'
    snp_parser(args.i.name, args.g.name, file_fwd, strand = '+', 
               snp_type = args.t, min_depth = args.min_depth, 
               npct = args.min_pct)
    # reverse
    file_rev = args.o + '.rev.tmp'
    snp_parser(args.i.name, args.g.name, file_rev, strand = '-', 
               snp_type = args.t, min_depth = args.min_depth, 
               npct = args.min_pct)
    # merge two files
    os.system('cat {} {} > {}'.format(file_fwd, file_rev, args.o))
    os.remove(file_fwd)
    os.remove(file_rev)


if __name__ == '__main__':
    main()


## EOF
