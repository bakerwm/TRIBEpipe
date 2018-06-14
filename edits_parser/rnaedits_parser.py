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
"""

import os, sys, argparse, tempfile
import shlex, subprocess


def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'edits_parser', description = 'Parsing editing events',
        epilog = 'Example: ')
    parser.add_argument('-i', required = True, metavar = 'BAM', 
        type = argparse.FileType('r'),
        help = 'BAM file')
    parser.add_argument('-g', required = True, metavar = 'Genome',
        type = argparse.FileType('r'),
        help = 'Reference sequence in FASTA format, indexed')
    parser.add_argument('-o', required = True, metavar = 'OUTPUT',
        help = 'path to save the results')
    parser.add_argument('-t', default = 'RNA', choices = ['RNA', 'DNA'],
        help = 'SNP type (RNA|DNA) for input BAM file')
    parser.add_argument('--min_depth', default = 2, type = int, 
        metavar = 'depth',
        help = 'minimum read depth at editing position, default: 1')
    parser.add_argument('--min_pct', default = 10, type = int,
        metavar = 'percentage',
        help = 'minimum editing percentage (1-100), default: 10')
    args = parser.parse_args()
    return args


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
    fout = [fs[0], str(start), fs[1], name, '200', strand, fs[2]] + fs[3:]
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
    fout = [fs[0], str(start), fs[1], name, '200', strand, fs[2]] + fs[3:]
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
        px = subprocess.run(shlex.split(cs))
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
