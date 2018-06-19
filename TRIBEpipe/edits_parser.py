#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract edit sites from BAM file

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
import pysam


logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'edits_parser', description = 'Parsing TRIBE editing sites',
        epilog = 'Output: ')
    parser.add_argument('-t', default = 'RNA', choices = ['RNA', 'DNA'],
        help = 'SNP type (RNA|DNA) for input BAM file, default: RNA')
    parser.add_argument('-i', required = True, metavar = 'BAM', 
        type = argparse.FileType('r'),
        help = 'alignments in BAM format')
    parser.add_argument('-g', required = True, metavar = 'Genome',
        type = argparse.FileType('r'),
        help = 'Reference sequence in FASTA format')
    parser.add_argument('-o', required = True, metavar = 'OUTPUT',
        help = 'file to save the results')    
    parser.add_argument('--depth_cutoff', default = 2, type = int, 
        metavar = 'depth',
        help = 'minimum read depth at editing position, default: 1')
    parser.add_argument('--pct_cutoff', default = 10, type = int,
        metavar = 'percentage',
        help = 'minimum editing percentage [1-100]%%, default: 10')
    parser.add_argument('--overwrite', action = 'store_true',
        help = 'if specified, overwrite existing file')
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


def rna_snp_filter(fs, depth_cutoff, pct_cutoff, strand = '+'):
    """
    filt bases: count>depth_cutoff, pct>pct_cutoff
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
        if not refbase == "T" or freqC < pct_cutoff or int(numTotal) < depth_cutoff:
            return None
    else:
        freqHit = freqG
        if not refbase == "A" or freqG < pct_cutoff or int(numTotal) < depth_cutoff:
            return None
    # BED output
    start = int(fs[1]) - 1
    name = fs[0] + '_' + fs[1] + '_{}_{}%'.format(numTotal, freqHit)
    f_out = [fs[0], start, fs[1], freqHit, name, fs[2]] + fs[3:8] + [str(numN)]
    f_out = list(map(str, f_out))
    return f_out


def dna_snp_filter(fs, depth_cutoff, pct_cutoff, strand = '+'):
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
        if not refbase == "T" or freqT < pct_cutoff or freqC > 0:
            return None
    else:
        freqHit = freqA
        if not refbase == "A" or freqA < pct_cutoff or freqG > 0:
            return None
    # BED output
    start = int(fs[1]) - 1
    name = fs[0] + '_' + fs[1] + '_{}_{}%'.format(numTotal, freqHit)
    f_out = [fs[0], start, fs[1], freqHit, name, fs[2]] + fs[3:8] + [str(numN)]
    f_out = list(map(str, f_out))
    return f_out


def snp_parser(bam, genome, out_file, strand = '+', snp_type = 'rna', 
               depth_cutoff = 1, pct_cutoff = 10, append = False):
    """
    extract nucleotide count table at each position
    """
    # parameters
    # out_path
    out_path = os.path.dirname(out_file)
    if out_path != '' and not os.path.exists(out_path):
        os.makedirs(out_path)
    # indexed genome
    if not os.path.exists(genome + '.fai'):
        tmp = pysam.faidx(genome) # make faidx
    # dna or rna
    if snp_type.lower() == 'dna':
        snp_filter = dna_snp_filter
    elif snp_type.lower() == 'rna':
        snp_filter = rna_snp_filter
    else:
        return None
    # strand
    if strand == '+':
        c1 = 'samtools view -f 16 -bhS {}'.format(bam)
    elif strand == '-':
        c1 = 'samtools view -F 16 -bhS {}'.format(bam)
    else:
        return None
    # values
    assert isinstance(depth_cutoff, int) and depth_cutoff > 0
    assert isinstance(pct_cutoff, int) and pct_cutoff >= 0 and pct_cutoff <= 100

    # write mode
    write_mode = 'at' if append else 'wt'

    # run commands
    c1 = 'samtools view -f 16 -bhS {}'.format(bam)
    c2 = 'samtools mpileup -f {} -'.format(genome)
    p1 = subprocess.Popen(shlex.split(c1), stdout = subprocess.PIPE)
    p2 = subprocess.Popen(shlex.split(c2), stdin = p1.stdout, 
                          stdout = subprocess.PIPE, stderr = subprocess.PIPE, 
                          universal_newlines = True)
    
    # parsing output
    with open(out_file, write_mode) as fo:
        while True:
            line = p2.stdout.readline().strip()
            if not line: break
            chr, pos, ref, depth, pileup, quality = line.strip().split('\t')
            ref = ref.upper()
            depth = int(depth)
            if depth >= depth_cutoff and depth_cutoff > 0 and ref != 'N':
                acgt_res = mpileup2acgt(pileup, quality, depth, ref, qlimit = 25, 
                                        noend = False, nostart = False)
                f_out = [chr, pos, ref, depth, acgt_res['A'], acgt_res['C'], 
                         acgt_res['G'], acgt_res['T']]
                # filtering
                f_out2 = snp_filter(f_out, depth_cutoff, pct_cutoff, 
                                    strand = strand)
                if not f_out2 is None:
                    fo.write('\t'.join(list(map(str, f_out2))) + '\n')
    return True



def edits_parser(bam, genome, outfile, snp_type, depth_cutoff, pct_cutoff,
                 overwrite = False):
    if os.path.exists(outfile) and overwrite is False:
        logging.info('edits file exists, skipping: ' + outfile)
    else:
        try:
            # forward strand    
            snp_parser(bam, genome, outfile, strand = '+', snp_type = snp_type,
                       depth_cutoff = depth_cutoff, pct_cutoff = pct_cutoff)
            # reverse strand
            snp_parser(bam, genome, outfile, strand = '-', snp_type = snp_type,
                       depth_cutoff = depth_cutoff, pct_cutoff = pct_cutoff, 
                       append = True)
        except IOError:
            logging.info('fail to call edits')
    return outfile




def main():
    logging.info('calling edits')
    args = get_args() #
    assert args.depth_cutoff > 0
    assert args.pct_cutoff >= 0 and args.pct_cutoff <= 100
    edits_parser(args.i.name, args.g.name, args.o, args.t, args.depth_cutoff,
                 args.pct_cutoff, overwrite = args.overwrite)
    # # forward strand
    # snp_parser(args.i.name, args.g.name, args.o, strand = '+', 
    #            snp_type = args.t, depth_cutoff = args.depth_cutoff,
    #            pct_cutoff = args.pct_cutoff)
    # # reverse strand
    # snp_parser(args.i.name, args.g.name, args.o, strand = '-', 
    #            snp_type = args.t, depth_cutoff = args.depth_cutoff,
    #            pct_cutoff = args.pct_cutoff, append = True)

if __name__ == '__main__':
    main()


## EOF
