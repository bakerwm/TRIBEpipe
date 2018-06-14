#!/usr/bin/env python
"""
parsing samtools mpileup output

"""


def acgt(pileup, quality, depth, reference, qlimit=53,
         noend=False, nostart=False):
    '''
    Parse the mpileup format and return the occurrence of
    each nucleotides in the given positions.

    pileup format:

    1    chr
    2    1-based coordinate
    3    reference base
    4    depth
    5    base qualities
    6    mapping qualities
    '''
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


def pileup2acgt():
    fn = 'out.mpileup'
    outfile = 'out.matrix'
    depth_cutoff = 2
    quality_cutoff = 25
    nostart = False
    noend = False
    header = ['chr', 'position', 'ref_base', 'read.depth', 'A', 'C', 'G', 'T']
    with open(fn, 'rt') as fi, open(outfile, 'wt') as fo:
        fo.write('\t'.join(header) + '\n')
        for line in fi:
            try:
                chromosome, position, reference, \
                    depth, pileup, quality = line.strip().split('\t')
                depth = int(depth)
                reference = reference.upper()
                if depth >= depth_cutoff and depth_cutoff > 0 and reference != 'N':
                    acgt_res = acgt(pileup, quality, depth,
                                    reference, qlimit = quality_cutoff,
                                    noend=noend, nostart=nostart)
                    fo.write('\t'.join(map(str, [chromosome, position, reference, depth, 
                        acgt_res['A'], acgt_res['C'], acgt_res['G'], acgt_res['T']])) + '\n')
            except ValueError:
                pass

pileup2acgt()