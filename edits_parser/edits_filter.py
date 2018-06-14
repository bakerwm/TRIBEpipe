#!/usr/bin/env python
"""filt edits by threshold, >= 20 reads"""

import os, sys

usage = """
Usage: python threshold_filter.py <N|20> <pct%|10> [<file1>, <file2>, ...]

Criteria:
1. >= 20 reads (column 9)
2. >= 10% edits pct (comumn 15 / column 16) editRNA / totalRNA

conversion: 
forward strand: A -> G
reverse strand: T -> C

"""

def base_filter3(fs, ncount, npct):
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
    refbase, numTotal, numA, numC, numG, numT, strand = fs[2:]
    freqA = int(float(numA) / float(numTotal) * 100)
    freqC = int(float(numC) / float(numTotal) * 100)
    freqG = int(float(numG) / float(numTotal) * 100)
    freqT = int(float(numT) / float(numTotal) * 100)
    freqHit = 0
    if strand == "-":
        freqHit = freqC
        if not refbase == "T" or freqC < npct or int(numTotal) >= ncount or freqC > 90:
            return None
    else:
        freqHit = freqG
        if not refbase == "A" or freqG < npct or int(numTotal) >= ncount or freqG > 90:
            return None
    # BED output
    start = int(fs[1]) - 1
    name = fs[0] + '_' + fs[1] + '_{}_{}%'.format(numTotal, freqHit)
    fout = [fs[0], str(start), fs[1], name, '200', fs[2]] + fs[3:]
    return fout


def base_filter(fs, ncount, npct):
    """
    fit: HyperTRIBE sam2matrix output
    filt bases: count>ncount, pct>npct

    example:

    column  header
    1       chr
    2       coord
    3       gene
    4       type
    5       A
    6       T
    7       C
    8       G
    9       Total
    10      AgDNA
    11      TgDNA
    12      CgDNA
    13      GgDNA
    14      TotalgDNA
    15      editRNAcount *
    16      totalRNAcount *
    17      editgDNAcount
    18      totalgDNAcount
    """
    numRNAEdit, numRNATotal = fs[14:16]
    freqEdit = int(float(nbase) / float(ntotal) * 100)
    if int(numRNATotal) == 0 or int(numRNATotal) < ncount or freqEdit < npct:
        return None
    else:
        end = int(fs[1]) + 1
        name = '{:d}%_{}r'.format(freqEdit, numRNATotal)
        tag = fs[0] + '_' + str(fs[1])
        fout = fs[0:2] + [end, str(freqEdit)] + fs + tag
        return fout


def base_filter2(fs, ncount, npct):
    """
    fit: Rsamtools pileup output
    filt bases: count>ncount, pct>npct
        
    example:

    column  header
    1       chr
    2       coord
    3       strand
    4       total
    5       A
    6       C
    7       G
    8       T

    """
    strand, numTotal, numA, numC, numG, numT = fs[2:]
    freqA = int(float(numA) / float(numTotal) * 100)
    freqC = int(float(numC) / float(numTotal) * 100)
    freqG = int(float(numG) / float(numTotal) * 100)
    freqT = int(float(numT) / float(numTotal) * 100)
    freqHit = 0
    if strand == "-":
        freqHit = freqC
        if int(numT) < 1 or freqC < npct or int(numTotal) < ncount or freqC > 90:
            return None
    else:
        freqHit = freqG
        if int(numA) < 1 or freqG < npct or int(numTotal) < ncount or freqG > 90:
            return None
    # BED output
    start = int(fs[1]) - 1
    name = fs[0] + '_' + fs[1] + '_{}_{}%'.format(numTotal, freqHit)
    fout = [fs[0], str(start), fs[1], name, '200', fs[2]] + fs[3:]
    return fout


# main
if len(sys.argv) < 4:
    print("parameters not enough")
    print(usage)
    sys.exit(1)


ncount = int(sys.argv[1])
npct = int(sys.argv[2])
mfiles = sys.argv[3:]

for f in mfiles:
    outfile = os.path.splitext(f)[0] + '.cutoff_' + str(ncount) + '.freq_' + str(npct) + '.bed'
    with open(f) as fi, open(outfile, 'w') as fo:
        for line in fi:
            line = line.strip()
            if 'Edit_coord' in line: # remove header
                # fo.write(line + '\n')
                continue
            elif 'read.depth' in line: # remove header
                continue
            else:
                tabs = line.split('\t')
                p = base_filter3(tabs, ncount, npct)
                if p:
                    fo.write('\t'.join(p) + '\n')

# EOF