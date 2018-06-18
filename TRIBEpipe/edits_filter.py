#!/usr/bin/env python
"""
filt edit sites 

Compute each sample respectively:
TRIBE - gDNA - wt_mRNA = sites

TRIBE sample:
1. count > N (20)
2. freqG > 10%

genomic DNA:
A% > 80%, G% = 0

background mRNA:
count > N (10)
freqG > 10%

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

import os, sys, re, argparse, tempfile, logging
import pybedtools



def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'edits_filter', description = 'Fitering editing events',
        epilog = 'Output: ')
    parser.add_argument('-i', required = True, metavar = 'TRIBE',
        help = 'edit events of TRIBE sample')
    parser.add_argument('-c', nargs = '+', required = True,  metavar = 'Control',
        type = argparse.FileType('r'),
        help = 'edit events of control samples, genomic DNA, wildtype mRNA, etc.')
    parser.add_argument('-g', required = True, metavar = 'GTF',
        help = 'gene annotation in GFF format')
    parser.add_argument('-o', required = True, metavar = 'Output',
        help = 'file to save the results')
    args = parser.parse_args()
    return args



def fetch_ensembl_anno(build = 'dm3', feature = 'gene', outfmt = 'bed'):
    """
    retrieve the features
    output: BED, GTF

    BED:
    chr  start  end  name  score  strand

    GTF:
    chr  source  start  end  .  strand  .  description
    """
    pass


def fetch_geneid(fn):
    """
    fetch gene id from GTF records
    gene_id "FBgn0267431"; gene_name "Myo81F"; gene_source "FlyBase"; gene_biotype "protein_coding";
    """
    assert isinstance(fn, str)
    ps = fn.split(';')
    dn = {}
    for i in fn.split(';'):
        if len(i) < 8: # shortest field
            continue
        id = i.strip().split(' ')
        id = re.sub('"', '', id[1])
        if 'gene_id' in i:
            dn['gene_id'] = id
        elif 'gene_name' in i:
            dn['gene_name'] = id
        elif 'transcript_id' in i:
            dn['transcript_id'] = id
        elif 'transcript_name' in i:
            dn['transcript_name'] = id
        elif 'exon_id' in i:
            dn['exon_id'] = id
        elif 'gene_biotype' in i:
            dn['gene_biotype'] = id
        else:
            pass
    return dn



def gtf2bed(in_gtf, out_bed, feature = True):
    """
    convert GTF to BED
    feature: gene  transcript  CDS  exon  five_prime_utr  three_prime_utr  start_codon
    """
    # n = 1
    groups = ['gene', 'transcript', 'CDS', 'exon', 'five_prime_utr', 
              'three_prime_utr', 'start_codon']
    try:
        with open(in_gtf, 'rt') as fi, open(out_bed, 'wt') as fo:
            for line in fi:
                fs = line.strip().split('\t')
                if len(fs) < 9: continue
                if not fs[2] in groups: continue
                chr = fs[0] if 'chr' in fs[0] else 'chr' + fs[0]
                start = int(fs[3]) - 1
                dt = fetch_geneid(fs[8])
                name = dt['gene_name'] if 'gene_name' in dt else dt['gene_id']
                fout = [chr, start, fs[4], name, 100, fs[6]]
                fout = list(map(str, fout))
                if feature is True or fs[2] == feature:
                    fo.write('\t'.join(fout) + '\n')
                else:
                    continue
                # if n > 10: break
                # n += 1
    except IOError:
        logging.error('fail to parse gtf and bed files')




def bed_filter(bed_in, bed_excludes):
    """
    exclude bed_exclude records from bed_in
    bed_in: BedTool
    bed_excludes: files
    """
    assert isinstance(bed_in, pybedtools.bedtool.BedTool)
    assert isinstance(bed_excludes, list)
    if len(bed_excludes) >= 1:
        b = pybedtools.BedTool(bed_excludes.pop())
        b2 = bed_in.intersect(b, v = True)
        if len(bed_excludes) == 0:
            return b2
        else:
            return bed_filter(b2, bed_excludes)
    else:
        return None



def bed_anno(bed_in, bed_info, bed_out):
    """
    Annotate bed records by another bed file
    default: bed6 and bed6
    """
    assert isinstance(bed_in, pybedtools.bedtool.BedTool)
    # b_in = pybedtools.BedTool(bed_in)
    b_info = pybedtools.BedTool(bed_info)
    t = bed_in.intersect(b_info, wa = True, wb = True, stream = True)
    try:
        with open(bed_out, 'wt') as fo:
            for b in t:
                fs = str(b).strip().split('\t')
                fs2 =  fs[:-6] + [fs[-3]] # add feature to the end
                fo.write('\t'.join(fs2) + '\n')
    except IOError:
        loging.info('fail, processing BED annotation')



def edits_filter(edits_tribe, edits_control, gtf, edits_tribe_filt):
    """filtering edits"""
    # convert GTF to BED, filt feature
    gene_bed = tempfile.NamedTemporaryFile()
    gtf2bed(gtf, gene_bed.name, feature = 'gene')
    
    # annotate
    ex_files = [i.name for i in edits_control]
    b1 = bed_filter(pybedtools.BedTool(edits_tribe), ex_files)
    bed_anno(b1, gene_bed.name, edits_tribe_filt)

    # remove temp files
    os.remove(gene_bed.name)

    return edits_tribe_filt


def main():
    args = get_args()
    edits_filter(args.i, args.c, args.g, args.o)
    # # convert GTF to BED, filt feature
    # gene_bed = tempfile.NamedTemporaryFile()
    # gtf2bed(args.g, gene_bed.name, feature = 'gene')
    # # gene_bed = 'Drosophila_melanogaster.BDGP6.92.bed'

    # # annotate
    # ex_files = [i.name for i in args.c]
    # b1 = bed_filter(pybedtools.BedTool(args.i), ex_files)
    # bed_anno(b1, gene_bed.name, args.o)

    # # remove temp files
    # os.remove(gene_bed.name)


if __name__ == '__main__':
    main()


## EOF

