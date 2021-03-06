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

__author__ = "Ming Wang <wangm08@hotmail.com>"
__copyright__ = "2018 by Ming Wang <wangm08@hotmail.com>"
__license__ = "MIT"
__email__ = "wangm08@hotmail.com"
__version__ = "0.1"


import os, sys, re, argparse, logging
import pybedtools

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'edits_filter', description = 'Fitering editing events',
        epilog = 'Output: ')
    parser.add_argument('-i', required = True, 
        type = argparse.FileType('r'), metavar = 'TRIBE',
        help = 'editing events of TRIBE sample')
    parser.add_argument('-gDNA', nargs = '+', required = False, 
        type = argparse.FileType('r'), metavar = 'gDNA',
        help = 'editing events in genomic DNA, included')
    parser.add_argument('-wtRNA', nargs = '+', required = False,
        type = argparse.FileType('r'), metavar = 'wtRNA',
        help = 'editing events in wildtype RNA-seq, excluded')
    parser.add_argument('-g', required = True, metavar = 'GTF',
        help = 'gene annotation in GFF format')
    parser.add_argument('-o', required = True, metavar = 'Output',
        help = 'file to save the results')
    parser.add_argument('--remove_tmp', action = 'store_true',
        help = 'if specified, remove temp file, anno.bed')
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
    except IOError:
        logging.error('fail to parse gtf and bed files')



def bed_filter(bed_in, bed_excludes, exclude = True):
    """
    exclude bed_exclude records from bed_in
    bed_in: BedTool
    bed_excludes: files
    return: BedTool (or None)
    """
    assert isinstance(bed_in, pybedtools.bedtool.BedTool)
    assert isinstance(bed_excludes, list)
    if len(bed_excludes) >= 1:
        b = pybedtools.BedTool(bed_excludes.pop())
        if exclude:
            b2 = bed_in.intersect(b, v = True)
        else:
            b2 = bed_in.intersect(b, u = True)
        if len(bed_excludes) == 0:
            return b2
        else:
            return bed_filter(b2, bed_excludes, exclude = exclude)
    else:
        return None



def bed_anno(bed_in, bed_info, bed_out):
    """
    Annotate bed records by another bed file
    default: bed6 and bed6
    """
    assert isinstance(bed_in, pybedtools.bedtool.BedTool)
    b_info = pybedtools.BedTool(bed_info)
    t = bed_in.intersect(b_info, f = 0.9, loj = True, stream = True)
    t.sort().saveas(bed_out + '.tmp')
    # # why tail()?
    # # pybedtools issue #246 
    # # https://github.com/daler/pybedtools/issues/246
    # for b in t.sort().tail(1000000000): 

    try:
        dd = {}
        with open(bed_out + '.tmp') as fi, open(bed_out, 'wt') as fo:
            for line in fi:
                fs = line.strip().split('\t')
                fout = fs[:-6] + [fs[-3]]
                pos = ':'.join(fs[0:3])
                if pos in dd:
                    fo.write(',' + fs[-3])
                else:
                    if dd == {}:
                        fo.write('\t'.join(fout))
                    else:
                        fo.write('\n' + '\t'.join(fout))
                    dd[pos] = 1
            fo.write('\n')
        os.remove(bed_out + '.tmp')
    except IOError:
        logging.error('fail, processing BED annotation')



def edits_filter(edits_tribe, edits_gDNA, edits_wtRNA, 
                 gtf, edits_tribe_filt, remove_tmp = False):
    """filtering edits"""
    assert isinstance(edits_tribe, str)
    # assert isinstance(edits_gDNA, str)
    # assert isinstance(edits_wtRNA, str)
    gene_bed = os.path.splitext(gtf)[0] + '.bed'
    if not os.path.exists(gene_bed):
        gtf2bed(gtf, gene_bed, feature = 'gene') # convert annotation

    # outdir
    out_path = os.path.dirname(edits_tribe_filt)
    if not os.path.exists(out_path) and not out_path == '':
        os.makedirs(out_path)
    
    ## exclude gDNA
    if isinstance(edits_gDNA, str):
        if os.path.exists(edits_gDNA):
            b1 = bed_filter(pybedtools.BedTool(edits_tribe), 
                            [edits_gDNA, ], exclude = True)
        else:
            b1 = pybedtools.BedTool(edits_tribe)
    elif isinstance(edits_gDNA, list):
        b_in = pybedtools.BedTool(edits_tribe)
        b_out = b_in
        for gDNA in edits_gDNA:
            if os.path.exists(gDNA):
                b_out = bed_filter(b_in, [gDNA, ], exclude = True)
            else:
                b_out = b_in
            # cycle
            b_in = b_out
        b1 = b_out
    else:
        b1 = pybedtools.BedTool(edits_tribe)

    # exclude wtRNA
    if isinstance(edits_wtRNA, str):
        if os.path.exists(edits_wtRNA):
            b2 = bed_filter(b1, [edits_wtRNA, ], exclude = True)
        else:
            b2 = b1
    elif isinstance(edits_wtRNA, list):
        b_in = b1
        b_out = b1
        for wt in edits_wtRNA:
            if os.path.exists(wt):
                b_out = bed_filter(b_in, [wt, ], exclude = True)
            else:
                b_out = b_in
            # cycle
            b_in = b_out
        b2 = b_out
    else:
        b2 = b1


    # annotation
    bed_anno(b2, gene_bed, edits_tribe_filt)

    # remove temp files
    if remove_tmp:
        os.remove(gene_bed)

    return edits_tribe_filt


def main():
    args = get_args()
    # edits_tribe = [i.name for i in args.i]
    edits_tribe = args.i.name
    edits_gDNA = None if args.gDNA is None else [i.name for i in args.gDNA]
    edits_wtRNA = None if args.wtRNA is None else [i.name for i in args.wtRNA]
    for i in args.i:
        edits_filter(edits_tribe, edits_gDNA, edits_wtRNA, args.g,
                     args.o, args.remove_tmp)


if __name__ == '__main__':
    main()


## EOF
