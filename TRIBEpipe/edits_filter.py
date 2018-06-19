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
    parser.add_argument('-i', nargs = '+', required = True, metavar = 'TRIBE',
        help = 'editing events of TRIBE sample')
    parser.add_argument('-gDNA', required = True, metavar = 'genomic_DNA',
        help = 'editing events in genomic DNA, included in results')
    parser.add_argument('-wtRNA', required = True, metavar = 'wildtype_RNA',
        help = 'editing events in wildtype RNA-seq, excluded from results')
    # parser.add_argument('-c', nargs = '+', required = True,  metavar = 'Control',
    #     type = argparse.FileType('r'),
    #     help = 'edit events of control samples, genomic DNA, wildtype mRNA, etc.')
    parser.add_argument('-g', required = True, metavar = 'GTF',
        help = 'gene annotation in GFF format')
    parser.add_argument('-o', required = True, metavar = 'Output',
        help = 'file to save the results')
    parser.add_argument('--remove_tmp', action = 'strore_true',
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
        logging.info('fail, processing BED annotation')



def edits_filter(edits_tribe, edits_gDNA, edits_wtRNA, 
                     gtf, edits_tribe_filt, remove_tmp = False):
    """filtering edits"""
    assert isinstance(edits_tribe, str)
    assert isinstance(edits_gDNA, str)
    assert isinstance(edits_wtRNA, str)
    gene_bed = os.path.splitext(gtf)[0] + '.bed'
    if not os.path.exists(gene_bed):
        gtf2bed(gtf, gene_bed, feature = 'gene')

    # outdir
    if not os.path.exists(os.path.dirname(edits_tribe_filt)):
        os.makedirs(os.path.dirname(edits_tribe_filt))
    
    # include gDNA
    b1 = bed_filter(pybedtools.BedTool(edits_tribe), [edits_gDNA, ], exclude = False)

    # exclude wtRNA
    b2 = bed_filter(b1, [edits_wtRNA, ], exclude = True)

    # annotation
    bed_anno(b2, gene_bed, edits_tribe_filt)

    # remove temp files
    if remove_tmp:
        os.remove(gene_bed)

    return edits_tribe_filt


def main():
    args = get_args()
    edits_tribe = [i.name for i in args.i]
    # edits_control = [i.name for i in args.c]
    edits_gDNA = args.gDNA
    edits_wtRNA = args.wtRNA
    # edits_filter(edits_tribe, edits_control, args.g, args.o, args.remove_tmp)
    for i in args.i:
        edits_filter(i.name, edits_gDNA, edits_wtRNA, args.g,
                     args.o, args.remove_tmp)


if __name__ == '__main__':
    main()


## EOF

