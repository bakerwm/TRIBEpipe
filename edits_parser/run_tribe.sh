#!/bin/bash

# TRIBE, HyperTRIBE analysis
#
# Extract A to I (G) conversion positions from RNA-seq data set
#
# Criteria:
#
# TRIBE sample:
# 1. count > N (20)
# 2. freqG > 10%
# 
# genomic DNA:
# 1. count > N (?)
# 2. A% > 80%, G% = 0
# 
# background mRNA (wt)
# 1. count > N (10)
# 2. freqG > 10%
#
#

genome_dir="/data/genome"
multi_cores=16 # number of threads to use


## command line
## group: genomic_DNA, wt_mRNA, TRIBE
## ref_fa: 
## ref_gtf:
## fastq:

[[ $# != 4 ]] && echo "Usage: run_tribe.sh <group:gDNA|wt_mRNA|TRIBE> <refname> <outdir> <fastq>" && exit 1
group=$1
refname=$2
outdir=$3
infq=$4

## parameters
# <grup>
group=$(echo $group | tr '[A-Z]' '[a-z]')
case $group in
    gdna)        smp_type=DNA; min_depth=10; min_pct=80 ;;
    wt_mrna)     smp_type=RNA; min_depth=10; min_pct=20 ;;
    tribe)       smp_type=RNA; min_depth=20; min_pct=20 ;;
    *)           echo "error, unknown group $group" && exit 1 ;;
esac

# <refname>
case $refname in
    dm3|dm6|mm10|hg19|hg38)     genome="${genome_dir}/$refname/bigZips/$refname.fa"; 
                                g_index="${genome_dir}/$refname/STAR_index"; 
                                g_gtf="${genome_dir}/$refname/annotation_and_repeats/$refname.refseq.gtf" ;;
    *)      echo "error, unknown genome: $refname" && exit 1 ;;
esac

# <outdir> <infq>
[[ ! -d $outdir ]] && mkdir -p $outdir
[[ ! -f $infq ]] && echo "error, input fastq not exists" && exit 1
prefix=$(basename $infq | sed -e 's/\.f[astq]*//' -e 's/\.gz$//')

# <scripts>
bin_dir=$(dirname $(readlink -fs $0))
trim="${bin_dir}/trim.py"
map="${bin_dir}/map.py"
edits_parser="${bin_dir}/edits_parser.py"

## trimming
trim_dir="${outdir}/input_reads"
python $trim -i $infq -o $trim_dir -m 19 -q 25 --cut 7 --threads $multi_cores 
clean_fq="${trim_dir}/${prefix}.clean.fq"

## mapping
mapping_dir="${outdir}/mapping"
python $map -i $clean_fq -n demo -o $mapping_dir -g $refname -x $g_index -t STAR -p $multi_cores --rmdup
bam_in=${mapping_dir}/${prefix}/${prefix}.nodup.bam

## extract edits
edits_dir="${outdir}/edits"
[[ ! -d $edits_dir ]] && mkdir -p $edits_dir
python $edits_parser -t $smp_type -i $bam_in -g $genome -o $edits_dir/${prefix}.edits.bedgraph --min_depth $min_depth --min_pct $min_pct

## filtering



