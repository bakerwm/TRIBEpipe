#!/bin/bash

# convert bam to acgt count matrix

[[ $# -lt 3 ]] &&  echo "Usage: get_edits.sh <in.bam> <genome.fa> <cutoff|20>" && exit 1
inbam=$1
gfa=$2
cutoff=$3
[[ $cutoff -lt 1 ]] && echo "cutoff should be >0" && exit 1
[[ -z $4 ]] && strand=1 || strand=$4


# genome index
[[ ! -f ${gfa}.fai ]] && samtools faidx $gfa


if [[ $strand == '1' ]] ; then
    # strandness
    # forward strand: A->G
    samtools view -f 16 -bhS $inbam | \
        samtools mpileup -f $gfa - #| \
        #sequenza-utils pileup2acgt -n $cutoff -p - 

    # reverse strand: T->C
    samtools view -F 16 -bhS $inbam | \
        samtools mpileup -f $gfa - #| \
        #sequenza-utils pileup2acgt -n $cutoff -p - 
else
    # strandless
    samtools view -bhS $inbam | \
        samtools mpileup -f $gfa - | \
        sequenza-utils pileup2acgt -n $cutoff -p -
fi
##
