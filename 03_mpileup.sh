#!/bin/bash
ref=/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta
mkdir -p mpileup
for sample in bamfiles/*.bam; do
	output=$(basename $sample | sed 's/\.bam$/.vcf/')
	sem -j-3 samtools mpileup -B -d 300000 -q 20 -Q 20 -r MT -f $ref $sample '>' mpileup/$output
done
sem --wait
echo PILEUP DONE;mail.py 'pileup done'
