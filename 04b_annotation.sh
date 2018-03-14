#!/bin/bash
ANNOTATE='python ~/Library/ysj/SNP_annotation/01_SNPannotate_ju_withGeneStrand.py'
for sample in varscan/*.vcf; do
	names=$(basename $sample | sed 's/.vcf//')
	echo "python ~/Library/ysj/SNP_annotation/01_SNPannotate_ju_withGeneStrand.py $sample > log/04b_$names.out"
	python ~/Library/ysj/SNP_annotation/01_SNPannotate_ju_withGeneStrand.py $sample > log/04b_$names.out
done
wait
echo DONE



# change name
cd varscan
for annot in *.annot; do
	out=$(echo $annot|sed 's/vcf.annot/snp/')
	mv $annot $out
done
