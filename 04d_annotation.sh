#!/bin/bash
annotator='python /home/users/kjyi/Library/ysj/SNP_annotation/02_CNVannotate_ju.161117.py'
for sample in varscan/*.indel; do
	name=$(basename $sample|sed 's/.indel//')
	$annotator $sample 2 yes > log/04d_$name.out
	#$annotator $sample > log/04b_$names.out
done
echo DONE


# change name
cd varscan
for anno in ./*_annoHg19.indel ;do
    target=$(echo $anno | sed 's/_annoHg19//')
    mv $anno $target
done
