#!/bin/bash
### varscan indel 

for mpileup in mpileup/*.vcf; do
	output=$(echo $mpileup | sed 's/mpileup/varscan/' | sed 's/vcf/indel/')
       	java -jar `which varscan` mpileup2indel $mpileup --min-var-freq 0.03 --output-vcf 1 > $output
done
echo DONE

