#!/bin/bash
### varscan
mkdir -p varscan
for sample in mpileup/*.vcf; do
	output=$(basename $sample)
	#sem -j-3 '/usr/java/jre1.8.0_91/bin/java -jar /home/users/tools/varscan2.4.2/VarScan.v2.4.2.jar mpileup2snp '$sample' --min-var-freq 0.03 --output-vcf 1 > varscan/'$output ' ; echo ' $output ' DONE'
	sem -j-3 'java -jar /home/users/tools/varscan2.4.2/VarScan.v2.4.2.jar mpileup2snp '$sample'\
		--min-var-freq 0.03 \
		--output-vcf 1 > varscan/'$output ' ; echo ' $output ' DONE'
done
sem --wait
echo VARSCAN calling DONE
