#!/bin/bash
# Basic Mapping process, made by AYH
# /home/users/kjyi/bin/ayhaligner
echo START
ayhaligner \
	-I filenames.txt \
	-d sampleids.txt \
	-O bamfiles \
	-@ 5 \
	-n 4 \
	-m 8g \
	-g \
	-a kijong.yi@kaist.ac.kr

