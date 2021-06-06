#!/usr/bin/env bash

# conda activate dev-env

bamdir=/home/eric/Carpenter_Lab/Projects/COPD/BAMs/
python3 rMATS/rmats2sashimiplot-master/src/rmats2sashimiplot/rmats2sashimiplot.py \
	--b1 ${bamdir}chr_air_ctl_er051.bam,${bamdir}chr_air_ctl_er052.bam,${bamdir}chr_air_ctl_er053.bam \
	--b2 ${bamdir}chr_cig_ctl_er054.bam,${bamdir}chr_cig_ctl_er055.bam,${bamdir}chr_cig_ctl_er056.bam \
	-t RI \
	-e /home/eric/Carpenter_Lab/Projects/COPD/SashimiPlots/RI/RI.txt \
	--l1 Air \
	--l2 Smoke \
	--exon_s 1 \
	--intron_s 1 \
	-o /home/eric/Carpenter_Lab/Projects/COPD/SashimiPlots/RI
