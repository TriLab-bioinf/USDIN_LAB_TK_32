#!/usr/bin/bash

SAMPLE=$1

module load samtools

samtools view -h ./04-realign/${SAMPLE}_realigned.bam |grep -v 'chrM'  | samtools view -hb - > ./04-realign/${SAMPLE}_realigned_nochrM.bam
samtools index -o ./04-realign/${SAMPLE}_realigned_nochrM.bam.bai  ./04-realign/${SAMPLE}_realigned_nochrM.bam


