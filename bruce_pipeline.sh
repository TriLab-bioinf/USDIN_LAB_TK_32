#!/usr/bin/bash
#SBATCH --mem=32g --cpus-per-task=16 --gres=lscratch:10 --time=1-00:00:00
#From http://seqanswers.com/forums/showthread.php?t=42776

set -o errexit 

# Command example:
# sbatch --time=24:00:00 bruce_pipeline.sh K3F_WES config.cfg

SAMPLE=$1
CONFIG=$2

if [[ ! -e ${CONFIG} ]]; then
    CONFIG=config.cfg
fi        

source ${CONFIG}

module load  bbtools/38.96 samtools GATK/3.8-1 picard

#####################################
# Trim adapters from paired end reads
#####################################

if [[ ! -d "01-trimming" ]]; then
    mkdir "01-trimming"
fi

# bbduk.sh -Xmx1g in1=read1.fq in2=read2.fq out1=clean1.fq out2=clean2.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

bbtools bbduk in1=${READ}/${SAMPLE}_1.fastq.gz in2=${READ}/${SAMPLE}_2.fastq.gz ref=${ADAPTERS} out1=./01-trimming/${SAMPLE}_R1_trimmed.fastq.gz out2=./01-trimming/${SAMPLE}_R2_trimmed.fastq.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo -Xmx1g > ./01-trimming/${SAMPLE}_adapter_trim.log 2>&1  

# Quality trimming

# bbduk.sh -Xmx1g in=read1.fq in=read2.fq out1=clean1.fq out2=clean2.fq qtrim=rl trimq=10

bbtools bbduk in1=./01-trimming/${SAMPLE}_R1_trimmed.fastq.gz in2=./01-trimming/${SAMPLE}_R2_trimmed.fastq.gz out1=./01-trimming/${SAMPLE}_R1_trimmed_qc.fastq.gz out2=./01-trimming/${SAMPLE}_R2_trimmed_qc.fastq.gz qtrim=rl trimq=20 -Xmx1g > ./01-trimming/${SAMPLE}_quality_trim.log 2>&1

#############################
# Merge PE reads with bbmerge
#############################

if [[ ! -d "02-merge" ]]; then
    mkdir "02-merge"
fi

# bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>

bbtools bbmerge in1=./01-trimming/${SAMPLE}_R1_trimmed_qc.fastq.gz in2=./01-trimming/${SAMPLE}_R2_trimmed_qc.fastq.gz out=./02-merge/${SAMPLE}_merged.fastq.gz outu1=./02-merge/${SAMPLE}_R1_unmerged.fastq.gz outu2=./02-merge/${SAMPLE}_R2_unmerged.fastq.gz strict=t > ./02-merge/${SAMPLE}_merge.log 2>&1

REF=`basename ${REFERENCE}`
REF_BASENAME=`echo ${REF%.f*a}`

############################################################################################
# Ensure that reference fasta file is in current directory, otherwise create a symbolic link
############################################################################################

if [[ ! -e ${REF} ]]; then
    ln -s ${REFERENCE} ${REF}
fi        

#########################################
# Ensure that reference dict is present 
#########################################

if [[ ! -e  ${REF_BASENAME}.dict ]]; then
    samtools dict ${REFERENCE} -o {REF_BASENAME}.dict       
fi

##############################################
# Ensure that reference index files is present
##############################################

if [[ ! -e  ${REF}.faidx ]]; then
    samtools faidx ${REF} -o ${REF}.fai
fi        

##################################
# Make mm10 reference genome file
##################################

if [[ ! -d ref ]]; then
    bbtools bbmap -Xmx24g ref=${REF}
fi            

#############################################################################
# Mapping (no need to specify genome as mm10 the only one available)
#############################################################################

if [[ ! -d "03-mapping" ]]; then
    mkdir "03-mapping"
fi

if [[ ! ${MERGE} ]]; then
    bbtools bbmap -Xmx24g in1=./01-trimming/${SAMPLE}_R1_trimmed_qc.fastq.gz in2=./01-trimming/${SAMPLE}_R2_trimmed_qc.fastq.gz out=./03-mapping/${SAMPLE}.sam ref=${REFERENCE} > ./03-mapping/${SAMPLE}_map.log 2>&1 
else
    bbtools bbmap -Xmx24g in=./02-merge/${SAMPLE}_merged.fastq.gz out=./03-mapping/${SAMPLE}.sam ref=${REFERENCE} > ./03-mapping/${SAMPLE}_map.log 2>&1
fi

#########################################
# Convert to bam and index for IGV viewer
#########################################

samtools sort -o ./03-mapping/${SAMPLE}_sorted.bam ./03-mapping/${SAMPLE}.sam
samtools index -o ./03-mapping/${SAMPLE}_sorted.bam.bai  ./03-mapping/${SAMPLE}_sorted.bam

###############
# Add RG field
###############

java -jar $PICARD_JAR AddOrReplaceReadGroups -I ./03-mapping/${SAMPLE}_sorted.bam -O ./03-mapping/${SAMPLE}_sorted_RG.bam -ID ${SAMPLE} -LB lib1 -PL ILLUMINA -PU unit1 -SM ${SAMPLE} > ./03-mapping/${SAMPLE}_add_group_id.log 2>&1

samtools index -o ./03-mapping/${SAMPLE}_sorted_RG.bam.bai ./03-mapping/${SAMPLE}_sorted_RG.bam

###########################
# Generate realigned table
###########################

if [[ ! -d "04-realign" ]]; then
    mkdir "04-realign"
fi

GATK --mem=16g RealignerTargetCreator -R mm10.fa -I ./03-mapping/${SAMPLE}_sorted_RG.bam -o ./04-realign/${SAMPLE}_realigner.intervals > ./04-realign/${SAMPLE}_realign_target_creator.log 2>&1

###########################
# Run indel realigner
###########################

GATK --mem=16g IndelRealigner -R mm10.fa -I ./03-mapping/${SAMPLE}_sorted_RG.bam --targetIntervals ./04-realign/${SAMPLE}_realigner.intervals -o ./04-realign/${SAMPLE}_realigned.bam > ./04-realign/${SAMPLE}_indel_realigner.log 2>&1

#############################################
# Extract variants and search for INDELs only
#############################################

if [[ ! -d "05-mpileup" ]]; then
    mkdir "05-mpileup"
fi

module load bcftools

bcftools mpileup -f mm10.fa ./04-realign/${SAMPLE}_realigned.bam | bcftools call -mv -Ob -o ./05-mpileup/${SAMPLE}_calls.bcf

bcftools mpileup -f  mm10.fa ./04-realign/${SAMPLE}_realigned.bam --indel-size 10000 | bcftools call -mv -Ob -o ./05-mpileup/${SAMPLE}_10kb_indels.bcf

# bcftools view -i '%QUAL>=20' -v indels  -q ./05-mpileup/${SAMPLE}_calls.bcf | bedtools intersect -a ~/Data/Ref_data/mm10_unique_RefGenes_NM_exon_coords.bed -b stdin -wa -wb

############################
# Run ExomeDepth on bam file
############################

module load R/4.2.2

# Remove chromosome M from bam file
sh ./aux-1_remove_chrM.sh ${SAMPLE}
# Run ExomeDepth
Rscript ./run_exomeDepth.R -i ${SAMPLE}_realigned_nochrM.bam -c ${CONTROL_BAMS} -o 07-exomeDepth -r ${REFERENCE}

