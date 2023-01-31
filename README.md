# USDIN_LAB_TK_32

### Clone git repository with all the required scripts. To do so, run the following command from your working directory 
```
git clone https://github.com/TriLab-bioinf/USDIN_LAB_TK_32.git
```
### Copy required script into your working directory by running the following command
```
# /path/to/the/command/USDIN_LAB_TK_32/copy_scrips.sh
# For instance:

./USDIN_LAB_TK_32/copy_scrips.sh
```

### Edit the configuration file config.cfg in your current working directory with a text editor and setup the following variables:

- READ : path to the folder containing the raw fastq files in a compress form (the pipeline expects the following suffix= fastq.gz).
- ADAPTERS : full path to file containing the sequencing adapter sequences in fasta format.
- REFERENCE : path to the reference genome in fasta format
- MERGE : set to 'true' or 'false' depending if you want to merge paired-end reads before mapping (true) or not.
- CONTROL_BAMS : File containing a list of bam files used as control for detecting structural variants or CNVs.

Example of config.cfg file
```
> cat config.cfg
# path ro reads directory
READ=/data/lorenziha/KAREN_USDIN/TK_32/READS

# path to adapters fasta file
ADAPTERS=/data/lorenziha/KAREN_USDIN/TK_32/adapters.fa

# Reference genome
REFERENCE=/data/lorenziha/KAREN_USDIN/TK_32/mm10.fa

# Merge reads before mapping? [true (default) or false]
MERGE=true

# File with names of bam files from control samples for exomeDepth (one file per line)
CONTROL_BAMS=control_bam_files.txt
```
Control bam files and their corresponding bam index files (.bai) should be in the same working directory as the one you are running the pipeline from.

### Run the pipeline as so:
```
sbatch ./bruce_pipeline.sh <sample name prefix> <configuration file>
```

Example:
```
sbatch ./bruce_pipeline.sh A20F_WES config.cfg
```

