#!/usr/bin/env Rscript
library("optparse")

# Parse options 
option_list = list(
    make_option(c("-i", "--input_bam_file"), type="character", default=NULL, 
              help="input bam file", metavar="character"),
    make_option(c("-c", "--control_bam_file"), type="character", default=NULL, 
              help="file containing names of control bam files", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default="ExomeDepth_DIR", 
              help="ExomeDepth output directory [default= %default]", metavar="character"),
    make_option(c("-r", "--reference"), type="character", default=NULL, 
              help="genome reference fasta file", metavar="character")
    )
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load libraries
suppressMessages(library("TxDb.Mmusculus.UCSC.mm10.knownGene"))
suppressMessages(library("plyranges"))
suppressMessages(library("ExomeDepth"))
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library("AnnotationHub"))
suppressMessages(library("AnnotationDbi"))

# Load main exomedepth function
source("aux1-exomedepth.R")

# Load annotation data from UCSC
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
mykeys <- keys(txdb, keytype = "EXONID")

exons.mm10 <- AnnotationDbi::select(txdb,
       keys = mykeys,
       keytype="EXONID",
       columns = c("EXONCHROM", "EXONSTART", "EXONEND","GENEID", "EXONID")
      )[,c(2,3,4,1,5)]

# Remove non-chromosomal contigs and mitochondrial chromosome. chrX and chrY should be removed as well if the controls have a different gender to the actual sample
keep <- grep(pattern = '^(chr[0-9]+|chrX|chrY)$', exons.mm10$EXONCHROM)
exons.mm10 <- exons.mm10[keep,]

# Generate GR object with annotation (used later...)
exons.mm10.GRanges <- GenomicRanges::GRanges(seqnames = exons.mm10$EXONCHROM,
                                                IRanges::IRanges(
                                                 start=exons.mm10$EXONSTART,
                                                 end=exons.mm10$EXONEND
                                                 ),
                                                names = exons.mm10$EXONID
                                            )

# Upload control bam files (to contrast with sample bam)
control_bams <- read.delim(file = opt$control_bam_file, header = FALSE)

myBams <- c(paste0('./04-realign/',opt$input_bam_file), control_bams[,1])

ExomeCount <- getBamCounts(bed.frame=exons.mm10,
                           bam.files=myBams,
                           include.chr=FALSE,
                           referenceFasta=opt$reference
                           )

### Identify the reference set of bam files (note the name change due to importing into R)
my.ref.samples <- c(control_bams[,1])

# Remove bam extension from input file name
sample_name <- gsub(".bam","",opt$input_bam_file)

dir.create(path = opt$outdir, showWarnings = FALSE)

# Run ExomeDepth
exomeDepth(my.test.num = opt$input_bam_file,
        outdir = opt$outdir,
        ExomeCount = ExomeCount,
        my.ref.samples = my.ref.samples,
        exons.mm10.GRanges = exons.mm10.GRanges,
        sample=sample_name
        ) 

print(paste("Analysis complete for", opt$input_bam_file))
