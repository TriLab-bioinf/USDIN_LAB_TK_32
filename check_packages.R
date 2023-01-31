#!/usr/bin/env Rscript

my_repo = "https://cloud.r-project.org"

if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager", repos=my_repo)
}

if (!require("plyranges", quietly = TRUE)){
    BiocManager::install("plyranges")
}

if (!require("IRkernel", quietly = TRUE)){
    install.packages("IRkernel", repos=my_repo)
}

if (!require("RIdeogram", quietly = TRUE)){
    install.packages("RIdeogram", repos=my_repo)
}



