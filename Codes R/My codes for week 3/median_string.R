#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")

#Implement DistanceBetweenPatternAndStrings
library("Biostrings")
source("d_pattern_dna.R")
source("motifenumerator.R")

#Implement MEDIANSTRING
median_string <- function(k, sequences, max_d=3){
    if(length(sequences)==1){
        sequences <- strsplit(sequences, " ")[[1]]
    }
    
    all_kmers <- unlist(sequence_kmers(sequences, k))
    all_kmers <- unique(all_kmers)

    kmers <- sapply(all_kmers, function(x){sum(d_pattern_dna(x, sequences, max_d))})

    kmers <- sort(kmers)
    return(kmers)
}



# sequences <- "AAATTGACGCAT GACGACCACGTT CGTCAGCGCCTG GCTGAGCACCGG AGTACGGGACAG"
# k <- 3
# 
# file <- readLines("data/dataset_158_9.txt")
# k <- as.numeric(file[1])
# sequences <- file[2:length(file)]
# median_string(k,sequences)


