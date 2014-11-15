#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")
source("aa_nomeclature.R")
source("cyclopeptide_mass_spectrum.R")

peptide_score <- function(sequence, spectrum){
    if(length(sequence) == 1){
        if(grepl("-", sequence)){
            sequence <- aaa_to_a(sequence)
        }
        sequence <- cyclopeptide_mass_spectrum(sequence)
    }
    
    if(length(spectrum) == 1){
        spectrum <- strsplit(spectrum, "\\s" )[[1]]
    }
    
    score <- sum(table(spectrum[spectrum %in% sequence]))
    return(score)
}


file <- readLines("data/dataset_4913_1 (1).txt")
sequence <- file[1]
spectrum <- file[2]
score <- peptide_score(sequence, spectrum)
cat(score)