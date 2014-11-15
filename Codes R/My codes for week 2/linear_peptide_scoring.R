#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")
source("aa_nomeclature.R")
source("linear_spectrum.R")
source("peptide_score.R")
source("cyclopeptide_mass_spectrum.R")

linear_peptide_scoring <-function(sequence, spectrum){
    if(length(spectrum) == 1){
        spectrum <- strsplit(spectrum, "\\s" )[[1]]
    }
    
    sequence <- linear_spectrum(sequence)
    score <- peptide_score(sequence, spectrum)
    return(score)
}

file <- readLines("data/dataset_4913_1 (1).txt")
sequence <- file[1]
spectrum <- file[2]
score <- linear_peptide_scoring(sequence, spectrum)
cat(score)