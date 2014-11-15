source("aa_nomeclature.R")

#load table
RNA_Codon <- read.table("RNA_Codon.txt")

#Count the number of different RNA sequences possible with a given AA sequence.
n_DNAStrings_from_AA <- function(sequence){
    #type aaa refers to when you have a 3 letters nomenclature for the AA spaced by "-".
    #type a refers to when you have a 1 letter nomenclature for the AA not spaced.
    
    if(grepl("-", sequence)){
        type = "aaa"
    }
    else{
        type = "a"
    }
    
    if(type=="aaa"){
        sequence <- aaa_to_a(sequence)
        sequence <- strsplit(sequence, NULL)[[1]]
        sequence <- sapply(sequence, function(pep){    
            length(grep(pep, RNA_Codon[,2]))    
        })
        return(prod(sequence))
        
    }
    
    if(type=="a"){
        sequence <- strsplit(sequence, NULL)[[1]]
        sequence <- sapply(sequence, function(pep){    
            length(grep(pep, RNA_Codon[,2]))    
        })
        return(prod(sequence))
    }
}