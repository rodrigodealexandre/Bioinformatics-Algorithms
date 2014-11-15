source("aa_nomeclature.R")

#load table
RNA_Codon <- read.table("RNA_Codon.txt")

protein_to_RNA <- function(sequence){
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
        if(length(sequence) <= 12){
            sequence <- expand.grid(lapply(sequence, function(pep){    
                RNA_Codon[grep(pep, RNA_Codon[,2]),1]
            }))
            sequence <- apply(sequence, 1,function(seq){paste(seq, collapse = "")})
            return(sequence)
        }
        else{
            return("sequence is to big")
        }
    }
    
    if(type=="a"){
        sequence <- strsplit(sequence, NULL)[[1]]
        if(length(sequence) <= 12){
            sequence <- expand.grid(lapply(sequence, function(pep){    
                RNA_Codon[grep(pep, RNA_Codon[,2]),1]
            }))
            sequence <- apply(sequence, 1,function(seq){paste(seq, collapse = "")})
            return(sequence)
        }
        else{
            return("sequence is to big")
        }
    }
}