#load table
aa_nomeclature <- read.table("aa_nomeclature.txt")

#Convert a 3 letters to a 1 letter peptide
aaa_to_a <- function(sequence){
    sequence <- strsplit(sequence, "-")[[1]]
    sequence <- sapply(sequence, function(pep){    
        aa_nomeclature[grep(pep, aa_nomeclature[,2]),1]})
    sequence <- paste(sequence, collapse = "")
    return(sequence)
}


#Convert a 1 letter to a 3 letters peptide
a_to_aaa <- function(sequence){
    sequence <- strsplit(sequence, NULL)[[1]]
    sequence <- sapply(sequence, function(pep){    
        aa_nomeclature[grep(pep, aa_nomeclature[,1]),2]})
    sequence <- paste(sequence, collapse = "-")
    return(sequence)
}