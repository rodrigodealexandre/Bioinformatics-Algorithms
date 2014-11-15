source("protein_to_RNA.R")

library("Biostrings")

search_AA_on_genome <- function(sequence, genome){
    genome <- DNAString(genome)
    sequence <- RNAStringSet(protein_to_RNA(sequence))
    r_sequence <- reverseComplement(sequence)
    sequence <- RNAStringSet(c(sequence, r_sequence))
    sequence <- as.character(DNAStringSet(sequence))
    
    x <- startIndex(matchPDict(DNAStringSet(sequence), genome))
    names(x) <- sequence
    
    sequence <- vapply(x,FUN.VALUE = 1, FUN=function(seq){
        
        if(is.integer(seq[1])){
            return(length(seq)[[1]])
        }
        else{
            return(0)
        }
    }
    )
    sequence <- sequence[sequence >0]
    sequence <- sequence[order(sequence,decreasing = T)]
    return(sequence)
    
}
