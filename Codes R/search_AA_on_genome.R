library("Biostrings")

genome <- readLines("B_brevis.txt")

k = 9
L = 500
t = 3

kmer <- oligonucleotideFrequency(genome, k)
kmer <- kmer[kmer>=t]



search_AA_on_genome <- function(genome, k, L, t){
    
    genome <- DNAString(genome)
    sequence <- oligonucleotideFrequency(genome, k)
    sequence <- sequence[sequence>=t]
    sequence <- names(sequence)
    
   
    
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