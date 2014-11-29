#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")

#Implement MOTIFENUMERATION.

library("gtools")
library("Biostrings")

count_difference <- function(kmer1, kmer2, d){
    if(length(kmer1) == length(kmer2)){
        count <- length(kmer1)
        for(k in 1:length(kmer1)){
            if(kmer1[k] != kmer2[k]){
                count <- count - 1
            }
        }
        difference <- length(kmer1) - count
        if(difference <= d){
            return(TRUE)
        }
        else{
            return(FALSE)
        }
    }
    else(cat("K-mer1 and K-mer2 has different lenghts"))
}

#finding all kmers in each sequence
sequence_kmers <- function(sequence, k){
    k_mers <- lapply(sequence,function(x){
        seq_loop_size <- length(DNAString(x))-k+1
        
        kmers <- sapply(1:seq_loop_size, function(z){
            y <- z + k -1
            kmer <- substr(x=x, start=z, stop=y)
            return(kmer)
        })
        return(kmers)
    })
    return(k_mers)
}

#apply doesn't seems to impruve speed
motifenumerator <- function(DNAs, k, d){
    bases <- c("A","T","G","C")
    all_combinations <- permutations(4,k,bases, repeats.allowed=T)
    seq_kmers <- sequence_kmers(DNAs, k)
    
    motifs <- lapply(seq_kmers, function(x){
        is_true <- matrix(ncol=k)
        count <- 0
        lapply(x, function(y){
            y <- strsplit(y, split="")[[1]]
            count <<- count +1
            is_true <<- rbind(all_combinations[apply(all_combinations,1, count_difference, y, d),],is_true)
        })
        all_combinations <<- unique(is_true[!is.na(is_true[,1]),])
        return(all_combinations)
    })[[length(seq_kmers)]]
    
    return(apply(motifs,1,paste,collapse=""))
}

#work as the same speed as the previeous function
motifenumerator <- function(DNAs, k, d){
    bases <- c("A","T","G","C")
    all_combinations <- permutations(4,k,bases, repeats.allowed=T)
    seq_kmers <- sequence_kmers(DNAs, k)
    for(x in seq_kmers){
        is_true <- matrix(ncol=k)
        count <- 0
        for(y in x){
            y <- strsplit(y, split="")[[1]]
            count <- count +1
            is_true <- rbind(all_combinations[apply(all_combinations,1, count_difference, y, d),],is_true)
        }
        all_combinations <- unique(is_true[!is.na(is_true[,1]),])
    }
    return(apply(all_combinations,1,paste,collapse=""))
}


# example <- readLines("data/motif_enumeration_data.txt")
# numbers <- as.numeric(strsplit(example[2], split=" ")[[1]])
# k <- numbers[1]
# d <- numbers[2]
# DNAs <- example[3:8]
# 
# 
# example <- readLines("data/dataset_156_7.txt")
# numbers <- as.numeric(strsplit(example[1], split=" ")[[1]])
# k <- numbers[1]
# d <- numbers[2]
# DNAs <- example[2:7]

# DNAs <- DNAStringSet(c("ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"))
# DNAs <- c("ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT")
# k=3
# n=1
# 
# cat(motifenumerator2(DNAs, k, d))
