#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")


library("kebabs")
source("most_prob_kmer.R")

greedy_motif_search <- function(sequences, k, t){
    if(length(sequences)==1){
        sequences <-gsub("\\s+", " ", sequences)
        sequences <- strsplit(sequences, " ")[[1]]
    }
    
    kernel <- spectrumKernel(k=k, normalized=FALSE)
    s1_kmers <- drop(getExRep(DNAString(sequences[1]), kernel))
    s1_kmers <- names(s1_kmers)
    most_probable_k <- sapply(s1_kmers, function(x){
        profile <- as.character(subject(pairwiseAlignment(rep(x,t), sequences,type="global-local", 
                                                          gapOpening=-Inf)))
        hamming <- stringDist(unique(profile), method = "hamming")
        c(profile, sum(hamming))
        
    })
    
    most_probable_k <- t(most_probable_k[,order(as.numeric(most_probable_k[t+1,]))])
    
    size <- length(most_probable_k[,1])
    
    table <- sapply(sequences, function(seq){
        profiles <- apply(most_probable_k, 1, function(z){
            matrix <- consensusMatrix(DNAStringSet(z[1:t]), as.prob=T, baseOnly=TRUE)
            list(matrix[1:4,])
        })
        
        names <- sapply(1:length(profiles), function(p){
            most_probable <- most_prob_kmer(seq, profiles[[p]][[1]], k)
        })
        return(c(names(names),names))
    })
    the_kmers <- table[1:size,]
    
    probs_only <- table[(size+1):(size+size),]
    t_prob <- cbind(the_kmers, t_prob=apply(probs_only, 1, function(w){prod(as.numeric(w))}))
    return(t_prob[order(as.numeric(t_prob[,t+1]), decreasing=T),])
}


# sequences <- "GGCGTTCAGGCA     AAGAATCAGTCA     CAAGGAGTTCGC     CACGTCAATCAC     CAATAATATTCG"
# k <- 3
# t <- 5
# greedy_motif_search(sequences,k,t)
# 

data <- readLines("data/greedy_data.txt")
t <- as.numeric(strsplit(data[2], " ")[[1]][2])
k <- as.numeric(strsplit(data[2], " ")[[1]][1])
sequences <- data[3:27]
greedy_motif_search(sequences,k,t)
