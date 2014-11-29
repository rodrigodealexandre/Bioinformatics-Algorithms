#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")

library("gtools")
library("Biostrings")
library("kebabs")

#Implement Profile-most probable k-mer in Text
most_prob_kmer <- function(sequence, profile, k){
    bases <- c("A","C", "G", "T")
    k <- as.numeric(k)
    if(length(profile)==4){
        profile <- matrix(as.numeric(unlist(lapply(profile, function(x){strsplit(x, " ")[[1]]}))), nrow=4, byrow=T)
    }
    rownames(profile) <- bases
        
    #all_kmers <- permutations(4,k,bases, repeats.allowed=T)
    #all_kmers <- unique(sequence_kmers(sequence,k)[[1]])
    
    kernel <- spectrumKernel(k=k, normalized=FALSE)
    all_kmers <- drop(getExRep(DNAString(sequence[1]), kernel))
    all_kmers <- unique(names(all_kmers))
    
    all_kmers <- matrix(all_kmers, nrow=length(all_kmers))
    all_kmers <- t(apply(all_kmers,1,function(x){strsplit(x, NULL)[[1]]}))
    
    all_k_profile <- t(apply(all_kmers, 1, function(x){
        sapply(1:k, function(y){
            profile[profile=x[y],y]
        })
    }))
    all_k_profile <- data.frame(all_k_profile)
    all_k_profile[,k+1] <- apply(all_k_profile,1,prod)
    
    sorted_profile <- order(all_k_profile[,k+1], decreasing=T)
    sorted_kmers_names <- all_kmers[sorted_profile,]
    sorted_profile <- all_k_profile[sorted_profile,]
    sorted_kmers <- sorted_profile[,k+1]
    names(sorted_kmers) <- apply(sorted_kmers_names,1,paste,collapse="")
    return(head(sorted_kmers,length(sequence)))
}



# data <- readLines("data/profile_most_1.txt")
# profile <- matrix(data[4:7], nrow=4)
# k <- as.numeric(data[3])
# sequence <- data[2]
# most_prob_kmer(sequence, profile, k)

# 
# data <- readLines("data/dataset_159_3.txt")
# profile <- matrix(data[3:length(data)], nrow=4)
# k <- as.numeric(data[2])
# sequence <- data[1]
# most_prob_kmer(sequence, profile, k)


