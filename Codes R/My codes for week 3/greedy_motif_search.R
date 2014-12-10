library("kebabs")
library("stringdist")
source("most_prob_kmer.R")


greedy_motif_search <- function(sequences, k, t){
    if(length(sequences)==1){
        sequences <-gsub("\\s+", " ", sequences)
        sequences <- strsplit(sequences, " ")[[1]]
    }
    best_motifs <- as.vector(sapply(sequences, substring, 1, k))
    kernel <- spectrumKernel(k=k, normalized=FALSE)
    s1_kmers <- drop(getExRep(DNAString(sequences[1]), kernel))
    s1_kmers <- names(s1_kmers)
    n_of_k <- length(s1_kmers)
    
    kmers_set <- sapply(s1_kmers, function(x){
        
        kmers_set <- x
        profile <- consensusMatrix(DNAStringSet(x), as.prob=T, baseOnly=TRUE)
        profile <- list(profile[1:4,])[[1]]
        
        sapply(2:t,function(y){
            
            s2_kmers <- drop(getExRep(DNAString(sequences[y]), kernel))
            s2_kmers <- names(s2_kmers)
            
            motifs <- sapply(s2_kmers, function(z){
                zz <- strsplit(z, NULL)[[1]]
                prod <- sapply(1:k, function(w){
                    profile[profile=zz[w],w]
                })
                prod(prod)
            })
            motifs <- sort(motifs[motifs>0],decreasing = T)
            if(length(motifs)==0){
                b_motif<- best_motifs[y]
                
            }
            else{
                b_motif<- names(motifs)[1]
                
            }
            
            
            kmers_set[y] <<- b_motif
            profile <<- consensusMatrix(DNAStringSet(kmers_set), as.prob=T, baseOnly=TRUE)
            kmers_set
        })[[t-1]]
        
        
    })
    
    most_common <- as.vector(apply(kmers_set, 2, function(x){names(sort(table(x),decreasing=T))[1]}))
    
    hamming <- sapply(1:n_of_k, function(i){
        sum(stringdist(most_common[i], kmers_set[,i]))
    })
    return(kmers_set[,order(hamming)][,1])
}



sequences <- c("GGCGTTCAGGCA",     "AAGAATCAGTCA" ,    "CAAGGAGTTCGC",     "CACGTCAATCAC",     "CAATAATATTCG")
k <- 3
t <- 5
greedy_motif_search(sequences, k , t)
GreedyMotifSearch(sequences, k , t)

system.time({greedy_motif_search(sequences, k , t)})
system.time({GreedyMotifSearch(sequences, k , t)})

# greedy_motif_search(sequences,k,t)


data <- readLines("data/greedy_data.txt")
t <- as.numeric(strsplit(data[2], " ")[[1]][2])
k <- as.numeric(strsplit(data[2], " ")[[1]][1])
sequences <- data[3:27]
greedy_motif_search(sequences,k,t)


data <- readLines("data/dataset_159_5.txt")
t <- as.numeric(strsplit(data[1], " ")[[1]][2])
k <- as.numeric(strsplit(data[1], " ")[[1]][1])
sequences <- data[2:26]
GreedyMotifSearch(sequences,k,t)