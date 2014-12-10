#------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 4")
source("composition_k.R")
source("overlap_graph.R")

de_bruijn_graph <- function(sequence, k){
    k <- as.numeric(k)-1
    kmers <- composition_k(sequence, k)
    kmers <- sort(kmers)
    overlap <- overlap_graph(kmers,sequence)
    
    return(unique(overlap))
    
}

# sequence <- "AAGATTCTCTAAGA"
# k <- 4

# data <- readLines("data/De_Bruijn_Graph_from_a_String.txt")
# sequence <- data[3]
# k <- as.numeric(data[2])


# data <- readLines("data/dataset_199_6.txt")
# sequence <- data[2]
# k <- as.numeric(data[1])
# 
# 
# result <- de_bruijn_graph(sequence, k)
# write.table(result, file="output/debruijngraph_result.txt", quote=F, 
#             sep="\n", row.names=F, col.names=F)

