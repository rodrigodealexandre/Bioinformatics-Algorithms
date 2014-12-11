#------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 4")
source("composition_k.R")

de_bruijn_graph_k <- function(kmers){
    if(length(kmers)==1){
        kmers <- strsplit(kmers, "\\s+")[[1]]
    }
    kmers <- sort(kmers)
    len <- nchar(kmers[1])
    first_k <- unique(substr(kmers,start=1,stop=(len-1)))
    result <- unlist(lapply(first_k, function(x){
        paste(x, paste0(substr(kmers[grepl(paste0("^", x),kmers)], 
                               start=2, stop=len), collapse=","), sep=" -> ")}))
    
    return(sort(result))    
}

# kmers <- c("GAGG", "CAGG", "GGGG", "GGGA", "CAGG", "AGGG", "GGAG")
# 
# data <- readLines("data/De_Bruijn_Graph_from_kmer.txt")
# kmers <- data[2:1122]
# 
# kmers <- readLines("data/dataset_200_7.txt")
# 
# result <- de_bruijn_graph_k(kmers)
# write.table(result, file="output/debruijngraph_k_result.txt", quote=F, 
#             sep="\n", row.names=F, col.names=F)
