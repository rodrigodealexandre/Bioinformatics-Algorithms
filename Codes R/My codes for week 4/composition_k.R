setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 4")

composition_k <- function(sequence, k){
    k_mers <- lapply(sequence,function(x){
        seq_loop_size <- nchar(x)-k+1
        
        kmers <- sapply(1:seq_loop_size, function(z){
            y <- z + k -1
            kmer <- substr(x=x, start=z, stop=y)
            return(kmer)
        })
        return(kmers)
    })
    
    uniq <- (unlist(k_mers))

    return(sort(uniq))
}

# 
# file <- readLines("data/string_com.txt")
# sequence <- file[3]
# k <- as.numeric(file[2])
# 
# 
# file <- readLines("data/dataset_197_3 (6).txt")
# sequence <- file[2]
# k <- as.numeric(file[1])
# 
# kmers <- composition_k(sequence, k)
# 
# write.table(kmers, file="output/composite_result.txt", quote=F,
#             sep="\n", row.names = F, col.names = F)
