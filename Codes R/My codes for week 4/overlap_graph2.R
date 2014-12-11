#------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 4")
library("Biostrings")

overlap_graph <- function(kmers, sequence=""){
    if(length(kmers)==1){
        kmers <- strsplit(kmers, " ")[[1]]
    }
    kmers <- sort(kmers)
    len <- nchar(kmers[1])
    result <- sapply(1:length(kmers),function(x){ 
        the_k <- kmers[x]
        n_next_k <- substr(kmers[x],start=2,stop=len)
        next_k <- kmers[grepl(paste0("^",n_next_k),kmers)]
        
        if(length(next_k) != 0){
            if(length(next_k) > 1){
                if(sequence != ""){
                    positions <- list()
                    for(n in 1:length(next_k)){
                        next_k_position <- gregexpr2(next_k[n],sequence)[[1]]
                        the_k_position <- gregexpr2(the_k,sequence)[[1]]
                        test <- list()
                        counting <- 0
                        for(y in 1:length(the_k_position)){
                            counting <- counting +1
                            if(length(next_k_position[next_k_position>the_k_position[y]])>0){
                                test[[counting]] <- next_k_position[next_k_position>the_k_position[y]]
                            }
                        }
                        next_k_position<-test
                        next_k_position <- next_k_position[!is.na(next_k_position)]
                        if(length(next_k_position)>0){
                            positions[[n]] <- paste0(next_k[n],",",unlist(next_k_position))
                        }
                    }
                    positions <- unique(unlist(positions))
                    real_kmers <- substr(positions,start=1,stop=len)
                    real_pos <- as.numeric(substr(positions,start=len+2,stop=nchar(positions)))
                    
                    count <- 0
                    r_kmers <- NULL
                    for(i in 1:length(real_pos)){
                        for(o in 1:length(the_k_position)){
                            if(real_pos[i]-the_k_position[o]==1){
                                count <- count + 1
                                r_kmers[count] <- real_kmers[i]
                                }
                            }
                    }
                    
                    next_k <- r_kmers[!is.na(r_kmers)]
                }
                next_k <- paste0(next_k, collapse=",")
            }
            k_vector <- paste(the_k, "->",next_k)
            
        }
        else{
            k_vector <- NA
        }
        k_vector
    })
    result <- result[!is.na(result)]
    return(result)
}


# data <- readLines("data/overlap_graph_1.txt")
# kmers <- data[2:9977]


# kmers <- readLines("data/dataset_198_9.txt")
# result <- overlap_graph(kmers)

# write.table(result, file="output/overlapgraph_result.txt", quote=F, 
#             sep="\n", row.names=F, col.names=F)


