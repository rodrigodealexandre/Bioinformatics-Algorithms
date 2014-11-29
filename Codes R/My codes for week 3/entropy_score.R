#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")
#Implement Log2 Score(Motifs)

library("Biostrings")

profile <- function(sequences){
    
    n <- length(sequences)
    sequences <- matrix(sequences)
    sequences <- lapply(sequences,function(y){paste(strsplit(y, "   ")[[1]], collapse="")})
    sequences <- DNAStringSet(toupper(unlist(sequences)))
    count<- consensusMatrix(sequences, baseOnly=TRUE)
    profile <- count/n
    return(profile)
}

entrophy <-  function(sequences){
    profile <- profile(sequences)
    log_profile <- -log2(profile)*profile
    log_profile[is.na(log_profile)] <- 0
    
    
    entrophy <-apply(log_profile,2,sum)
    entrophy
    
    return(entrophy)
}

# sequences <- c("T   C   G   G   G   G   g   T   T   T   t   t",          
#                "c   C   G   G   t   G   A   c   T   T   a   C",
#                "a   C   G   G   G   G   A   T   T   T   t   C",
#                "T   t   G   G   G   G   A   c   T   T   t   t",
#                "a   a   G   G   G   G   A   c   T   T   C   C",
#                "T   t   G   G   G   G   A   c   T   T   C   C",
#                "T   C   G   G   G   G   A   T   T   c   a   t",
#                "T   C   G   G   G   G   A   T   T   c   C   t",
#                "T   a   G   G   G   G   A   a   c   T   a   C",
#                "T   C   G   G   G   t   A   T   a   a   C   C")
# 
# round(sum(entrophy(sequences)),4)












##### OLD MANUAL ######

# library("tau")
# 
# sequences <- c("T   C   G   G   G   G   g   T   T   T   t   t",          
#                "c   C   G   G   t   G   A   c   T   T   a   C",
#                "a   C   G   G   G   G   A   T   T   T   t   C",
#                "T   t   G   G   G   G   A   c   T   T   t   t",
#                "a   a   G   G   G   G   A   c   T   T   C   C",
#                "T   t   G   G   G   G   A   c   T   T   C   C",
#                "T   C   G   G   G   G   A   T   T   c   a   t",
#                "T   C   G   G   G   G   A   T   T   c   C   t",
#                "T   a   G   G   G   G   A   a   c   T   a   C",
#                "T   C   G   G   G   t   A   T   a   a   C   C")
# n <- length(sequences)
# 
# bases <- c("A","C","G","T")
# sequences <- matrix(sequences)
# sequences <- lapply(sequences,function(y){strsplit(y, "   ")[[1]]})
# sequences <- matrix(unlist(sequences), ncol=length(sequences[[1]]), byrow=T)
# sequences <- toupper(sequences)
# 
# sequences <- apply(sequences,2,function(y){textcnt(y, n=1L, method="string", tolower = F)})
# 
# count <- list(NULL)
# c = 0
# for(x in sequences){
#     for(y in bases){
#         if(grepl(paste(names(x), collapse=" "),pattern = y)==F){
#             x[length(x)+1] <- 0
#             names(x)[length(x)] <- y
#         }
#     }
#     c <- c + 1
#     x <- x[order(names(x))]
#     count[[c]] <- x
# }
# 
# count <- matrix(unlist(count), ncol=length(count[[1]]), byrow=T)
# colnames(count) <- bases
# count <- t(count)
# 
# profile <- count/n
# 
# log_profile <- -log2(profile)*profile
# log_profile[is.na(log_profile)] <- 0
# 
# 
# entrophy <-apply(log_profile,2,sum)
# entrophy
# round(sum(entrophy),4)





