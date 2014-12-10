#------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 4")
library("Biostrings")

string_spelled <- function(kmers){
    if(length(kmers)==1){
        kmers <- strsplit(kmers, " ")[[1]]
    }
    start_sequence <- kmers[1]
    len <- nchar(kmers[1])
    rest <- sapply(2:length(kmers),function(x){ 
        substr(kmers[x],len,len)
    })
    sequence <- paste(start_sequence, paste(rest, sep="", collapse=""), 
                      sep="", collapse="")
    return(sequence)
}

# kmers<- c("ACCGA",
#       "CCGAA",
#       "CGAAG",
#       "GAAGC",
#       "AAGCT")
# 
# data <- readLines("data/GenomePathString.txt")
# kmers <- data[2:(length(data)-2)]

# kmers <- readLines("data/dataset_198_3.txt")
# 
# cat(string_spelled(kmers))