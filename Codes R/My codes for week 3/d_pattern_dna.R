#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")

#Implement DistanceBetweenPatternAndStrings
library("Biostrings")

d_pattern_dna <- function(pattern, sequence, max_d=3){
    if(length(sequence)==1){
        sequence <- strsplit(sequence, " ")[[1]]
    }
    pattern_number <- length(strsplit(pattern, NULL)[[1]])
    
    pattern <- DNAString(pattern)
    matched <- NULL
    n <- 0
    for(x in sequence){
        count <- 0
        n <- n + 1
        sequence_string <- DNAString(x)
        match_pattern <- NULL
        while(count <= max_d){
            match_pattern <- matchPattern(pattern, sequence_string, max.mismatch=count)
            count <- count + 1
            if(length(match_pattern) > 0){
                break
            }
        }
        if(length(match_pattern)> 0){
            for (i in 1:length(match_pattern)) {
                possibleError <- tryCatch(
                    as.character(match_pattern[i]),
                    error=function(e) e
                )
                
                if(!inherits(possibleError, "error")){
                    matched[n] <- as.character(match_pattern[i])
                } 
            }
        }
        else{
            matched[n] <- paste(rep("N",pattern_number),collapse="")
        }
        
    }
    matched[is.na(matched)] <- paste(rep("N",pattern_number),collapse="")
    hamming_distance <- sapply(matched, function(x){stringDist(c(as.character(pattern), x))})
    return(hamming_distance)
}


# 
# pattern <- "AAA"
# sequence <- "TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT"
# sum(d_pattern_dna(pattern, sequence))
# 
# 
# file <- readLines("data/dataset_5164_1.txt")
# pattern <- file[1]
# sequence <- file[2]
# sum(d_pattern_dna(pattern, sequence))
