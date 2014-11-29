#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 2")
pep_mass <- read.table("integer_mass_table.txt")
library(plyr)
source("linear_spectrum.R")


leaderboard_cyclopeptide_sequencing <-function(spectrum, n){
    if(length(spectrum) == 1){
        spectrum <- strsplit(spectrum, "\\s" )[[1]]
    }
    if(length(spectrum) == 1){
        
        spectrum <- cyclopeptide_mass_spectrum(spectrum)
    }
    spectrum <- sort(as.numeric(spectrum))
    max_mass <- max(as.numeric(spectrum))
    
    if(length(spectrum) > 2){    
        unique_values <- unique(pep_mass[,2])
        compose <- unique(pep_mass[,2])
        score <- 0
        spec <- NULL
        count <- 1
        compose <- as.data.frame(compose)
        names(compose) <- paste("Var", count, sep="")
        
        while(TRUE){
            print(paste("Loop number", count, sep=" ")
                  count <- count + 1
                  unique_values <- as.data.frame(unique_values)
                  names(unique_values) <- paste("Var", count, sep="")
                  compose <- as.data.frame(compose)
                  compose <- lapply(1:length(compose[,1]), function(x){
                      do.call('expand.grid', c(compose[x,],unique_values))})
                  compose <- ldply(compose)
                  compose <- compose[apply(compose, 1, sum) <= max_mass,]
                  
                  if(length(compose[,1])==0){
                      break
                  }
                  print("scoring")
                  scores <- apply(compose, 1, function(x){
                      lin <- linear_spectrum(x)
                      sum(table(lin[lin %in% spectrum]))})
                  compose <- compose[order(scores, decreasing=TRUE)[1:n],]
                  scores <- scores[order(scores, decreasing=TRUE)[1:n]]
                  scores <- scores[!is.na(scores)]
                  compose <- compose[!is.na(compose[,1]),]
                  for(z in 1:length(compose[,1])){
                      if(sum(compose[z,])==max_mass){
                          if(scores[z] > score){
                              score <- scores[z]
                              spec <- compose[z,]
                          }
                      }
                  }
                  test <- compose
                  test[,length(test)+1] <- scores 
        }
        result <- paste(spec, collapse="-")
        result <- paste(result, "score =", score)
        return(result)
    }
}


file <- readLines("data/leaderboard.txt")
spectrum <- file[3]
n <- file[2]
leaderboard_cyclopeptide_sequencing(spectrum, n)


file <- readLines("data/dataset_102_7 (1).txt")
spectrum <- file[2]
n <- file[1]
leaderboard_cyclopeptide_sequencing(spectrum, n)
