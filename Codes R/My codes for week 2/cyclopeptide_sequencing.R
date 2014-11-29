#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 2")
pep_mass <- read.table("integer_mass_table.txt")
library(plyr)
source("cyclopeptide_mass_spectrum.R")

cyclopeptide_sequencing <-function(spectrum){
    if(length(spectrum) == 1){
        spectrum <- strsplit(spectrum, "\\s" )[[1]]
    }
    if(length(spectrum) == 1){
        
        spectrum <- cyclopeptide_mass_spectrum(spectrum)
    }
    
    if(length(spectrum) > 2){
        in_pep_mass <- spectrum[which(spectrum %in% pep_mass[,2])]
        in_pep_mass <- as.numeric(in_pep_mass)
        
        unique_values <- unique(in_pep_mass)
        compose <- unique(in_pep_mass)
        
        loopCount = 2
        compose <- lapply(1:length(compose), function(x){expand.grid(compose[x],unique_values)})
        compose <- ldply(compose)
        compose <- compose[sapply(1:length(compose[,1]), function(x){
            length(in_pep_mass[which(in_pep_mass %in% compose[x,])])>=loopCount & 
                length(which(spectrum %in% sum(compose[x,])))>0}),]
        if(length(in_pep_mass)>2){
            for(i in 3:length(in_pep_mass)){
                loopCount = loopCount +1
                compose <- lapply(1:length(compose[,1]), function(x){
                    compose <- do.call('expand.grid', c(as.data.frame(compose)[x,],as.data.frame(unique_values)))})
                compose <- ldply(compose)
                names(compose)[i] <- paste("Var",i,sep="")
                compose <- compose[sapply(1:length(compose[,1]), function(x){
                    is_true <- length(in_pep_mass[which(in_pep_mass %in% compose[x,])])>=loopCount
                    for(z in 1:length(compose[x,])){
                        l <- length(compose[x,])
                        is_true <- is_true & length(which(spectrum %in% sum(compose[x,z:l])))>0
                    }
                    return(is_true)
                }),]}}
        compose <- as.vector(apply(compose, 1, paste, collapse="-"))
        return(compose)
    }
    else{
        return(cat("Invalid imput"))
    }
}


#spectrum <- readLines("data/cycloseq_data.txt")
spectrum = "0 113 128 186 241 299 314 427"
cyclo <- cyclopeptide_sequencing(spectrum)
cat(cyclo)
