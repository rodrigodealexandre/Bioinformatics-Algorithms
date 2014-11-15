#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")
source("cyclopeptide_mass_spectrum.R")

linear_spectrum <- function(sequence){
    mass <- 0
    if(length(sequence) == 1){
        if(grepl("-", sequence)){
            sequence <- aaa_to_a(sequence)
        }
        split_seq <- strsplit(sequence, NULL)[[1]]
        n <- length(split_seq)
        count <- 1
        for(i in 1:n){
            for(c in 1:i){
                count <- count + 1
                mass[count] <- aa_mass(paste(split_seq[c:i], collapse=""))
            }
        }
        
    }
    else{
        spec <- sequence        
        n <- length(spec)
        count <- 1
        for(i in 1:n){
            for(c in 1:i){
                count <- count + 1
                mass[count] <- sum(as.numeric(paste(spec[c:i])))
            }
        }
    }
    return(sort(mass))
    
}

sequence <- "NQEL"
lin_spec <- linear_spectrum(sequence)
cat(lin_spec)