#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 2")
pep_mass <- read.table("integer_mass_table.txt")
source("mass_spectrum.R")


cyclo_order <-function(spectrum){
    if(length(spectrum) == 1){
        spectrum <- strsplit(spectrum, "\\s" )[[1]]
    }
    if(length(spectrum) == 1){
        
        spectrum <- linear_mass_spectrum(spectrum)
    }
    
    
    if(length(spectrum) > 2){
        
        spectrum <- spectrum[which(spectrum %in% pep_mass[,2])]
        spectrum <- as.numeric(spectrum)
        
        cyclo_order <- list(NULL)
        n <- 1
        for(i in 1:length(spectrum)){
            cyclo_order[[n]] <- spectrum[c(i:length(spectrum),1:i-1)]
            cyclo_order[[n+1]] <- rev(spectrum[c(i:length(spectrum),1:i-1)])
            n <- n + 2
        }
        cyclo_order <- unlist(lapply(cyclo_order, paste, collapse="-"))
        return(cyclo_order)
    }
    else{
        return(cat("Invalid imput"))
    }
}
