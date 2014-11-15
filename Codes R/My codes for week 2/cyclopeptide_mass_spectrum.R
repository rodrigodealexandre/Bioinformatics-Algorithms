setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")
pep_mass <- read.table("integer_mass_table.txt")

aa_seqence <- function(sequence){
    s <- strsplit(sequence, NULL)[[1]]
    n <- length(s)
    split_seq <- sequence
    
    if(n>1){
        s_extended <- c(s , s[1:(n-2)])
        
        for(x in 1:n){
            for(z in 1:(n-1)){
                split_seq <- c(split_seq, paste(s_extended[x:c(x+z-1)],collapse=""))
            }
        }
    }
    return(split_seq)
}

aa_mass <- function(sequence){
    count <- 0
    split_seq <- strsplit(sequence, NULL)[[1]]
    for(x in split_seq){
        count <- count + pep_mass[pep_mass==x,2]
    }
    return(count)
}

cyclopeptide_mass_spectrum <- function(sequence){
    count <- 0
    split_seq <- aa_seqence(sequence)
    for(s in split_seq){
        count <- c(count, aa_mass(s))
    }
    return(count[order(count)])
}