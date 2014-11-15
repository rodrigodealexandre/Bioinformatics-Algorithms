setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")


#-------------------------------------------------





#-------------------------------------------------
#load table
aa_nomeclature <- read.table("aa_nomeclature.txt")

#Convert a 3 letters to a 1 letter peptide
aaa_to_a <- function(sequence){
    sequence <- strsplit(sequence, "-")[[1]]
    sequence <- sapply(sequence, function(pep){    
        aa_nomeclature[grep(pep, aa_nomeclature[,2]),1]})
    sequence <- paste(sequence, collapse = "")
    return(sequence)
}


#Convert a 1 letter to a 3 letters peptide
a_to_aaa <- function(sequence){
    sequence <- strsplit(sequence, NULL)[[1]]
    sequence <- sapply(sequence, function(pep){    
        aa_nomeclature[grep(pep, aa_nomeclature[,1]),2]})
    sequence <- paste(sequence, collapse = "-")
    return(sequence)
}

#-------------------------------------------------
#load table
RNA_Codon <- read.table("RNA_Codon.txt")

#Count the number of different RNA sequences possible with a given AA sequence.
n_DNAStrings_from_AA <- function(sequence){
    #type aaa refers to when you have a 3 letters nomenclature for the AA spaced by "-".
    #type a refers to when you have a 1 letter nomenclature for the AA not spaced.
    
    if(grepl("-", sequence)){
        type = "aaa"
    }
    else{
        type = "a"
    }
    
    if(type=="aaa"){
        sequence <- aaa_to_a(sequence)
        sequence <- strsplit(sequence, NULL)[[1]]
        sequence <- sapply(sequence, function(pep){    
            length(grep(pep, RNA_Codon[,2]))    
        })
        return(prod(sequence))
        
    }
    
    if(type=="a"){
        sequence <- strsplit(sequence, NULL)[[1]]
        sequence <- sapply(sequence, function(pep){    
            length(grep(pep, RNA_Codon[,2]))    
        })
        return(prod(sequence))
    }
}


#-------------------------------------------------
#load table
RNA_Codon <- read.table("RNA_Codon.txt")

protein_to_RNA <- function(sequence){
    #type aaa refers to when you have a 3 letters nomenclature for the AA spaced by "-".
    #type a refers to when you have a 1 letter nomenclature for the AA not spaced.
    if(grepl("-", sequence)){
        type = "aaa"
    }
    else{
        type = "a"
    }
    
    if(type=="aaa"){
        sequence <- aaa_to_a(sequence)
        sequence <- strsplit(sequence, NULL)[[1]]
        if(length(sequence) <= 12){
            sequence <- expand.grid(lapply(sequence, function(pep){    
                RNA_Codon[grep(pep, RNA_Codon[,2]),1]
            }))
            sequence <- apply(sequence, 1,function(seq){paste(seq, collapse = "")})
            return(sequence)
        }
        else{
            return("sequence is to big")
        }
    }
    
    if(type=="a"){
        sequence <- strsplit(sequence, NULL)[[1]]
        if(length(sequence) <= 12){
            sequence <- expand.grid(lapply(sequence, function(pep){    
                RNA_Codon[grep(pep, RNA_Codon[,2]),1]
            }))
            sequence <- apply(sequence, 1,function(seq){paste(seq, collapse = "")})
            return(sequence)
        }
        else{
            return("sequence is to big")
        }
    }
}


#-------------------------------------------------
library("Biostrings")

search_AA_on_genome <- function(sequence, genome){
    genome <- DNAString(genome)
    sequence <- RNAStringSet(protein_to_RNA(sequence))
    r_sequence <- reverseComplement(sequence)
    sequence <- RNAStringSet(c(sequence, r_sequence))
    sequence <- as.character(DNAStringSet(sequence))
    
    x <- startIndex(matchPDict(DNAStringSet(sequence), genome))
    names(x) <- sequence
    
    sequence <- vapply(x,FUN.VALUE = 1, FUN=function(seq){
        
        if(is.integer(seq[1])){
            return(length(seq)[[1]])
        }
        else{
            return(0)
        }
    }
    )
    sequence <- sequence[sequence >0]
    sequence <- sequence[order(sequence,decreasing = T)]
    return(sequence)
    
}


#-------------------------------------------------
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


#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")

peptide_score <- function(sequence, spectrum){
    if(length(sequence) == 1){
        if(grepl("-", sequence)){
            sequence <- aaa_to_a(sequence)
        }
        sequence <- cyclopeptide_mass_spectrum(sequence)
    }
    
    if(length(spectrum) == 1){
        spectrum <- strsplit(spectrum, "\\s" )[[1]]
    }
    
    score <- sum(table(spectrum[spectrum %in% sequence]))
    return(score)
}

file <- readLines("data/dataset_4913_1 (1).txt")
sequence <- file[1]
spectrum <- file[2]
score <- peptide_score(sequence, spectrum)
cat(score)

#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")
pep_mass <- read.table("integer_mass_table.txt")
library(plyr)

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


#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")

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


#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")

linear_peptide_scoring <-function(sequence, spectrum){
    if(length(spectrum) == 1){
        spectrum <- strsplit(spectrum, "\\s" )[[1]]
    }
    print("linear_spec")
    sequence <- linear_spectrum(sequence)
    print("pep_score")
    score <- peptide_score(sequence, spectrum)
    return(score)
}

file <- readLines("data/dataset_4913_1 (1).txt")
sequence <- file[1]
spectrum <- file[2]
score <- linear_peptide_scoring(sequence, spectrum)
cat(score)


multiple_linear_peptide_scoring <-function(sequence, spectrum, n){
    if(length(sequence) == 1){
        sequence <- strsplit(sequence, "\\s" )[[1]]
    }
    sequence <- as.data.frame(t(sequence))
    scores <- sapply(sequence, function(x){linear_peptide_scoring(x, spectrum)})
    scores <- sort(scores, decreasing=T)
    return(scores[1:n])
}

file <- readLines("data/dataset_4913_3.txt")
sequence <- file[1]
spectrum <- file[2]
n <- file[3]

score <- multiple_linear_peptide_scoring(sequence, spectrum, n)
score
cat(names(score))

#-------------------------------------------------




#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Codes R/My codes for week 2")
pep_mass <- read.table("integer_mass_table.txt")
library(plyr)
linear_spectrum


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



#-------------------------------------------------

#-------------------------------------------------
