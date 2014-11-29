#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")

#Identifying the evening element


#Kmer expected number of occurrences in a DNA string
k_expected_n <- function(L, NBases, k, NOccurrence, NSequences) {
    return((choose(L - NOccurrence*(k - 1), NOccurrence)/NBases^(NOccurrence*k))*NSequences)
}


#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")

#Implement MOTIFENUMERATION.

library("gtools")
library("Biostrings")
library("stringr")

#counting number of mismatches
count_difference <- function(kmer1, kmer2, d){
    if(length(kmer1) == length(kmer2)){
        count <- length(kmer1)
        for(k in 1:length(kmer1)){
            if(kmer1[k] != kmer2[k]){
                count <- count - 1
            }
        }
        difference <- length(kmer1) - count
        if(difference <= d){
            return(TRUE)
        }
        else{
            return(FALSE)
        }
    }
    else(cat("K-mer1 and K-mer2 has different lenghts"))
}

#finding all kmers in each sequence
sequence_kmers <- function(sequence, k){
    k_mers <- lapply(sequence,function(x){
        seq_loop_size <- length(DNAString(x))-k+1
        
        kmers <- sapply(1:seq_loop_size, function(z){
            y <- z + k -1
            kmer <- substr(x=x, start=z, stop=y)
            return(kmer)
        })
        return(kmers)
    })
    return(k_mers)
}

#apply doesn't seems to impruve speed
motifenumerator <- function(DNAs, k, d){
    bases <- c("A","T","G","C")
    all_combinations <- permutations(4,k,bases, repeats.allowed=T)
    seq_kmers <- sequence_kmers(DNAs, k)
    
    motifs <- lapply(seq_kmers, function(x){
        is_true <- matrix(ncol=k)
        count <- 0
        lapply(x, function(y){
            y <- strsplit(y, split="")[[1]]
            count <<- count +1
            is_true <<- rbind(all_combinations[apply(all_combinations,1, count_difference, y, d),],is_true)
        })
        all_combinations <<- unique(is_true[!is.na(is_true[,1]),])
        return(all_combinations)
    })[[length(seq_kmers)]]
    
    return(apply(motifs,1,paste,collapse=""))
}


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



#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")

#Implement DistanceBetweenPatternAndStrings
library("Biostrings")

#Implement MEDIANSTRING
median_string <- function(k, sequences, max_d=3){
    if(length(sequences)==1){
        sequences <- strsplit(sequences, " ")[[1]]
    }
    
    all_kmers <- unlist(sequence_kmers(sequences, k))
    all_kmers <- unique(all_kmers)
    
    kmers <- sapply(all_kmers, function(x){sum(d_pattern_dna(x, sequences, max_d))})
    
    kmers <- sort(kmers)
    return(kmers)
}


#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")

#Implement Profile-most probable k-mer in Text

library("gtools")
library("Biostrings")
source("motifenumerator.R")

most_prob_kmer <- function(sequence, profile, k){
    bases <- c("A","C", "G", "T")
    k <- as.numeric(k)
    profile <- matrix(as.numeric(unlist(lapply(profile, function(x){strsplit(x, " ")[[1]]}))), nrow=4, byrow=T)
    rownames(profile) <- bases
    
    #all_kmers <- permutations(4,k,bases, repeats.allowed=T)
    all_kmers <- sequence_kmers(sequence,k)
    all_kmers <- matrix(all_kmers, nrow=length(all_kmers))
    all_kmers <- t(apply(all_kmers,1,function(x){strsplit(x, NULL)[[1]]}))
    
    all_k_profile <- t(apply(all_kmers, 1, function(x){
        sapply(1:k, function(y){
            profile[profile=x[y],y]
        })
    }))
    all_k_profile <- data.frame(all_k_profile)
    all_k_profile[,k+1] <- apply(all_k_profile,1,prod)
    
    sorted_profile <- order(all_k_profile[,k+1], decreasing=T)
    sorted_kmers_names <- all_kmers[sorted_profile,]
    sorted_profile <- all_k_profile[sorted_profile,]
    sorted_kmers <- sorted_profile[,k+1]
    names(sorted_kmers) <- apply(sorted_kmers_names,1,paste,collapse="")
    return(head(sorted_kmers))
}




#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")






#-------------------------------------------------
setwd("D:/Dropbox/Courses/Coursera Courses/Bioinformatics Algorithms (Part 1)/Bioinformatics-Algorithms/Codes R/My codes for week 3")
