skew <- function(genome){
    output = as.vector(0)
    count = 0
    n = 2
    genome <- strsplit(genome, NULL)[[1]]
    for(x in genome){
        if(x == "G"){
            count = count + 1
        }
        else if(x == "C"){
            count = count - 1
        }
        
        output[n] <- count
        n = n+1
    }    
    return(cat(output))
}


library("stringr")

minskew <- function(genome){
    output = as.vector(0)
    count = 0
    n = 2
    genome <- strsplit(genome, NULL)[[1]]
    for(x in genome){
        if(x == "G"){
            count = count + 1
        }
        else if(x == "C"){
            count = count - 1
        }
    
        output[n] <- count
        n = n+1
    }    
    minSkew = min(output)
    minSkewIndexes <- which(output==minSkew, arr.ind=T)
    #note that the -1 is there because the answer was in python
    return(minSkewIndexes-1)
}
minskew(genome)


maxskew <- function(genome){
    output = as.vector(0)
    count = 0
    n = 2
    genome <- strsplit(genome, NULL)[[1]]
    for(x in genome){
        if(x == "G"){
            count = count + 1
        }
        else if(x == "C"){
            count = count - 1
        }
        
        output[n] <- count
        n = n+1
    }    
    maxSkew = max(output)
    maxSkewIndexes <- which(output==maxSkew, arr.ind=T)
    #note that the -1 is there because the answer was in python
    return(maxSkewIndexes-1)
}