motif_update <- function(s2_kmers, profile){
    motifs <- sapply(s2_kmers, function(z){
        zz <- strsplit(z, NULL)[[1]]
        prod <- sapply(1:k, function(w){
            profile[profile=zz[w],w]
        })
        prod(prod)
    })
    motifs <- sort(motifs[motifs>0],decreasing = T)
    if(length(motifs)==0){
        return(best_motifs[y])
    }
    else{
        return(names(motifs)[1])
    }
}