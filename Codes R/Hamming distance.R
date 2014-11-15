seq1= "CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT"
seq2= "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"


library("stringdist")
stringdist(seq1, seq2, method="h")


library("Biostrings")
compare <- DNAStringSet(c(seq1, seq2))
stringDist(compare, method = "hamming")
