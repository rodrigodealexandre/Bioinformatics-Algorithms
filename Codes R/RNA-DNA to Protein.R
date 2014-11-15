library("Biostrings")

sequence <- RNAString(readLines("dataset_96_5.txt"))
as.character(translate(sequence))



y <- c("ATGCATTGGACGTTAG") # Creates sample DNA sequence.
AAdf <- read.table(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/AA.txt", header=T,"\t") # Imports genetic code.
AAv <- AAdf[,2]; names(AAv) <- AAdf[,1] # Creates named vector with genetic code
y <- gsub("(...)", "\\1_", y) # Inserts "_" after each triplet.
y <- unlist(strsplit(y, "_")) # Splits on "_" and returns vectorized triplets.
y <- y[grep("^...$", y)] # Removes incomplete triplets.
AAv[y] # Translation into protein by name-based subsetting.