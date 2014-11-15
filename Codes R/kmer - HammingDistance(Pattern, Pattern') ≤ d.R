library("Biostrings")

cat(start(matchPattern(pattern, genome,max.mismatch=d))-1)

genome = DNAString("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT")
pattern = DNAString("ATTCTGGA")
d = 3

