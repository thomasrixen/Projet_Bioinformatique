library(Biostrings)

# Sequences
seq1 <- DNAString("ATTCCGA")
seq2 <- DNAString("TATCG")

substitution_matrix <- matrix(c(2, -1, -1, -1,
                                -1,  2, -1, -1,
                                -1, -1,  2, -1,
                                -1, -1, -1,  2), 
                              nrow = 4, 
                              dimnames = list(c("A", "C", "G", "T"), 
                                              c("A", "C", "G", "T")))

# Perform global alignment
alignment <- pairwiseAlignment(pattern=seq1, subject=seq2, 
                               substitutionMatrix=substitution_matrix,
                               gapOpening=0, gapExtension=4, 
                               type="overlap")


# Extract aligned sequences
aligned_pattern <- as.character(alignedPattern(alignment))
aligned_subject <- as.character(alignedSubject(alignment))

# Initialize mismatch count
mismatch_count <- 0
gaps_count <- 0

# Iterate through the aligned sequences and count mismatches
for (i in 1:nchar(aligned_pattern)) {
  if (substr(aligned_pattern, i, i) != substr(aligned_subject, i, i) &&
      substr(aligned_pattern, i, i) != "-" &&
      substr(aligned_subject, i, i) != "-") {
    mismatch_count <- mismatch_count + 1
  }
}

# Iterate through the aligned sequences and count mismatches
for (i in 1:nchar(aligned_pattern)) {
  if (substr(aligned_pattern, i, i) == "-" ||
      substr(aligned_subject, i, i) == "-") {
    gaps_count <- gaps_count + 1
  }
}

alignment_score <- score(alignment)
alignment_match <- nmatch(alignment)
# Output the mismatch count
cat("Number of Score:", alignment_score, "\n")
cat("Number of Match:", alignment_match, "\n")
cat("Number of Mismatches:", mismatch_count, "\n")
cat("Number of Gaps:", gaps_count, "\n")

print(c(alignment_score, alignment_match, mismatch_count, gaps_count))