library(ape)

#=================================================#
#                   QUESTION 1                    #
#=================================================#
mat <- read.dna("NC_005337.fasta", format="fasta", as.character=TRUE)
seq <- as.vector(mat)

sequence_length = length(seq)
print("sequence lenght is equal to :")
sequence_length

#=================================================#
#                   QUESTION 2                    #
#=================================================#

dimer_matrix <- matrix(c(
    0.0975, 0.0589, 0.0489, 0.0863,
    0.0660, 0.0422, 0.0503, 0.0506,
    0.0626, 0.0447, 0.0430, 0.0588,
    0.0656, 0.0634, 0.0668, 0.0945
), nrow = 4, byrow = TRUE)


rownames(dimer_matrix) <- c("A*", "C*", "G*", "T*")
colnames(dimer_matrix) <- c("*A", "*C", "*G", "*T")

saveRDS(dimer_matrix, "dimer_frequencies.rds")


#=================================================#
#                   QUESTION 3                    #
#=================================================#



observed_dimer_frequencies <- matrix(0, nrow = 4, ncol = 4)
rownames(observed_dimer_frequencies) <- c("A*", "C*", "G*", "T*")
colnames(observed_dimer_frequencies) <- c("*A", "*C", "*G", "*T")


count_dimer <- function(nuc1, nuc2, seq){
    index <- 1
    count <- 0
    len_seq <- length(seq)

    while(index < len_seq){
        if (index + 1 > len_seq) {
           break
        }
        if (!is.na(seq[index]) && !is.na(seq[index + 1]) && seq[index] == nuc1 && seq[index + 1] == nuc2) {
            count <- count + 1
        }
        index <- index + 1

    }
    count
}


for (colonne in colnames(observed_dimer_frequencies)) {
    for (ligne in rownames(observed_dimer_frequencies)) {
      observed_dimer_frequencies[ligne, colonne] = count_dimer(get_lower_letter(ligne), get_lower_letter(colonne), seq)/len
      #print(observed_dimer_frequencies)
    }
}


get_lower_letter <- function(input_string) {
  if (grepl("[A-Z]\\*", input_string)) {
    lower_letter <- tolower(sub(".*([A-Z])\\*", "\\1", input_string))
    return(lower_letter)
  } else if (grepl("\\*[A-Z]", input_string)) {
    lower_letter <- tolower(sub("\\*([A-Z]).*", "\\1", input_string))
    return(lower_letter)
  } else {
    return(input_string)
  }
}

print("This is the observed dimer frequencies in the sequence :")
observed_dimer_frequencies

get_top_labels <- function(matrix, get_largest = TRUE) {
  values <- as.vector(matrix[ , -1])
  labels <- matrix[ , 1]

  sorted_indices <- ifelse(get_largest, order(values, decreasing = TRUE), order(values))

  top_indices <- sorted_indices[1:3]
  top_labels <- labels[top_indices]

  return(top_labels)
}

saveRDS(observed_dimer_frequencies, "observed_dimer_frequencies.rds")