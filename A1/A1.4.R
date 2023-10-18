library(ape)

#=================================================#
#                   QUESTION 1                    #
#=================================================#
mat <- read.dna("NC_005337.fasta", format="fasta", as.character=TRUE)
seq <- as.vector(mat)
len = length(seq)


odds_ratios <- matrix(0, nrow = 4, ncol = 4)
rownames(odds_ratios) <- c("A*", "C*", "G*", "T*")
colnames(odds_ratios) <- c("*A", "*C", "*G", "*T")

count_nuc <- function(sequence, nuc) {
    position <- which(sequence == nuc)
    num_nuc <- length(position)
    return(num_nuc/len)
}

count_a <- count_nuc(seq, "a")
count_t <- count_nuc(seq, "t")
count_c <- count_nuc(seq, "c")
count_g <- count_nuc(seq, "g")

background_prob <- c("a" = count_a, "t" = count_t, "c" = count_c, "g" = count_g)
background_prob



observed_dimer_frequencies <- readRDS("observed_dimer_frequencies.rds")
observed_dimer_frequencies

compute_odds_ratios <- function(nuc1, nuc2, background, odf) {
    up_nuc1 <- uper_star_char(nuc1, TRUE)
    up_nuc2 <- uper_star_char(nuc2, FALSE)
    observed_frequency <- odf[up_nuc1, up_nuc2]

    expected_frequency_nuc1 <- background[nuc1]
    expected_frequency_nuc2 <- background[nuc2]

    return (observed_frequency/(expected_frequency_nuc1*expected_frequency_nuc2))
}

uper_star_char <- function(char, add_star) {
  if (add_star) {
    # If add_star is TRUE, add "*" to the right of the uppercase character without a space
    result <- paste0(toupper(char), "*")
  } else {
    # If add_star is FALSE, add "*" to the left of the uppercase character without a space
    result <- paste0("*", toupper(char))
  }
  return(result)
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

for (colonne in colnames(odds_ratios)) {
    for (ligne in rownames(odds_ratios)) {
        odds_ratios[ligne, colonne] = compute_odds_ratios(get_lower_letter(ligne), get_lower_letter(colonne), background_prob, observed_dimer_frequencies)
    }
}



odds_ratios

saveRDS(odds_ratios, file = "odds_ratios.rds")