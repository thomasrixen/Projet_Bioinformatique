# Load the 'ape' library for working with DNA sequences
library(ape)

# Specify the filename for the DNA sequence in FASTA format
sequence_filename <- "bacterial_sequence.fasta"

# Read the DNA sequence from the file
mat <- read.dna(sequence_filename, format = "fasta", as.character = TRUE)
genome_sequence <- as.vector(mat)
sequence_length <- length(sequence)  # Should be 'genome_sequence' instead of 'sequence'

# Define start and stop codons
start_codons <- c("ttg", "ctg", "ata", "att", "atc", "atg", "gtg")
stop_codons <- c("tga", "taa", "tag")

# Define a function to find the start positions of a given codon in a DNA sequence
find_frame_start <- function(sequence, frame_size = 3, increment, target_string) {
    indices <- c() 

    target_length <- nchar(target_string)
    sequence_length <- length(sequence)
    
    for (i in 1:(sequence_length - frame_size + 1)) {
        frame <- paste(sequence[i:(i + frame_size - 1)], collapse = "")

        if (frame == target_string) {
            indices <- c(indices, i)
            i <- i + target_length - 1
        } else {
            i <- i + increment - 1
        }
        
        if (i > sequence_length) {
            break
        }
    }

    return(indices)
}

# Define a function to find the indexes of multiple codons in a DNA sequence
get_codons_index <- function(sequence, codons) {
    list <- c()
    for (codon in codons) {
        list <- c(list, find_frame_start(sequence, 3, 1, codon))
    }
    sort_list <- sort(list)
    return(sort_list)
}

# Define a function to count the length of open reading frames (ORFs)
count_len <- function(start_indexes, stop_indexes, min_len, len_seq) {
    count <- 0
    start_index <- 1
    while (start_index < length(start_indexes)) {
        index <- start_indexes[start_index] + 3
        while (!is.element(index, stop_indexes) && index <= len_seq - 2) {
            index <- index + 3
        }
        if ((index + 3) - start_indexes[start_index] >= min_len) {
            count <- count + 1
        }
        while (start_indexes[start_index] < index + 3 && start_index <= length(start_indexes)) {
            start_index <- start_index + 1
        }
    }
    return(count)
}

# Get the indexes of all start codons in the DNA sequence
list_start <- get_codons_index(genome_sequence, start_codons)
cat("Indexes of all the start codons:", list_start, "\n")

# Get the indexes of all stop codons in the DNA sequence
list_stop <- get_codons_index(genome_sequence, stop_codons)
cat("Indexes of all the stop codons:", list_stop, "\n")

# Initialize a vector to store the counts of ORFs of different lengths
ORFs_result <- c(0, 0, 0, 0, 0, 0)
labels <- c("0", "10", "50", "100", "300", "500")
k_nuc <- c(0, 10, 50, 100, 300, 500)
names(ORFs_result) <- labels

# Count the number of ORFs of different lengths
for (k in k_nuc) {
    ORFs_result[as.character(k)] <- count_len(list_start, list_stop, k, length(genome_sequence))
}

# Display the results
ORFs_result

# Save the results as an RDS file
saveRDS(ORFs_result, file = "ORFs_result.rds")
