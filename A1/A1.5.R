library(seqinr)

sequence_filename <- "bacterial_sequence.fasta"

mat <- read.dna(sequence_filename, format = "fasta", as.character = TRUE)
genome_sequence <- as.vector(mat)
sequence_length <- length(sequence)

start_codons <- c("ttg", "ctg", "ata", "att", "atc", "atg", "gtg")
stop_codons <- c("tga", "taa", "tag")

find_frame_start <- function(sequence, frame_size, increment, target_string) {
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

# test_sequence <- c("t", "t", "c", "t", "g", "t", "a", "c", "t", "g")
# target_string <- "ctg"
# frame_size <- 3
# increment <- 1

# result <- find_frame_start(test_sequence, frame_size, increment, target_string)
# print(result)


get_codons_index <- function(sequence, codons){
    list <- c()
    for (codon in codons) {
       list <- c(list, find_frame_start(sequence, 3, 1, codon))
    }
    sort_list <- sort(list)
    return(sort_list)
}


# frame_size <- 3
# increment <- 1
# target_string <- "ctg"

# start_index <- find_frame_start(genome_sequence, frame_size, increment, target_string)
# cat("Indexes of the start of the frame:", start_index, "\n")

list_start <- get_codons_index(genome_sequence, start_codons)
cat("Indexes of all the start codons:", list_start, "\n")

list_stop <- get_codons_index(genome_sequence, stop_codons)
cat("Indexes of all the stop codons:", list_stop, "\n")

