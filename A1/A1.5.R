library(ape)

sequence_filename <- "bacterial_sequence.fasta"

mat <- read.dna(sequence_filename, format = "fasta", as.character = TRUE)
genome_sequence <- as.vector(mat)
sequence_length <- length(sequence)

start_codons <- c("ttg", "ctg", "ata", "att", "atc", "atg", "gtg")
stop_codons <- c("tga", "taa", "tag")

#Argument:
# - sequence: genome input
# - frame_size: taille du codon (mit à 3 par défaut)
# - increment: Permet de bouger la frame (on utilisera 1 pour trouver tout les start et stop et 3 une fois qu'on veut trouver le premier stop qui suit un start)
# - target_string: correspond à une string de taille 3, des codons start ou stop

#Return
#Retourne une list avec tout les indexes de la sequence où commence un codon (target_string)

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

#Exemple d'utilisation

test_sequence <- c("t", "t", "c", "t", "g", "t", "a", "c", "t", "g", "t", "t", "a", "g", "a", "t")
# target_string <- "ctg"
# frame_size <- 3
# increment <- 1

# result <- find_frame_start(test_sequence, frame_size, increment, target_string)
# print(result) #return 3 et 8

#

#Argument:
# - sequence: genome input
# - codons: list de codons dont on veut trouver les indexes

#Return
#Une liste triée avec tout les indexes de codons start ou stop

get_codons_index <- function(sequence, codons){
    list <- c()
    for (codon in codons) {
       list <- c(list, find_frame_start(sequence, 3, 1, codon))
    }
    sort_list <- sort(list)
    return(sort_list)
}

count_len <- function(start_indexes, stop_indexes, min_len, len_seq){
    count <- 0
    start_index <- 1
    while(start_index < length(start_indexes)){
        index <- start_indexes[start_index] + 3
        while(!is.element(index, stop_indexes) && index <= len_seq-2){
            index <- index + 3
        }
        if ((index + 3)-start_indexes[start_index] >= min_len) {
           count <- count + 1
        # uncomment if you want to have details about all the sequence   
        #    cat("The codon nr :", count, "\n")
        #    cat("start_index = ", start_indexes[start_index], "\n")
        #    cat("end_index = ", index, "\n")
        #    cat("longueur = ", (index + 3)-start_indexes[start_index], "\n")
        }
        while(start_indexes[start_index] < index + 3 && start_index <= length(start_indexes)){
            start_index <- start_index + 1
        }
    }
    return(count)
}

list_start <- get_codons_index(genome_sequence, start_codons)
cat("Indexes of all the start codons:", list_start, "\n")

list_stop <- get_codons_index(genome_sequence, stop_codons)
cat("Indexes of all the stop codons:", list_stop, "\n")

ORFs_result <- c(0, 0, 0, 0, 0, 0)
labels <- c("0", "10", "50", "100", "300", "500")
k_nuc = c(0, 10, 50, 100, 300, 500)
names(ORFs_result) <- labels

for (k in k_nuc) {
   ORFs_result[as.character(k)] <- count_len(list_start, list_stop, k, length(genome_sequence))
}

ORFs_result

saveRDS(ORFs_result, file = "ORFs_result.rds")
