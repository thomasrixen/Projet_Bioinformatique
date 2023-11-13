# Load the rentrez package
library(rentrez)

# Function to retrieve sequences from GenBank
retrieve_sequences <- function(accession_numbers) {
  sequences <- list()
  count = 0
  for (accession in accession_numbers) {
    # Remove any potential whitespace or new line characters
    accession <- trimws(accession)
    # Fetch the sequence
    sequence_data <- entrez_fetch(db="nucleotide", id=accession, rettype="fasta")
    sequences[[accession]] <- sequence_data
    count <- count + 1
    cat(count,"%\n")
  }
  return(sequences)
}

# Read the accession numbers from a file (replace 'path_to_your_file.txt' with the actual file path)
accession_numbers <- readLines('AccessionNumbers.txt')

# Retrieve sequences
sequences <- retrieve_sequences(accession_numbers)

# Optionally, save the sequences to a file
#writeLines(unlist(sequences), 'genbank_sequences.fasta')


#======================================================================================================



# Function to filter and process the sequences
filter_sequences <- function(sequences) {
  # Extracting the accession number and sequence
  parsed_sequences <- lapply(names(sequences), function(name) {
    sequence_data <- unlist(strsplit(sequences[[name]], "\n", fixed = TRUE))
    sequence_header <- sequence_data[1]
    sequence <- paste(sequence_data[-1], collapse = "")
    list(accession_number = name, sequence = sequence)
  })

  # Convert the list to a dataframe for easier processing
  sequence_df <- do.call(rbind.data.frame, parsed_sequences)

  # Sort by accession number
  sequence_df <- sequence_df[order(sequence_df$accession_number), ]

  # Keep sequences with only A, T, C, G
  sequence_df <- subset(sequence_df, grepl("^[ATCG]+$", sequence))

  # Remove duplicates (keeping the first occurrence)
  sequence_df <- sequence_df[!duplicated(sequence_df$sequence), ]

  return(sequence_df)
}

# Example usage with previously retrieved sequences
# Assume 'sequences' is the list of sequences retrieved in the previous step
filtered_sequences <- filter_sequences(sequences)

# The 'filtered_sequences' dataframe now contains the processed sequences
# Function to report the number of filtered sequences
count_filtered_sequences <- function(filtered_sequences) {
  return(nrow(filtered_sequences))
}

# Function to report the average length of the filtered sequences
average_length_of_sequences <- function(filtered_sequences) {
  if (nrow(filtered_sequences) == 0) {
    return(0)
  }
  return(mean(nchar(filtered_sequences$sequence)))
}

# Usage Example
# Assuming 'filtered_sequences' is the dataframe obtained from the 'filter_sequences' function
num_sequences <- count_filtered_sequences(filtered_sequences)
avg_length <- average_length_of_sequences(filtered_sequences)

print(paste("Number of Filtered Sequences:", num_sequences))
print(paste("Average Length of Sequences:", avg_length))

