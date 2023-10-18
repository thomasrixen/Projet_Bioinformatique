seq_length <- function(fasta_file){
  list <- read.dna(fasta_file,format = "fasta", as.character=TRUE)
  return(length(list))
}
nb_sliding_windows <- function(l, size, shift){
  high<-size
  ret=1
  while(high<l){
    high<-high+shift
    ret<-ret+1
  }
  return(ret)
}
frequency_in_sliding_window <- function(sequence, ntd, size, shift) {
  nb_windows <- nb_sliding_windows(length(sequence), size, shift)
  frequencies <- numeric(nb_windows)
  
  for (i in 1:nb_windows) {
    start <- 1 + (i - 1) * shift
    end <- start + size - 1
    window <- sequence[start:min(end, length(sequence))]
    
    # Calculate the count of 'ntd' in the window without using "window == ntd"
    ntd_count <- sum(window %in% ntd)
    
    # Calculate the relative frequency
    frequencies[i] <- ntd_count / length(window)
  }
  
  return(frequencies)
}
gc_content <- function(sequence, size, shift) {
  gc_frequencies <- frequency_in_sliding_window(sequence, c("g", "c"), size, shift)
  return(gc_frequencies)
}
at_content <- function(sequence, size, shift) {
  at_frequencies <- frequency_in_sliding_window(sequence, c("a", "t"), size, shift)
  return(at_frequencies)
}
require(ape)
data<-read.dna("AF086833.2.fasta", format="fasta", as.character = TRUE)
seq <- as.vector(data)
print(length(seq))
gc_data<- gc_content(seq, 500, 250)
at_data<- at_content(seq, 500, 250)
print(gc_data)
print(at_data)
x <- seq(1,length(seq),length(seq)/length(gc_data))
df <- data.frame(x, gc_data, at_data)

require(ggplot2)

g <- ggplot(df,aes(x=x))
g <-  g + geom_line(aes(y=gc_data), colour="red", linetype="dashed")
g <- g + geom_line(aes(y=at_data), colour='green')
g