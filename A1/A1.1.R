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
  gc_frequencies <- frequency_in_sliding_window(sequence, c("G", "C"), size, shift)
  return(gc_frequencies)
}