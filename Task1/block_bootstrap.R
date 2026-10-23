block_bootstrap <- function(data, B) {
  n <- nrow(X) # assuming that rows are observations and columns are variables
  
  # sample indices from {1, ..., n-B+1}
  # for the starting points of the blocks
  h = ceiling(n/B)
  i_s <- sample(seq(1, n-B+1), h, replace = TRUE)
  
  # generate the blocks based on the indices
  # and concatenate into a matrix
  replicate <- do.call(rbind, lapply(i_s, function(i) {data[i:(i+B-1), ]}))     # -1 because +1 in the sequence where indices are sampled from
  replicate <- replicate[1:n,]    # we keep n first rows - in case n is not divisible by B and we need to cut the last values off
  
  return(replicate)
}