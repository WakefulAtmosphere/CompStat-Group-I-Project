block_bootstrap <- function(data, B) {
  n <- nrow(data) # olettaen ett? rivit on havaintoja ja sarakkeet eri muuttujia
  
  # assuming that the length of the data is divisible by B
  # sample indices from h = {1, ..., n-B+1}
  # for the starting points of the blocks
  i_s <- sample(seq(1, n-B+1), n/B)
  
  # generate the blocks based on the indices
  # and concatenate into a matrix/data frame
  replicate <- do.call(rbind, lapply(i_s, function(i) {data[i:(i+B), ]}))
  
  # do I repeat that x times??
  # or for all possible mixings??
  
  return(replicate)
}