
find_degenerate <- function(input, forw_bind_st, forw_bind_end) {
  
  nexts <- read.csv("next_MAFFT_21-0-920.csv", header = TRUE)
  x = 1:29903
  y = rep(0, 29903)
  data = cbind(x, y)
  data[nexts[, 1], 2] = nexts[, 7]
  forw_mut = data[forw_bind_st:forw_bind_end,]
  forw = strsplit(input, "")
  for (k in 1:length(forw[[1]])) {
    if (forw_mut[k, 2] != 0) {
      ss = which(nexts[, 1] == as.integer(forw_bind_st) + k - 1)
      forw[[1]][k] = nexts[ss, 5]
    }
  }
  return(paste(unlist(forw), collapse = ""))
}