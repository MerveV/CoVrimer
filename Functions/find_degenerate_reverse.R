
find_degenerate_reverse <-function(input, rev_bind_st, rev_bind_end) {
  x = 1:29903
  y = rep(0, 29903)
  nexts <- read.csv("next_MAFFT_21-0-920.csv", header = TRUE)
  data = cbind(x, y)
  data[nexts[, 1], 2] = nexts[, 7]
  rev_mut = data[rev_bind_st:rev_bind_end,]
  rev = strsplit(input, "")
  rv_copy = rev
  for (m in 1:length(rev[[1]])) {
    rv_copy[[1]][m] <-
      switch(
        rev[[1]][m],
        "g" = "c",
        "c" = "g",
        "a" = "t",
        "t" = "a",
        "y" = "Y",
        "r" = "R",
        "k" = "K",
        "m" = "M",
        "b" = "B",
        "v" = "V",
        "d" = "D",
        "h" = "H",
        "s" = "S",
        "w" = "W",
        "n" = "N"
      )
  }
  
  for (k in 1:length(rev[[1]])) {
    if (rev_mut[k, 2] != 0) {
      rv_copy[[1]][k] = nexts[which(nexts[, 1] == as.integer(rev_bind_st) + k - 1), 5]
    }
  }
  
  rv_copy2 <- rv_copy
  
  for (m in 1:length(rv_copy[[1]])) {
    rv_copy2[[1]][m] <-
      switch(
        rv_copy[[1]][m],
        "g" = "c",
        "c" = "g",
        "a" = "t",
        "t" = "a",
        "Y" = "R",
        "R" = "Y",
        "K" = "M",
        "M" = "K",
        "B" = "V",
        "V" = "B",
        "D" = "H",
        "H" = "D",
        "S" = "W",
        "W" = "S",
        "N" = "N"
      )
    
    
  }
  return(reverse(paste(unlist(rv_copy2), collapse = "")))
}