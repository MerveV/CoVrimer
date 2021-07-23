
find_degenerate_reverse2 <-function(input, rev_bind_st, rev_bind_end) {
  source("find_degenerate.R")
  x = 1:29903
  y = rep(0, 29903)
  nexts <- read.csv("next_MAFFT_21-0-920.csv", header = TRUE)
  data = cbind(x, y)
  data[nexts[, 1], 2] = nexts[, 7]
  
  
  tcom=reverseComplement(DNAString(input))
  tcom= as.character(tcom)
  deg=find_degenerate(tolower(tcom), rev_bind_st, rev_bind_end)
 # deg=find_degenerate(tolower(tcom),562,579)

  rev = strsplit(deg, "")
  rv_copy=rev
  for (m in 1:length(rev[[1]])) {
    rv_copy[[1]][m] <-
      switch(
        rev[[1]][m],
        "g" = "g",
        "c" = "c",
        "a" = "a",
        "t" = "t",
        "Y" = "Y",
        "R" = "R",
        "K" = "K",
        "M" = "M",
        "B" = "B",
        "V" = "V",
        "D" = "D",
        "H" = "H",
        "S" = "S",
        "W" = "W",
        "N" = "N"
      )
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
        "W" = "W",
        "S" = "S",
        "N" = "N"
      )
    
  }
  
  r=reverse(paste(unlist(rv_copy2), collapse = ""))
      
  return(r)
}