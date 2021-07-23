calculate_GC <- function(input) {
  inp = strsplit(input, "")
  w = length(which(inp[[1]] == "a"))
  x = length(which(inp[[1]] == "t"))
  y = length(which(inp[[1]] == "g"))
  z = length(which(inp[[1]] == "c"))
  GC = (y + z) / (y + z + w + x) * 100
  return(GC)
}