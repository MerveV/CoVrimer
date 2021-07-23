find_Gc <- function(input) {
  primer <- strsplit(as.character(input), "")
  maxprimer <- primer
  minprimer <- primer
  for (m in 1:length(primer[[1]])) {
    # for max primer GC
    
    if (primer[[1]][m] %in% c("g","G",  "R", "S", "K", "B", "D")) {
      maxprimer[[1]][m] = "g"
    } else if (primer[[1]][m] %in% c("c", "C", "Y", "M", "V", "H")) {
      maxprimer[[1]][m] = "g"
    } else if (primer[[1]][m] %in% c("a","A", "W")) {
      maxprimer[[1]][m] = "a"
    } else {
      maxprimer[[1]][m] = "t"
    }
    
    # for min primer GC
    
    if (primer[[1]][m] %in% c("g", "G","S")) {
      minprimer[[1]][m] = "g"
    } else if (primer[[1]][m] %in% c("t","T", "Y", "K", "B", "D", "H")) {
      minprimer[[1]][m] = "a"
    } else if (primer[[1]][m] %in% c("a","A", "R", "M", "V")) {
      minprimer[[1]][m] = "a"
    } else {
      minprimer[[1]][m] = "c"
    }
    
  }
  
  max_Gc <-
    specify_decimal(calculate_GC(as.character(c2s(unlist(maxprimer)))), 2)
  min_Gc <-
    specify_decimal(calculate_GC(as.character(c2s(unlist(minprimer)))), 2)
  
  return(c(max_Gc, min_Gc))
  
}