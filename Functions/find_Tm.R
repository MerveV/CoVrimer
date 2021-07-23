find_Tm <- function(input) {
  primer <- strsplit(as.character(input), "")
  rv_copy2 <- primer
  rv_copy <- primer
  for (m in 1:length(primer[[1]])) {
    # for max primer Tm
    
    if (primer[[1]][m] %in% c("g", "R", "S", "K", "B", "D")) {
      rv_copy2[[1]][m] = "g"
    } else if (primer[[1]][m] %in% c("c", "Y", "M", "V", "H")) {
      rv_copy2[[1]][m] = "c"
    } else if (primer[[1]][m] %in% c("a", "W")) {
      rv_copy2[[1]][m] = "a"
    } else {
      rv_copy2[[1]][m] = "t"
    }
    
    
    
    # for min primer Tm
    
    if (primer[[1]][m] %in% c("g", "S")) {
      rv_copy[[1]][m] = "g"
    } else if (primer[[1]][m] %in% c("t", "Y", "K", "B", "D", "H")) {
      rv_copy[[1]][m] = "t"
    } else if (primer[[1]][m] %in% c("a", "R", "W", "M", "V")) {
      rv_copy[[1]][m] = "a"
    } else {
      rv_copy[[1]][m] = "c"
    }
    
  }
  
  max_Tm <-
    specify_decimal(as.numeric(calculate_tm(as.character(rv_copy2))), 2)
  min_Tm <-
    specify_decimal(as.numeric(calculate_tm(as.character(rv_copy))), 2)
  
  return(c(max_Tm, min_Tm))
  
}