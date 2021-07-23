specify_decimal <- function(x, k) {
  trimws(format(round(x, k), nsmall = k)) # make desired decimal in text part
}