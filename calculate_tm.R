calculate_tm <- function(input) {
  library(TmCalculator)
  #inp = strsplit(input, "")
  pri=input
  Tm=specify_decimal(Tm_NN(pri,ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = "DNA_NN3",
                           tmm_table = "DNA_NN3", imm_table = "DNA_NN3",de_table = "DNA_NN3", dnac1 = 25,
                           dnac2 = 25, selfcomp = FALSE, Na = 0, K = 50, Tris = 0, Mg = 1.5, dNTPs = 0.6, saltcorr = 5),2)
  return(Tm)
}