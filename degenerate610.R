library(tidyverse)
library(dplyr)
library(readr)
library(stringr)
library(peRspective)
library(Xmisc)

# All degenerate pairs from deg.R
degen <- read_csv("degforlineages_(04.05).csv")

# Vector of mutation frequencies
mutdata <-  read_csv("MutationDataMay12.csv")
mutdata[,1] <-  NULL
mutfreq <- c(rep(0,29903))
for (i in 1:nrow(mutdata)) {
  mutfreq[mutdata$position[i]] = as.numeric(mutdata$`total freq`[i])
}
rm(mutdata)

# SNP positions from lineages
lineagemut <- read_csv("lineagemutations.csv")
colnames(lineagemut) <- c("reg","aamut","nucmut","mutpos")

# 5 conserved nucleotides at the 3' end
degcons3 <-  degen  %>%
  rowwise() %>% 
  filter(all(mutfreq[(fwend-4) :fwend] <= 0.001)) %>%
  filter(all(mutfreq[revend :(revend+4)] <= 0.001))

# Tm, GC, ampliconlength, primerlength
degpairs1 <- degcons3 %>%
  filter(amplicon_length <= 180)%>%
  filter(fwGC >= 40)%>%
  filter(revGC >= 40)%>%
  filter(between(fwTm,55,65))%>%
  filter(between(revTm,55,65))%>%
  filter(fwlength <= 25) %>%
  filter(revlength <= 25)

# SARS Specificity
ddd <- read.csv("SARS-specifityamong44virus.csv",header = TRUE) # n*100/top

# 10th column
for ( i in 1: dim(ddd)[1]){
  ddd[i,10] <- ddd[i,7]/44
} # Not percentage, just difference/total

ss <- degpairs1
pub.for=c()
pub.rev=c()

for (i in 1:nrow(ss)){
  n = (ss$fwstart[i]):(ss$fwend[i]) # forward primer position
  lF = (ss$fwend[i]) - (ss$fwstart[i]) + 1 # forward primer length
  m = (ss$revend[i]):(ss$revstart[i]) # reverse primer position
  lR = (ss$revstart[i]) - (ss$revend[i]) + 1# reverse primer length

  conspercF = (sum(ddd[n[1:(lF-5)],10])/(lF-5))
  conspercR = (sum(ddd[m[6:lR],10])/(lR-5))
  
  s1 = as.double(specify_decimal(conspercF,4))
  s2 = as.double(specify_decimal(conspercR,4))
  
  last5forw = ddd[n[(lF-4): lF],10]  # last 5 nt 
  last5rev = ddd[m[1:5],10]
  
  last5f = sum(last5forw * (1 + last5forw ))/5
  last5r = sum(last5rev * (1 + last5rev ))/5
  
  s1 = s1 + last5f
  s2 = s2 + last5r
  
  pub.for <- append(pub.for,s1)
  pub.rev <- append(pub.rev,s2)
}

ss <- cbind(ss,pub.for,pub.rev)
sss <- ss[which(ss[,26]!=0 & ss[,27]!=0),] # remove if SARS-specifitcy is 0


# Lineage mutation positions
formutinc <- sss %>%
  rowwise() %>%
  mutate(fwincluding = str_flatten(lineagemut$mutpos[which(lineagemut$mutpos %in% (fwstart:fwend))],collapse = ",")) %>%
  mutate(revincluding = str_flatten(lineagemut$mutpos[which(lineagemut$mutpos %in% (revend:revstart))],collapse = ",")) %>%
  mutate(fwcons = if_else(all(mutfreq[fwstart:fwend] <= 0.001),"+","-"))%>%
  mutate(revcons = if_else(all(mutfreq[revend:revstart] <= 0.001),"+","-"))

# Mutated positions
newfreq <- mutfreq
newfreq[which(newfreq <= 0.001)] <- 0
newfreq[lineagemut$mutpos] <- 1

a <- c(1:29903)
newfreq <- bind_cols (a,newfreq)
colnames(newfreq) <- c("pos","freq")
rrr <- newfreq[newfreq$freq != 0,]

formutinc_ <- formutinc %>%
  rowwise() %>%
  mutate(fwincluding = str_flatten(rrr$pos[which(rrr$pos %in% (fwstart:fwend))],collapse = ",")) %>%
  mutate(revincluding = str_flatten(rrr$pos[which(rrr$pos %in% (revend:revstart))],collapse = ","))


ff1 <- formutinc_ %>%
  filter(fwincluding != "") %>%
  filter(revcons == "+") %>%
  group_by(fwincluding) %>%
  slice(which.max(pub.rev)) %>%
  ungroup(fwincluding)

ff2 <- formutinc_ %>%
  filter(revincluding != "") %>%
  filter(fwcons == "+") %>%
  group_by(revincluding) %>%
  slice(which.max(pub.for)) %>%
  ungroup(revincluding)

fpairs <- bind_rows(ff1,ff2)
fpairs <- unique(fpairs)

ff3 <- formutinc_ %>%
  filter(fwincluding != "") %>%
  filter(revcons == "-") %>%
  group_by(fwincluding) %>%
  slice(which.max(pub.rev)) %>%
  ungroup(fwincluding)

ff4 <- formutinc_ %>%
  filter(revincluding != "") %>%
  filter(fwcons == "-") %>%
  group_by(revincluding) %>%
  slice(which.max(pub.for)) %>%
  ungroup(revincluding)

fpairs1 <- bind_rows(ff3,ff4)
fpairs1 <- unique(fpairs1)


ffin <- unique(bind_rows(fpairs,fpairs1))
unique(ffin$revincluding)
unique(ffin$fwincluding)

write_csv(ffin,"degenerate610.csv")

