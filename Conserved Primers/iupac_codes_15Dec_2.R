
library(readr)

# date: 15 Dec 2022

#Note: Variation annotation computing is based on 6854423 high quality human sequences of 29769 variation sites.
df <- read_csv("MutationData15Dec.csv", 
                              col_types = cols(...1 = col_skip()))


mutasyon <-data.frame(df[,1:4])
#--------
# in covrimer 6 and further version, do not do this, use SNP counts, since some positions dont have snp but have deletion/insertion
# change snp count to total count
mutasyon[,4] <- df[,16]
#--------
mutasyon[which(mutasyon[,3]==""),4]=0  ### burayı her zaman kontrol et, bu defa mutasyon olmayan yer yoktu. 
mutasyon <- mutasyon[-(which(mutasyon[,3]=="")),]
# find reference (for double check) and altered nucleotides and add them to column 5 and 6
for (i in 1:dim(mutasyon)[1]){
  
  a <- mutasyon[i,3] # "C>G:249,C>B:3,C>A:17,C>M:3,C>K:1,C>H:1,C>V:3,C>T:80,C>S:18,C>Y:2"
  aa <- unlist(strsplit(a, ":")) # "C>G"     "249,C>B" "3,C>A"   "17,C>M"  "3,C>K"   "1,C>H"   "1,C>V"   "3,C>T"   "80,C>S"  "18,C>Y"  "2"      
  aaa <- unlist(strsplit(aa, ",")) # "C>G" "249" "C>B" "3"   "C>A" "17"  "C>M" "3"   "C>K" "1"   "C>H" "1"   "C>V" "3"   "C>T" "80"  "C>S" "18"  "C>Y" "2"  
  ind <- which(aaa %in% c(1:100000)) # find places for numbers  #  2  4  6  8 10 12 14 16 18 20
  
  aaa <- aaa[-ind] # remove them # "C>G" "C>B" "C>A" "C>M" "C>K" "C>H" "C>V" "C>T" "C>S" "C>Y"
  C <- sum(startsWith(aaa,"C>"))
  A <- sum(startsWith(aaa,"A>"))
  TT <- sum(startsWith(aaa,"T>"))
  G <- sum(startsWith(aaa,"G>"))
  if (C>TT && C>A && C>G ){ ref <- "C"}
  if (A>TT && A>C && A>G ){ ref <- "A"}
  if (G>TT && G>A && G>C ){ ref <- "G"}
  if (TT>C && TT>A && TT>G ){ ref <- "T"}
  
  altC<- sum(endsWith(aaa,">C"))
  altA <- sum(endsWith(aaa,">A"))
  altT <- sum(endsWith(aaa,">T"))
  altG <- sum(endsWith(aaa,">G"))
  
  # ambiguous nt are not used. 
  # altY <- sum(endsWith(aaa,">Y"))
  # altM <- sum(endsWith(aaa,">M"))
  # altW <- sum(endsWith(aaa,">W"))
  # altK <- sum(endsWith(aaa,">K"))
  # altS <- sum(endsWith(aaa,">S"))
  # 
  # altR <- sum(endsWith(aaa,">R"))
  # altH <- sum(endsWith(aaa,">H"))
  # altB <- sum(endsWith(aaa,">B"))
  # altV <- sum(endsWith(aaa,">V"))
  # altD <- sum(endsWith(aaa,">D"))
  # altN <- sum(endsWith(aaa,">N"))
  #altY,altM,altW,altK,altS,altR,altH,altB,altV,altD,altN
  #"Y","M","W","K","S","R","H","B","V","D","N"
  
  alts <- data.frame(c(altC,altA,altG,altT))
  rownames(alts) <- c("C","A","G","T")
  indd <- which(alts[,1]>0)
  al <- rownames(alts)[indd]
  alt <-paste(al,sep = ",", collapse = ",")
  
  
  mutasyon[i,5] <- ref
  mutasyon[i,6] <- alt
}

colnames(mutasyon)[5:6]=c("Reference","Altered")

# all possible nucleotides for each position (reference+altered)
for (i in 1:dim(mutasyon)[1]){
  mutasyon[i,7] <- paste0(paste(substring(mutasyon[i,5],1:2,2), collapse=","),mutasyon[i,6] ,sep = ",", collapse = ",") 
} # 7th column contains ref+altered


library(Biostrings)
# finding iupac
sub <- data.frame(mutasyon[,7])
sub[,2]=""
for (i in 1:dim(sub)[1]){
  cod <- c(gsub(',', '', sub[i,1]))
  sub[i,2] <- mergeIUPACLetters(cod)
  
}
#paste0(mutasyon[i,5],mutasyon[i,6], sep = ",", collapse = ",") 
#final[which(final$mutasyon...7.)]
# find number of changes in nt, 1:only 1 altered, 2: two altered, 3: 3 altered, so with reference it makes N
sub[,"number of changes"]=""
for (i in 1:dim(sub)[1]){
  sub[i,3] <- length(strsplit(sub[i,1],",")[[1]])-1
  
}
txt <- "a test of capitalizing"
gsub("(\\w)(\\w*)", "\\U\\1\\L\\2", txt, perl=TRUE)


withiupac <- cbind(mutasyon,sub)
#remove ,
t <- withiupac[which(withiupac$Altered==""),7] 
withiupac[which(withiupac$Altered==""),7] <-   gsub("\\,","",t) 



# write.csv(withiupac, "withIUPAC-25-08.csv", row.names = F)

mper <- 6854423
# to make it similar to previous nexts data
# its column names: Empty, ref, alt, iupac, change, mutation percentage
final <- withiupac[,c(1:2,7,9:10)]
sss <- mutasyon[,4]*100/mper # mutation percentage
final <- cbind(final,sss)
colnames(final) <- c("position" ,"reference.nt" ,"Altered nt.", "iupac", "number of changes" ,"mutation percentage"  )


write.csv(final,"allpositions-iupac-change-percentDec15.csv",row.names = F) # for all positions having mutations

# for only positions greater than 0.1% as threshold
final2 <- final[which(final[,6]>0.1),]  #1770

write.csv(final2,"Dec15_0.1threshold-iupac-change-percent.csv",row.names = F)

#THIS IS FINAL NEXTS DATA WHICH WILL BE USED IN APP!!!!!!




