library(readxl)
library(dplyr)
library(stringr)
library(tidyr) 
library(seqinr)
library(readxl)
getwd()
setwd("~/Desktop/Covrimer versions and others in desktop/CoVrimer folders/CoVrimer 7")
newvariation <- read_excel("variation-annotation.xlsx", 
                              col_types = c("text", "text", "numeric", 
                                            "text", "text", "text", "skip", "skip", 
                                            "skip", "skip", "skip", "skip", "skip"))

# date:15 Dec 2022

# Note: Variation annotation computing is based on 6854423 high quality human sequences of 29769 variation sites.
newvariation <- data.frame(newvariation)
library(readr)
tot_num <- 6854423
# newvariation$number <- NA
# for (r in 1:nrow(newvariation)) {
#   newvariation$number [r] <- parse_number(newvariation[r,6])
# }
newvariation2 <- newvariation%>%mutate(number=parse_number(Base.changes.Virus.Number))

df<- data.frame(unique(newvariation2[,1]))
colnames(df)[1] <- "position"
df[,1] <- as.numeric(df[,1])
# The reference nucleotides were added as a column, using the corresponding positions of the reference genome 
ref <- read.fasta("reference.fasta")
r <- ref$NC_045512
df$`reference nt` <- NA
for (n in 1:nrow(df)) {
  df$`reference nt`[n] <- toupper(r[as.numeric(df$position[n])])
}

df$`SNP` <- NA
df$`SNP count` <- NA
df$`SNP freq` <- NA
df$`insertion` <- NA
df$`insertion count` <- NA
df$`insertion freq` <- NA
df$`deletion` <- NA
df$`deletion count` <- NA
df$`deletion freq` <- NA
df$`indel` <- NA
df$`indel count` <- NA
df$`indel freq` <- NA

for (n in 1:nrow(df)) {
  
  pos <- df[n,1] 
  
  d <- newvariation2[which(newvariation2[,1]==pos),]
  
  snp <- paste0(d[which(d[,5]=="SNP"),6], sep = ", ",collapse="")
  df$`SNP`[n] <-snp
  df$`SNP count`[n] <- sum(d[which(d[,5]=="SNP"),7])

  inser <- paste0(d[which(d[,5]=="Insertion"),6], sep = ", ",collapse="")
  df$`insertion`[n] <-inser
  df$`insertion count`[n] <- sum(d[which(d[,5]=="Insertion"),7])
  
  dele <- paste0(d[which(d[,5]=="Deletion"),6], sep = ", ",collapse="")
  df$`deletion`[n] <-dele
  df$`deletion count`[n] <- sum(d[which(d[,5]=="Deletion"),7])
  
  ind <- paste0(d[which(d[,5]=="Indel"),6], sep = ", ",collapse="")
  df$`indel`[n] <-ind
  df$`indel count`[n] <- sum(d[which(d[,5]=="Indel"),7])

  summ <- sum(df$`SNP count`[n],df$`insertion count`[n],df$`deletion count`[n],df$`indel count`[n] )
  df$`SNP freq`[n]  <- df$`SNP count`[n]/summ
  df$`insertion freq`[n]  <- df$`insertion count`[n]/summ
  df$`deletion freq`[n]  <- df$`deletion count`[n]/summ
  df$`indel freq`[n]  <- df$`indel count`[n]/summ
}

library(stringr)

df2 <- df
#df2[,3] <- str_sub(df2[,3],1,nchar(df2[,3])-1)
for (i in 1:nrow(df2)) {
  df2[i,3] <- str_sub(df2[i,3],1,nchar(df2[i,3])-2)
  df2[i,6] <- str_sub(df2[i,6],1,nchar(df2[i,6])-2)
  df2[i,9] <- str_sub(df2[i,9],1,nchar(df2[i,9])-2)
  df2[i,12] <- str_sub(df2[i,12],1,nchar(df2[i,12])-2)
}
#%>% mutate(`SNP`=replace(`SNP`, cyl==4, NA))

df <- df2
df$`total freq` <- NA
df$`isolates` <- NA
for (n in 1:nrow(df)) {
  df$`isolates`[n] <-sum(  df$`SNP count`[n],  df$`insertion count`[n],  df$`deletion count`[n],  df$`indel count`[n])
}



df$`total freq`<-df$`isolates`/tot_num


write.csv(df, "MutationData15Dec.csv")
