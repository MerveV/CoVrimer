# downloaded data, total sequence is 954725
library(readxl)
library(seqinr)
source("calculate_GC.R")
source("calculate_tm.R")
source("find_Tm.R") # max/min
source("find_Gc.R") # max/min
source("forw_degenerate.R") # replaced source("find_degenerate.R")
source("reverse_degenerate.R") # replaced source("find_degenerate_reverse.R")
source("specify_decimal.R")
v <- read_excel("~/Downloads/variation-annotation (1).xlsx")

#  threshold 0.001 = 0.1%
vv <- data.frame(cbind(v[,c(1,3)], v[,3]/954725))
vv[,3] <- round(vv[,3],5)
vv[,1] <- as.numeric(vv[,1])
v0.1 <- vv[which(vv[,1]>265),] # remove 5'UTR
v0.1 <- v0.1[which(v0.1[,1]<29675),] # remove 3'UTR
# 266..21555 gene="ORF1ab"
# 29675..29903 3'UTR

# 0.001 den büyük olanları al
v0.1 <- v0.1[which(v0.1[,3]>0.00100),] 

diff=NULL
for (i in 1:(dim(v0.1)[1]-1)){
  diff= rbind.data.frame(diff,c(v0.1[i,1],v0.1[i+1,1]))
}

diff[,3]=diff[,2]-diff[,1]-1
# aradaki farkı veren mutasyonsuz pozisyonları yazdır
diff[,1]=diff[,1]+1
diff[,2]=diff[,2]-1
diff17 <- diff[which(diff[,3]>17),] # 17bp den uzun olanları seç
library(openPrimeR)
# OpenPrimeR supplies predefined XML files specifying default settings for different applications:
list.files(system.file("extdata", "settings", package = "openPrimeR"), pattern = "*\\.xml")
# select the high-stringency primer design conditions for Taq polymerase and load the DesignSettings object with read_settings()):
settings.xml <- system.file("extdata", "settings", "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
settings <- read_settings(settings.xml)

design.settings <- settings # change name 
#remove the requirement for a GC clamp in the following way:
constraints(design.settings) <- constraints(design.settings)[!grepl("gc_clamp", names(constraints(design.settings)))]
#design only primers of a specific length
constraints(design.settings)[["primer_length"]] <- c("min" = 18, "max" = 26)
# prevent any mismatch binding
conOptions(design.settings)[["allowed_mismatches"]] <- 0

#### View the active constraint limits
constraintLimits(design.settings)
####### change settings 
constraintLimits(design.settings)$gc_clamp <- c("min" = 0, "max" = 4)
PCR(design.settings)$annealing_temperature <- 50 # Evaluate primers with a fixed annealing temperature

constraintLimits(design.settings)$melting_temp_range<- c("min" = 55, "max" = 75)
constraintLimits(design.settings)$melting_temp_diff<- c("min" = 0, "max" = 5)
constraintLimits(design.settings)$no_runs<- c("min" = 0, "max" = 3)


library(DECIPHER)
dbConn <- dbConnect(SQLite(), ":memory:") # create database
Seqs2DB("/Users/mervevural/Desktop/CoVrimer Revision/reference_genome.fasta", type = "FASTA", dbFile = dbConn, identifier = "reference")
dna <- SearchDB(dbConn)



cons5=NULL # create amplicon 
for (i in 1:(dim(diff17)[1]-1)){
  for (j in 1:5){
    cons5=rbind.data.frame(cons5,c(diff17[i,],diff17[i+j,]))
    
  }}
# max amplicon lenght in column 7 
cons5[,5]=as.numeric(cons5[,5])
cons5[,1]=as.numeric(cons5[,1])
cons5[,7]=cons5[,5]-cons5[,1]-1

colnames(cons5) <- c("FW start","FW end", "FW cons","RV star", "RV end", "RV cons","Amplicon")
cons6<- cons5[ which(cons5[,7]<181),]
cons6<- cons6[ which(cons6[,7]>79),]


library(dplyr)
final <- NULL
t <- cons6
for (r in 1:dim(t)[1]) {
  newsequence=dna$`1`[t[r,1]:t[r,5]]
  write.fasta(newsequence, names= "cons", file="seq.FASTA") # first save seq
  seq <- read_templates("seq.FASTA") # load seq
  if (t[r,3]<26) { asg1=t[r,3]}else { asg1=26}
  
  if (t[r,6]<26){ asg2=t[r,6]} else{asg2=26}
  
  seq<- assign_binding_regions(seq, fw = c(1,asg1), rev = c(1,asg2))
  int_forward=get_initial_primers("InitialPrimers", seq, c(18,asg1),"fw") 
  
  # forward start and end 
  forw1=as.data.frame(int_forward[["Forward"]])
  forw1[,2]=int_forward[["primer_length_fw"]]
  allowf=seq[["Allowed_End_fw_initial"]]-seq[["Allowed_Start_fw_initial"]]+1
  if (t[r,3]<19){
    forw1[,3]=t[r,1]
    forw1[,4]=forw1[,3]+forw1[,2]-1
  }else {
    mf=allowf  # position of allowed 
    mmf=NULL
    for (i in 1:(mf-18)){
      for (k in 18:t[r,3]) {
        a=i+k-1
        if (a<=mf){
          mmf=c(mmf,a)
        }
      }
    }
    forw1[,3]=mmf-forw1[,2]+ t[r,1]# start of forward primer 
    forw1[,4]=forw1[,3]+forw1[,2]-1# end of forward primer 
    
  }
  colnames(forw1) <- c("fw","fwl","fwstart","fwend")
  
  forw1 <- forw1%>% rowwise() %>% 
    mutate(fwtm=specify_decimal(as.numeric(calculate_tm(fw)),2),
           fwgc=specify_decimal(as.numeric(calculate_GC(fw)),2)) 
  
  
  # REVERSE
  
  int_reverse=get_initial_primers("InitialPrimers", seq, c(18,asg2),"rev")
  rev1=as.data.frame(int_reverse[["Reverse"]])
  rev1[,2]=int_reverse[["primer_length_rev"]]
  allow=seq[["Allowed_End_rev_initial"]]-seq[["Allowed_Start_rev_initial"]]+1
  if (t[r,6]<19){
    rev1[,3]=t[r,5]-t[r,6]+1
    rev1[,4]=t[r,5]
  }else {
    m=length(newsequence)-allow +1
    mm=NULL
    for (i in length(newsequence):(m+18)){
      for (k in 18:t[r,6]) {
        a=i-k+1
        if (a>=m){
          mm=c(mm,a)
        }
      }
    }
    rev1[,3]=mm+t[r,1]-1 # to find start
    rev1[,4]=rev1[,3]+rev1[,2]-1 # to find end 
  }
  
  
  colnames(rev1) <- c("rv","rvl","rvstart","rvend")
  as_tibble(rev1)
  
  rev1 <- rev1%>% rowwise() %>% 
    mutate(rvtm=specify_decimal(as.numeric(calculate_tm(rv)),2),
           rvgc=specify_decimal(as.numeric(calculate_GC(rv)),2)) 
  
  cat(r,"r")
  print(dim(final))
  
  f <- forw1%>% 
    filter(fwtm>40)%>% 
    filter(fwgc>30)
  rev <- rev1%>% 
    filter(rvtm>40)%>% 
    filter(rvgc>30)
  
  final <- rbind(final,merge(f,rev))
  
  
}






final[,"amplicon"] <- final[,10]-final[,3]+1
final[,14]=specify_decimal((abs(as.numeric(final[,5])- as.numeric(final[,11]))),2)
colnames(final) <- c("fw_primer","fw_length","fw_start","fw_end", "fw_tm","fw_GC",
                     "rv_primer","rv_length","rv_start","rv_end", "rv_tm","rv_GC", "amplicon","TMdiff" )

dim(final) # 627514
# Tm, GC, ampliconlength, primerlength
final1 <- final %>%
  filter(fw_GC >= 40)%>%
  filter(rv_GC >= 40)%>%
  filter(TMdiff < 5)%>%
  filter(between(fw_tm,55,65))%>%
  filter(between(rv_tm,55,65))%>%
  filter(between(amplicon,60,180))

dim(final1) # 65232 
# final2 <- final1%>%
#   filter(fwlength <= 25) %>%
#   filter(revlength <= 25)

# --------------FUNCTION run / no more 4bp run------------
library(stringr)
run<- function(input,nt){
  a<- str_locate_all(input,nt)
  a=as.data.frame(a)
  a=a[,1]
  count=1
  if (length(a)>1){
    for (i in 1:(length(a)-1) ){
      if ( a[i+1]==(a[i]+1) ) { 
        count=count+1
        if (count>5){
          break
        }}
      else {count=1}
    }
  }
  else { count=1}
  
  
  
  return(count)
}

end <- final1
final1run <- NULL
for (i in 1:dim(end)[1]){
  if (run(end[i,1], "a")<5){
    if (run(end[i,1], "t")<5){
      if (run(end[i,1], "g")<5){
        if (run(end[i,1], "c")<5){
          if (run(end[i,7], "a")<5){
            if (run(end[i,7], "t")<5){
              if (run(end[i,7], "g")<5){
                if (run(end[i,7], "c")<5){
                  final1run=rbind(final1run,end[i,])
                }
              }
            }
          }
        }
      }
    }
  }
}
dim(final1run) #62729  

# --------------FUNCTION GC clamp in last 5 bases, not more than 3------------
# create another column that contains last 5 nt at 3' end
input <- final1run[1,1]
start <- final1run[1,3]
end <- final1run[1,2]
library(bioseq)
GCclamp <- function(input,length){
  pos <- seq_extract_position(dna(input), length-4, length)
  return(seq_count_pattern(pos, "C") + seq_count_pattern(pos, "G"))
}
colnames(final1run)

final1runclamp <- final1run%>% rowwise()%>% 
  mutate(forwGCclamp=GCclamp(fw_primer,fw_length),
         reveGCclamp=GCclamp(rv_primer,rv_length))

final1runclamp <- final1runclamp%>%
  filter(forwGCclamp<4)%>%
  filter(reveGCclamp<4)
dim(final1runclamp) # 51217 

# avoid ending T
final1runclampt <- final1runclamp%>%
  filter(!endsWith(fw_primer,"t"))%>%
  filter(!endsWith(rv_primer,"t"))

dim(final1runclampt) #24727    16
  
library(readr)

##### Read files
ddd=read.csv("SARS-specifityamong44virus.csv",header = TRUE) # n*100/top
for ( i in 1: dim(ddd)[1]){
  ddd[i,10]=ddd[i,7]/44
} # not percentage, just difference/total
# 0 means, that nt is same for all others, no specific

# function for finding SARS-Specificity (weighted)-----------
sarsspec <- function(s){
  pub.for=c()
  pub.rev=c()
  for (i in 1:dim(s)[1]){
    n=s[i,3]:s[i,4] # forward primer position
    lF=s[i,2] # forward primer length 
    m=s[i,9]:s[i,10] # reverse primer position
    lR= s[i,8]# reverse primer length 
    
    
    conspercF=(sum(ddd[n[1:(lF-5)],10])/(lF-5))
    conspercR=(sum(ddd[m[6:lR],10])/(lR-5))
    
    s1=as.double(specify_decimal(conspercF,4))
    s2=as.double(specify_decimal(conspercR,4))
    
    last5forw=ddd[n[(lF-4): lF],10]  # last 5 nt 
    last5rev= ddd[m[1:5],10]
    
    last5f=sum ( last5forw * (1+  last5forw ))/5
    last5r=sum ( last5rev * (1+  last5rev ))/5
    
    s1=s1+last5f
    s2=s2+last5r
    
    
    pub.for <- append(pub.for,s1)
    pub.rev <- append(pub.rev,s2)
  }
  
  ss <- cbind(s,pub.for,pub.rev)
  return(ss)
}
final1runclamptd <- data.frame(final1runclampt)
ss <- sarsspec(final1runclamptd) # add SARS-Specif


sss <- ss%>%
  filter(pub.for!=0 &pub.rev !=0)
  
  
dim(sss) #  24067 

sss <- sss%>%
  filter(`Fw GC`<60)%>%
  filter(`Rv GC`<60)

# removing mut count <652 as 0, conserved
data <- data.frame(1:29903,0)
data[vv[,1],2] <- vv[,3]
data[,3] <- 0
data[which(data[,2]>0.00100),3] <- 1

data[28776:28794,]
# 1760 positions have mutations >652
cc <- NULL;xx=sss
for (i in 1:dim(xx)[1]) {
  n <- xx[i,3]:xx[i,10]
  cc <- append(cc,sum(data[n,3]))
}
sss <- cbind(sss[,c(1:12,14,13)],cc,sss[,17:18])
colnames(sss)=c("Forward Primer","Fw Length", "Fw Start", "Fw End", "Fw Tm", "Fw GC",
                "Reverse Primer","Rv Length", "Rv Start", "Rv End", "Rv Tm", "Rv GC",
                "Tm Difference", "Amplicon Length", "Amplicon Mutations","Fw SARS-CoV-2 Specificity","Rv SARS-CoV-2 Specificity")


# final Filtering
library(dplyr)
length(unique(sss[,1]))

forreverse <- sss %>%group_by(`Reverse Primer`)%>%
  slice(which.max(`Fw SARS-CoV-2 Specificity`)) %>%
  ungroup()
forreverse <- forreverse %>%group_by(`Rv Start`)%>%
  slice(which.max(`Fw SARS-CoV-2 Specificity`)) %>%
  ungroup()
forreverse <- forreverse %>%group_by(`Rv End`)%>%
  slice(which.max(`Fw SARS-CoV-2 Specificity`)) %>%
  ungroup()
forfor <- sss %>%group_by(`Forward Primer`)%>%
  slice(which.max(`Rv SARS-CoV-2 Specificity`)) %>%
  ungroup()


forfor <- forfor %>%group_by(`Fw Start`)%>%
  slice(which.max(`Rv SARS-CoV-2 Specificity`)) %>%
  ungroup()

forfor <- forfor %>%group_by(`Fw End`)%>%
  slice(which.max(`Rv SARS-CoV-2 Specificity`)) %>%
  ungroup()

finalpairs <- bind_rows(forreverse,forfor)
finalpairs <- unique(finalpairs)



write.csv(finalpairs,"Conservedprimers_12May.csv",row.names = F)


