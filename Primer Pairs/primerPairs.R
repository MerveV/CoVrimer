#Primer Pairs - Unfiltered

library(dplyr)
library(rDNAse)
library(readr)

# Primer pairs were filtered so that their total mutation frequencies are below 0.001
# Last 5 nucleotides of the primers were scanned and the counts of nucleotides with 0 total frequency were added as a column
# The resulting primer pairs were saved as a single excel file 

nsp1 <- read.csv("nsp1primers.csv")
nsp2 <- read.csv("nsp2primers.csv")
nsp3_1 <- read.csv("nsp3primers1.csv")
nsp3_2 <- read.csv("nsp3primers2.csv")
nsp3_3 <- read.csv("nsp3primers3.csv")
nsp3_4 <- read.csv("nsp3primers4.csv")
nsp3_5 <- read.csv("nsp3primers5.csv")
nsp3_6 <- read.csv("nsp3primers6.csv")
nsp3_7 <- read.csv("nsp3primers7.csv")
nsp3_8 <- read.csv("nsp3primers8.csv")
nsp3_9 <- read.csv("nsp3primers9.csv")
nsp3_10 <- read.csv("nsp3primers10.csv")
nsp3_11 <- read.csv("nsp3primers11.csv")
nsp3_12 <- read.csv("nsp3primers12.csv")
nsp3_13 <- read.csv("nsp3primers13.csv")
nsp3_14 <- read.csv("nsp3primers14.csv")
nsp3_15 <- read.csv("nsp3primers15.csv")
nsp4 <- read.csv("nsp4primers.csv")
nsp5 <- read.csv("nsp5primers.csv")
nsp6 <- read.csv("nsp6primers.csv")
nsp7 <- read.csv("nsp7primers.csv")
nsp8 <- read.csv("nsp8primers.csv")
nsp9 <- read.csv("nsp9primers.csv")
nsp10 <- read.csv("nsp10primers.csv")
nsp11.12 <- read.csv("nsp11.12primers.csv")
nsp13 <- read.csv("nsp13primers.csv")
nsp14 <- read.csv("nsp14primers.csv")
nsp15 <- read.csv("nsp15primers.csv")
nsp16 <- read.csv("nsp16primers.csv")
orf3a <- read.csv("orf3aprimers.csv")
s_1 <- read.csv("Sprimers1.csv")
s_2 <- read.csv("Sprimers2.csv")
s_3 <- read.csv("Sprimers3.csv")
s_4 <- read.csv("Sprimers4.csv")
s_5 <- read.csv("Sprimers5.csv")
e <- read.csv("Eprimers.csv")
m <- read.csv("Mprimers.csv")
n <- read.csv("Nprimers.csv")
orf6 <- read.csv("ORF6primers.csv")
orf7a <- read.csv("ORF7aprimers.csv")
orf7b <- read.csv("ORF7bprimers.csv")
orf8 <- read.csv("ORF8primers.csv")
orf10 <- read.csv("ORF10primers.csv")


# Mutation data was loaded
mutdata <- read.csv("MutationDataMay12.csv")
mutdata$X <-  NULL
mutfreq <- c(rep(0,29903))
for (i in 1:nrow(mutdata)) {
  mutfreq[mutdata$position[i]] = as.numeric(mutdata$total.freq[i])
}


  merg_e <- e  %>%
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "e") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  

  merg_m <- m  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "m") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_n <- n  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "n") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp1 <- nsp1  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp1") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp10 <- nsp10  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp10") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp13 <- nsp13  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp13") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp14 <- nsp14  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp14") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp15 <- nsp15  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%  
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp15") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp16 <- nsp16  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp16") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>% 
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp2 <- nsp2  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp2") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.1 <- nsp3_1  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.2 <- nsp3_2  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.3 <- nsp3_3  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.4 <- nsp3_4  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.5 <- nsp3_5  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.6 <- nsp3_6  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.7 <- nsp3_7  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%  
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.8 <- nsp3_8  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.9 <- nsp3_9  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.10 <- nsp3_10  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.11 <- nsp3_11  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%  
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.12 <- nsp3_12  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  

  merg_nsp3.13 <- nsp3_13  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>% 
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.14 <- nsp3_14  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>% 
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


  merg_nsp3.15 <- nsp3_15  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp3") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_nsp4 <- nsp4  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp4") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_nsp5 <- nsp5  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp5") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_nsp6 <- nsp6  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp6") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_nsp7 <- nsp7  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp7") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
 
  merg_nsp8 <- nsp8  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp8") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_nsp9 <- nsp9  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%  
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "nsp9") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_orf3a <- orf3a  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "orf3a") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_orf6 <- orf6  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%  
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "orf6") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>% 
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
   
  merg_orf7a <- orf7a  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "orf7a") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))

  
  merg_orf7b <- orf7b  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "orf7b") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_orf8 <- orf8  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "orf8") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_s.1 <- s_1  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "s") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_s.2 <- s_2  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>%   
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "s") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_s.3 <- s_3  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%   
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "s") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_s.4 <- s_4  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "s") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  
  
  merg_s.5 <- s_5  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>% 
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "s") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))
  

  merg_orf10 <- orf10  %>% 
    mutate(revprimer = revchars(revprimer)) %>%  rowwise() %>%  
    filter(all(mutfreq[fwstart :fwend] <= 0.001)) %>% 
    filter(all(mutfreq[revstart :revend] <= 0.001)) %>%
    mutate(region = "orf10") %>%
    mutate(forwardmut = sum(mutfreq[fwstart:fwend]))%>%
    mutate(reversemut = sum(mutfreq[revend:revstart]))%>%
    mutate(total=forwardmut+reversemut) %>%
    mutate(fwfreq = paste(mutfreq[fwstart :fwend],collapse = ",")) %>%  
    mutate(revfreq = paste(mutfreq[revstart :revend],collapse = ",")) %>%
    mutate(fwlast5nuc0count = length(which(mutfreq[(fwend-4):fwend]==0))) %>%
    mutate(revlast5nuc0count = length(which(mutfreq[revend :(revend+4)]==0))) %>%
    mutate(fwmutavg = forwardmut/nchar(fwprimer))%>%
    mutate(revmutavg = reversemut/nchar(revprimer))%>%
    mutate(totalavg = (fwmutavg+revmutavg)/(nchar(fwprimer)+nchar(revprimer)))


under1000 <- bind_rows(merg_e,merg_m,merg_n,merg_nsp1,merg_nsp10,merg_nsp13,merg_nsp14,merg_nsp15,merg_nsp16,merg_nsp2,merg_nsp3.1,
                       merg_nsp3.10,merg_nsp3.11,merg_nsp3.12,merg_nsp3.13,merg_nsp3.14,merg_nsp3.15,merg_nsp3.2,merg_nsp3.3,
                       merg_nsp3.4,merg_nsp3.5,merg_nsp3.6,merg_nsp3.7,merg_nsp3.8,merg_nsp3.9,merg_nsp4,merg_nsp5,merg_nsp6,
                       merg_nsp7,merg_nsp8,merg_nsp9,merg_orf3a,merg_orf6,merg_orf7a,merg_orf8,merg_s.1,merg_s.2,merg_s.3,merg_s.4,
                       merg_s.5)


head(under1000)
write_csv(under1000,"under1000May12.csv")
