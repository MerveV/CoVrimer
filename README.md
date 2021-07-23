# CoVrimer
### **CoVrimer** is a Shiny R based application that allows searching for published/novel primer pairs and visualizing their alignment to the reference genome. 
### CoVrimer shows the amplicon and  any mutation found in the amplicon as well as in the primer binding regions. CoVrimer also provides a list of in-house designed conserved and degenerate primer pairs across the viral genome and presents information on occurrence of mutations, alignment to the genome, and conservation of the binding site among other bat genomes.

Variation data were taken from [CNCN-NGDC](https://bigd.big.ac.cn/ncov/) on May 12, 2021. Variation annotation computing was performed based on 954,725 high quality SARS-CoV-2 sequences.

In our analysis, positions exhibiting less than 0.1% mutation frequency were considered conserved. 

*Degenerate primers were designed to include 82 variants from 10 different SARS-CoV-2 strains (B.1.1.7, P.1, P.2, B.1.351, B.1.617.1, B.1.526, B.1.526.1, B.1.526.2, B.1.427 and B.1.429)*





-------------------------------------------------------------------------------------




**Country:** Country of the study

**Panel:** Panel names or given name for that primer pair (CDC Panel (Centers for Disease Control and Prevention), HKU panel, NMDC Panel (National Microbiology Data Center ), WHO Panel (World Health Organization)

**Target Gene:** The gene where primer pair binds to 

**Label:** Forward, Reverse and Probe labels

**Primer Sets:** Sequences of the oligonucleotides in the Primer set

**Degenerate Version:**  Sequences of the primers in the primer set with degenerate nucleotides incorporated at positions exhibiting mutation frequency over threshold

**Positions:** Start and end position of primers 

**Amplicon Size:** Length of the product 

**Mutation Number in Fw/Rv:** Number of mutated positions in forward/reverse primer binding region

**Mutation Frequency in Fw/Rv:** Frequency of mutation events at each position of forward/reverse primer binding region

**Reference Genome:** Genome used in the study for primer design

**Used Programs:** The programs used for sequence alignment and primer design 

**Type of Test:** Type of diagnostic test for which primers are designed and used

**PMID:** Article Identifiers for PUBMED  

**Date:** Publication date of the article

**Link:** Easy access to the article
