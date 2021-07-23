# CoVrimer
### **CoVrimer** is a Shiny R based application that allows searching for published/novel primer pairs and visualizing their alignment to the reference genome. 
### CoVrimer shows the amplicon and  any mutation found in the amplicon as well as in the primer binding regions. CoVrimer also provides a list of in-house designed conserved and degenerate primer pairs across the viral genome and presents information on occurrence of mutations, alignment to the genome, and conservation of the binding site among other bat genomes.

Variation data were taken from [CNCN-NGDC](https://bigd.big.ac.cn/ncov/) on May 12, 2021. Variation annotation computing was performed based on 954,725 high quality SARS-CoV-2 sequences.

In our analysis, positions exhibiting less than 0.1% mutation frequency were considered conserved. 


-------------------------------------------------------------------------------------


*Degenerate primers were designed to include 82 variants from 10 different SARS-CoV-2 strains (B.1.1.7, P.1, P.2, B.1.351, B.1.617.1, B.1.526, B.1.526.1, B.1.526.2, B.1.427 and B.1.429)*





-------------------------------------------------------------------------------------


 <span style="color:red"> **Red dots (above the threshold (0.1%))**</span> and <span style="color:orange"> **orange dots (below the threshold (0.1%))**</span> represent mutations in nucleotides.

 <span style="color:blue"> **Blue**</span> 
 and <span style="color:green"> **green**</span> letters refer to the nucleotides in the reference genome (NC_045512) and indicate primer binding sites and the region between primers, respectively.

 Next to the points are the degenerate IUPAC codes for each alteration.
                              
PhyloP values (44 Bat virus strains Basewise Conservation) were obtained from UCSC.


-------------------------------------------------------------------------------------

**IUPAC Codes**

**R:**	A or G
**Y:**	C or T
**S:** G or C
**W:**	A or T
**K:**	G or T
**M:**	A or C
**B:**	C or G or T
**D:**	A or G or T
**H:**	A or C or T
**V:**	A or C or G
**N:**	any base


-------------------------------------------------------------------------------------
