library(BiocManager)
options(repos = BiocManager::repositories())
# rsconnect::showLogs(appName="CoVrimer",streaming=TRUE) #to see logs






library(data.table)
# setwd("~/Desktop/CoVrimer 4")
library(shiny)
library(shinycssloaders)
library(dashboardthemes)
library(shinydashboard)
library(DT)
library(msaR)
#library(msa)
library(tidyverse)
library(openPrimeR)
library(readr)
library(ggplot2)
library(DECIPHER)
library(seqinr)
library(readxl)
library("htmltools")
library("bsplus")
library(readr)
library(ape)
library(ggpubr)
library(TmCalculator)
library(shinyBS)

library("dplyr")

source("calculate_GC.R")
source("calculate_tm.R")
source("find_Tm.R") # max/min
source("find_Gc.R") # max/min
source("forw_degenerate.R") # replaced source("find_degenerate.R")
source("reverse_degenerate.R") # replaced source("find_degenerate_reverse.R")
source("specify_decimal.R")

server <- function(input, output, session) {
  
  primers2 <- read_excel("publishedprimers_15Dec.xlsx") # published primers
  
  output$publishedprimerdownload <- downloadHandler(
    filename = function() {
      paste("Publishedprimerspairs", ".csv", sep = "")
    },
    content = function(con) {
      write.csv(primers2, row.names = F, con)
    }
  )
  
  output$details <- renderUI({
    if(input$go){
      includeMarkdown("mdfiles/introdetails.md")
    }
    
  })
  
  output$publishedprimerTable <- DT::renderDataTable ({
    primers2 <- read_excel("publishedprimers_15Dec.xlsx") # published primers
    
    DT::datatable(
      primers2,
      rownames = FALSE,
      escape = FALSE,
      extensions = c( 'Buttons', 'ColReorder'),
      selection = "single",
      
      options = list(
        colReorder = TRUE,
        lengthMenu = c(10, 30, 50),
        pageLength = 3,
        scrollX = TRUE,
        #autoWidth = TRUE,
        dom = 'Blfrtip',
        buttons = I('colvis'), initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#B2B6B9', 'color': '000000'});",
          "}")
      )
    ) %>%
      formatStyle(columns = c(1:18), fontSize = '90%')
    
  })  # primer table   tabPanel= "information of primers"
 
  
   msa1plot <- eventReactive(input$msa1, {
    if (input$publishedprimerTable_rows_selected) {
      input$publishedprimerTable_rows_selected[1]
      
      primers <- read_excel("publishedprimers_15Dec.xlsx")
      s = input$publishedprimerTable_rows_selected[1] # selected row
      primers = as.data.frame(primers)
      
      a = strsplit(as.character(primers[s, 5]), split = "\r\n") #forward primer
      
      reverseprimer = as.character(reverseComplement(DNAString(toString(a[[1]][2])))) #reverse primer
      forwardprimer = as.character(tolower(toString(a[[1]][1])))#forward primer
      
      stringsets <-DNAStringSet(c(toString(dna$`1`))) # whole Coronavirus genome
      
      fwalignment = pairwiseAlignment(stringsets, DNAStringSet(toString(forwardprimer)), type = "local-global")  # alignment
      start_pos = fwalignment@pattern@range@start # start position
      rvalignment = pairwiseAlignment(stringsets, DNAStringSet(toString(reverseprimer)), type = "local-global")  # alignment
      
      end_pos = rvalignment@pattern@range@start + rvalignment@subject@range@width - 1  # end position
      
      ref = dna$`1`[(start_pos):(end_pos)] # reference DNA with certain positions
      alg = paste(forwardprimer, sep = "", collapse = NULL)#forward primer
      alg1 = paste(reverseprimer, sep = "", collapse = NULL)#reverse primer
      
      if (primers[s, 7] != "NA") {
        
        
        if(length(a[[1]])<3){
          set = NULL
          set = c(set, paste(ref, collapse = ""))
          set = c(set, paste(alg, collapse = ""))
          set = c(set, paste(alg1, collapse = ""))
          set = DNAStringSet(set)
          set@ranges@NAMES = c("Reference", "Forward primer", "Reverse primer")
         
          msaR:: msaR(DECIPHER::AlignSeqs(set), menu = FALSE) # multiple sequence alignment
        
        }
        else if(length(a[[1]])==3){
          forwardprobe = as.character(tolower(toString(a[[1]][3])))#forward primer
          patterns <- c("c3-spacer|fam|bbq|bhq|thf|iabkfg|bqh|iabkfq|rox|zen|mgb|ibfq|tamra|quasar|iabrqsp|tex|hex|labkfq|-|[pdu]")
          reverset= str_replace_all(forwardprobe, patterns, '')
          
          reverset = gsub("[0-9]*", '',  reverset) 
          reverset = gsub("\\W", '',  reverset)  # remove not word
          
          probe1 = paste(reverset, sep = "", collapse = NULL)
          set = NULL
          set = c(set, paste(ref, collapse = ""))
          set = c(set, paste(alg, collapse = ""))
          set = c(set, paste(alg1, collapse = ""))
          set = c(set, paste(probe1, collapse = ""))
          set = DNAStringSet(set)
          set@ranges@NAMES = c("Reference", "Forward primer", "Reverse primer","Probe")
          msaR::msaR(DECIPHER::AlignSeqs(set), menu = FALSE) # multiple sequence alignment
         
        }
        else if(length(a[[1]])==4){
          
          if( s==88){
            
            
            reverseprimer2 = as.character(reverseComplement(DNAString(toString(a[[1]][4])))) #reverse primer
            forwardprimer2 = as.character(tolower(toString(a[[1]][3])))#forward primer
            
            stringsets <-DNAStringSet(c(toString(dna$`1`))) # whole Coronavirus genome
            
            fwalignment = pairwiseAlignment(stringsets, DNAStringSet(toString(forwardprimer2)), type = "local-global")  # alignment
            start_pos = fwalignment@pattern@range@start # start position
            rvalignment = pairwiseAlignment(stringsets, DNAStringSet(toString(reverseprimer2)), type = "local-global")  # alignment
            
            end_pos = rvalignment@pattern@range@start + rvalignment@subject@range@width - 1  # end position
            
            ref2 = dna$`1`[(start_pos):(end_pos)] # reference DNA with certain positions
            alg2 = paste(forwardprimer2, sep = "", collapse = NULL)#forward primer
            alg12 = paste(reverseprimer2, sep = "", collapse = NULL)#reverse primer
            
            
            
            set = NULL
            set = c(set, paste(c(ref,DNAString("NNNNNNNN"),ref2), collapse = ""))
            set = c(set, paste(alg, collapse = ""))
            set = c(set, paste(alg1, collapse = ""))
            
            set = c(set, paste(alg2, collapse = ""))
            set = c(set, paste(alg12, collapse = ""))
            set = DNAStringSet(set)
            cc=c("Reference","Forward","Reverse","Forward-del","Reverse-del")
            
            set@ranges@NAMES =cc
            msaR::msaR(DECIPHER::AlignSeqs(set), menu = FALSE) # multiple sequence alignment
            
            
          }
          else{
            probe1 = as.character(tolower(toString(a[[1]][3])))#probe 1
            patterns <- c("c3-spacer|fam|bbq|bhq|thf|iabkfg|bqh|iabkfq|rox|zen|mgb|ibfq|tamra|quasar|iabrqsp|tex|hex|labkfq|[pdu]")
            reverset= str_replace_all(probe1, patterns, '')
            reverset = gsub("[0-9]*", '',  reverset)
            reverset = gsub("\\W", '',  reverset)  # remove not word
            
            probe1 = paste(reverset, sep = "", collapse = NULL)
            probe2 = as.character(tolower(toString(a[[1]][4])))#probe 1
            reverset= str_replace_all(probe2, patterns, '')
            reverset = gsub("[0-9]*", '',  reverset) 
            reverset = gsub("\\W", '',  reverset)  # remove not word
            
            probe2 = paste(reverseComplement(DNAString(reverset)), sep = "", collapse = NULL)
            
            set = NULL
            set = c(set, paste(ref, collapse = ""))
            set = c(set, paste(alg, collapse = ""))
            set = c(set, paste(alg1, collapse = ""))
            set = c(set, paste(probe1, collapse = ""))
            set = c(set, paste(probe2, collapse = ""))
            set = DNAStringSet(set)
            cc=strsplit(as.character(primers[s, 4]), split = ",\r\n") #forward primer
            
            set@ranges@NAMES =c("Reference",unlist(cc))
            msaR::msaR(DECIPHER::AlignSeqs(set), menu = FALSE) # multiple sequence alignment
          }
          
          
          
          
          
        }  
      }
    }
  })
  
  output$publishedMSA <-  msaR::renderMsaR({
    
   msa1plot()
    
    
    
    
    
  }) # MSA result
  
  
  ranges <- reactiveValues(a1x = NULL, a1y = NULL)
  observeEvent(input$A1dblclick, {
    brush <- input$A1brush
    if (!is.null(brush)) {
      ranges$a1x <- c(brush$xmin, brush$xmax)
      ranges$a1y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$a1x <- NULL
      ranges$a1y <- NULL
    }
  })
   
  # graph for amplicon mutations
  output$publishedprimerplot <- renderPlot({
    
    shiny::validate(shiny::need(input$publishedprimerTable_rows_selected[1], "Select a primer from the table."))
    # plot for amplicon showing mutation places
    primers <- read_excel("publishedprimers_15Dec.xlsx")
    
    primers = as.data.frame(primers)
    s = input$publishedprimerTable_rows_selected[1]
    
    if (primers[s, 7] != "NA") {
      # Check if position of primers are given
      a = strsplit(as.character(primers[s, 7]), split = "\r\n") # primer position
      f = strsplit(as.character(a[[1]][1]), split = "-") #"13814" "13837"
      ff = strsplit(as.character(f[[1]][1]), split = "\t") # forward primer position
      ff = strsplit(as.character(ff), split = "\r")
      f2 = strsplit(as.character(f[[1]][2]), split = "\t")
      f2 = strsplit(as.character(f2), split = "\r")
      fstartpos = as.numeric(ff) # forward start
      fendpos = as.numeric(f2) # forward end
      
      r = strsplit(as.character(a[[1]][2]), split = "-") # reverse primer position
      rr = strsplit(as.character(r[[1]][1]), split = "\t") # reverse primer start position
      rr = strsplit(as.character(rr), split = "\r")
      r2 = strsplit(as.character(r[[1]][2]), split = "\t")
      r2 = strsplit(as.character(r2), split = "\r") # reverse primer end position
      
      rstartpos = as.numeric(rr) # reverse start
      rendpos = as.numeric(r2) # reverse end
      fromfw_torv = fstartpos:rendpos # amplicon
      
      y =as.numeric(allpos[fromfw_torv, 2]) # select the frequencies of amplicon positions
      al = as.data.frame(cbind(fromfw_torv, y))  # first column: positions of amplicon, second: frequencies
      referencenucleotides = as.data.frame(read.csv("reference nucleotides", row.names = 1)) # reading reference nucleotides
      rn = referencenucleotides[(fendpos+1):(rstartpos-1),] # ref nt in region between primers
      refforw <- referencenucleotides[fstartpos:fendpos,] # ref nt of forward primer
      refreverse <- referencenucleotides[rstartpos:rendpos,] # ref nt of reverse primer
      
      # for iupac codes on plot
      mutnexts = merge(allpos, data.frame(fromfw_torv), by.x = "X1", by.y = "fromfw_torv") # to find degenerate nt
      
      colnames(mutnexts) <- c("fromfw_torv","y" ,"iupac")
      # isolate count is added to column 3, it will be used to scale color on plot
      
      al <- mutnexts
      al[,"V4"]=0
      # V4 column of al contains ""(mut count=0), 0(mut count>=955),1(mut count<955)
      al[which(al[,2]>0),4] <- 2
      al[which(al[,2]<0.1),4] <- 1 # do 1 what's below the threshold 
      al[which(al[,2]==0),4] <-3
      titleofgraph=paste("Visualizing Primer Binding Region and Amplicon","\n",primers[s, 3],"_",fstartpos,"-",rendpos )
      
      
      ggplot2::ggplot(al, aes(x = fromfw_torv, y = y)) +
        ggplot2:: geom_segment(aes(x = fromfw_torv,xend = fromfw_torv,y = 0, yend = y, linetype="dotted")) + ggplot2::scale_linetype_identity() +
        ggplot2::geom_point(aes(color = (factor(V4)))) +
        ggplot2::scale_color_manual(values = c("orange","red","white")) +
        
        ggplot2::annotate("text",label = rn,x = ((fendpos+1):(rstartpos-1)) + 0.05,y = -0.001,size = 4,colour = "darkgreen" ,fontface =2) +
        ggplot2::annotate("text",label = refforw,x = (fstartpos:fendpos) + 0.05,y = -0.001,size = 4,colour = "blue" ,fontface =2) +
        ggplot2::annotate("text",label = refreverse,x = (rstartpos:rendpos) + 0.05,y = -0.001,size = 4,colour = "blue",fontface =2 ) +
        ggplot2::geom_text(
          data = subset(al, al[, 3]!= ""),
          aes(label = as.character(subset(al, al[, 3]!= "")[,3])),
          hjust = 0,
          vjust = 0)+ ggplot2::coord_cartesian(xlim = ranges$a1x,
                        ylim = ranges$a1y,
                        expand = FALSE) +
        ggplot2::xlim(fstartpos-0.2, rendpos +0.2)+
        ggplot2::ylim(-0.2, (max(allpos[fromfw_torv, 2])+1.2))+
        ggplot2::xlab("Position") +
        ggplot2::ylab("Percentages of Mutations")+ggtitle(titleofgraph) +
        ggplot2:: theme_light() +
        ggplot2::theme(plot.title = element_text(color="red", size=14, face="bold",hjust=0.5),
              legend.position = "none",
              panel.grid.major.x = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank())
      
      
    }
    
    
  })
  
  output$click_info <- DT::renderDataTable({
    shiny::validate( shiny::need(input$publishedprimerTable_rows_selected[1], "Select a primer from the table."))
    primers <- read_excel("publishedprimers_15Dec.xlsx")
    
    primers = as.data.frame(primers)
    s = input$publishedprimerTable_rows_selected[1]
   
     if (primers[s, 7] != "NA") {
      # Check if position of primers are given
      a = strsplit(as.character(primers[s, 7]), split = "\r\n") # primer position
      f = strsplit(as.character(a[[1]][1]), split = "-") #"13814" "13837"
      ff = strsplit(as.character(f[[1]][1]), split = "\t") # forward primer position
      ff = strsplit(as.character(ff), split = "\r")
      f2 = strsplit(as.character(f[[1]][2]), split = "\t")
      f2 = strsplit(as.character(f2), split = "\r")
      fstartpos = as.numeric(ff) # forward start position
      fendpos = as.numeric(f2)
      
      r = strsplit(as.character(a[[1]][2]), split = "-") # reverse primer position
      rr = strsplit(as.character(r[[1]][1]), split = "\t") # reverse primer start position
      rr = strsplit(as.character(rr), split = "\r")
      r2 = strsplit(as.character(r[[1]][2]), split = "\t")
      r2 = strsplit(as.character(r2), split = "\r") # reverse primer end position
      
     
      rstartpos = as.numeric(rr)
      rendpos = as.numeric(r2) #reverse end position
      
      n <- as.vector(fstartpos: rendpos)
    
    }
    Mutationdata <- data.frame(Mutationdata)
   # m<- data.frame(Mutationdata[is.element(Mutationdata[,1], n), ] )
    # m<- data.frame(Mutationdata[Mutationdata[,1] %in% n, ] )
    
    #print(m)
    m <- data.frame(Mutationdata  %>% dplyr::filter(position%in% n))
    DT::datatable(data.frame(m), caption = 'Table: Total sequence number is 5,345,994',
                  options = list(dom = "Blfrtip", buttons = list("copy","csv"),
                                 paging = TRUE,searching =TRUE,lengthMenu = c(10, 30, 50),
                                                                               pageLength = 5,
                                                                               scrollX = TRUE))

    
  })
  output$publishedphylop <- renderPlot({
    
    shiny::validate(
      shiny::need(input$publishedprimerTable_rows_selected[1], "Select a primer from the table.")
    )
    # plot for amplicon showing mutation places
    primers <- read_excel("publishedprimers_15Dec.xlsx")
    
    primers = as.data.frame(primers)
    s = input$publishedprimerTable_rows_selected[1]
    
    if (primers[s, 7] != "NA") {
      # Check if position of primers are given
      a = strsplit(as.character(primers[s, 7]), split = "\r\n") # primer position
      
      f = strsplit(as.character(a[[1]][1]), split = "-") #"13814" "13837"
      
      ff = strsplit(as.character(f[[1]][1]), split = "\t") # forward primer position
      ff = strsplit(as.character(ff), split = "\r")
      f2 = strsplit(as.character(f[[1]][2]), split = "\t")
      f2 = strsplit(as.character(f2), split = "\r")
      fstartpos = as.numeric(ff)
      fendpos = as.numeric(f2)
      
      r = strsplit(as.character(a[[1]][2]), split = "-") # reverse primer position
      
      rr = strsplit(as.character(r[[1]][1]), split = "\t") # reverse primer start position
      rr = strsplit(as.character(rr), split = "\r")
      r2 = strsplit(as.character(r[[1]][2]), split = "\t")
      r2 = strsplit(as.character(r2), split = "\r") # reverse primer end position
      
      
      rstartpos = as.numeric(rr)
      rendpos = as.numeric(r2)
      
      a <- fstartpos
      b <- fstartpos:rendpos
      
      PhyloP <- read.delim("PhyloP", header = FALSE)
      PhyloP[, 2] = as.double(specify_decimal(PhyloP[, 2], 2))
      
      
      new=(PhyloP[b,2])
      new=as.data.frame(cbind(b,new))
      colnames(new) <- c("pos", "p")
      
      
      ggplot2::ggplot(data=new, aes(x=pos, y=p)) + ggplot2::geom_segment(aes(x = pos,xend = pos,y = 0,yend = p)) +
        ggplot2::geom_bar(stat="identity", width=1,fill="steelblue")+ ggplot2::ylim((min(new[,2])-1), (max(new[, 2])+1))+
        ggplot2::geom_text(aes(label=p), vjust=1.6, color="black", size=4, angle=15)+ ggplot2::coord_cartesian(xlim = ranges$a1x,
                                                                                             ylim = ranges$a1y,
                                                                                             expand = FALSE) +
        ggplot2::xlab("Position") +
        ggplot2::ylab("PhyloP values")+
        ggplot2::theme_minimal()
     
    }
    
    
  })
  
  
  
  ##########  ##########  ##########  ##########  ##########  ##########
  
  # ----                                        Conserved primer table
  
  ranges <- reactiveValues(st = NULL, en = NULL)
  
  observeEvent(input$do2, {
    range_values <- switch (input$conservedRegion,
                            "all" = c(266, 29903),
                            "nsp1" = c(266, 805),
                            "nsp2" = c(806, 2719),
                            "nsp3" = c(2720, 8554),
                            "nsp4" = c(8555, 10054),
                            "3C-like proteinase" = c(10055, 10972),
                            "nsp6" = c(10973, 11842),
                            "nsp7" = c(11843, 12091),
                            "nsp8" = c(12092, 12685),
                            "nsp9" = c(12686, 13024),
                            "nsp10" = c(13025, 13441),
                            "RNA-dependent RNA polymerase" = c(13442, 16236),
                            "helicase" = c(16237, 18039),
                            "3'-to-5' exonuclease" = c(18040, 19620),
                            "endoRNAse" = c(19621, 20658),
                            "2'-O-ribose methyltransferase" = c(20659, 21552),
                            "surface" = c(21563, 25384),
                            "ORF3" = c(25393, 26220),
                            "envelope" = c(26245, 26472),
                            "membrane" = c(26523, 27191),
                            "ORF6" = c(27202, 27387),
                            "ORF7a" = c(27394, 27759),
                            "ORF7b" = c(27756, 27887),
                            "ORF8" = c(27894, 28259),
                            "capsid" = c(28274, 29533),
                            "ORF10" = c(29558, 29674)
    )
    ranges$st <-as.numeric(range_values[1])
    ranges$en <- as.numeric(range_values[2])
    
  })
  
  
  Conservedprimers <- eventReactive (input$do2, {
    # new primer table
    Conservedprimers = read_csv("Conservedprimers_15Dec.csv")
    Conservedprimers <- Conservedprimers[which(Conservedprimers[, 3] >= ranges$st), ]
    Conservedprimers <- Conservedprimers[which(Conservedprimers[, 10] <= ranges$en), ]
    Conservedprimers[,16] <- round(Conservedprimers[,16],3)
    Conservedprimers[,17] <- round(Conservedprimers[,17],3)
    return(as.data.frame(Conservedprimers))
  })
  output$conserverdprimerTable <- DT::renderDataTable ({
    shiny::validate(shiny::need(input$do2, "Click the 'visualize' button to see the table of primer pairs."))
    
    consprim =Conservedprimers()
    
    colnames(consprim)=c("Forward Primer","Fw Length", "Fw Start", "Fw End", "Fw Tm", "Fw GC",
                         "Reverse Primer","Rv Length", "Rv Start", "Rv End", "Rv Tm", "Rv GC",
                         "Tm Difference", "Amplicon Length", "Amplicon Mutations","Fw SARS-CoV-2 Specificity","Rv SARS-CoV-2 Specificity")
    colnames(consprim)[c(1, 7)] = c(
      '<span style="color:red">Forward Primer</span>',
      '<span style="color:red">Reverse Primer</span>'
    )
    
    consprim[, 1] = toupper(consprim[, 1])
    consprim[, 7] = toupper(consprim[, 7])
    
    datatable(
      consprim,
      rownames = FALSE,
      escape = FALSE,
      caption = 'Potential conserved primer pairs in your gene of interest. You can use the boxes below column names to filter.',
      filter = "top",
      selection = "single",
      options = list(initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#B2B6B9', 'color': 'fff'});",
        "}"), columnDefs=list(list(className='dt-center',targets="_all")),
        headerCallback = DT::JS(
          "function(thead) {",
          "  $(thead).css('font-size', '1.2em');",
          "}"
        ),
        language = list(zeroRecords = "No primer pairs found to display"),
        lengthMenu = c(3,10, 30, 50),
        pageLength = 3,
        colnames=c("Forward Primer","Fw Length", "Fw Start", "Fw End", "Fw Tm", "Fw GC",
                   "Reverse Primer","Rv Length", "Rv Start", "Rv End", "Rv Tm", "Rv GC",
                   "Tm Difference", "Amplicon Length", "Fw SARS-CoV-2 Specificity","Rv SARS-CoV-2 Specificity"),
        scrollX = TRUE
      )
    ) %>%
      formatStyle(
        '<span style="color:red">Forward Primer</span>',
        fontWeight = 'bold',
        color = 'red'
      ) %>%
      formatStyle(
        '<span style="color:red">Reverse Primer</span>',
        fontWeight = 'bold',
        color = 'red'
      )
    
    
  })
  
  ranges <- reactiveValues(cx = NULL, cy = NULL)
  observeEvent(input$condblclick, {
    brush <- input$conbrush
    if (!is.null(brush)) {
      ranges$cx <- c(brush$xmin, brush$xmax)
      ranges$cy <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$cx <- NULL
      ranges$cy <- NULL
    }
  })
  
  output$conservedplot <- renderPlot({
    
    shiny::validate(shiny::need(input$conserverdprimerTable_rows_selected[1], "Select a primer from the table."))
    
  
    Conservedprimers =data.frame(Conservedprimers())
    ss = input$conserverdprimerTable_rows_selected[1]
    
    fstartpos = Conservedprimers[ss, 3]
    fendpos=Conservedprimers[ss, 4]
    rstartpos = Conservedprimers[ss, 9]
    rendpos = Conservedprimers[ss, 10]
    fromfw_torv = fstartpos:rendpos# from forward start to reverse end
    
    y =as.numeric(allpos[fromfw_torv, 2]) # select the frequencies of amplicon positions
    al = as.data.frame(cbind(data.frame(fromfw_torv), y))  # first column: positions of amplicon, secon: frequencies
    
    referencenucleotides = as.data.frame(read.csv("reference nucleotides", row.names = 1)) # reading reference nucleotides
    rn = referencenucleotides[(fendpos+1):(rstartpos-1),] # ref nt in region between primers
    refforw <- referencenucleotides[fstartpos:fendpos,] # ref nt of forward primer
    refreverse <- referencenucleotides[rstartpos:rendpos,] # ref nt of reverse primer
 
    
    # for iupac codes on plot
    mutnexts = merge(allpos, data.frame(fromfw_torv), by.x = "X1", by.y = "fromfw_torv") # to find degenerate nt
    
    colnames(mutnexts) <- c("fromfw_torv","y" ,"iupac")
    # isolate count is added to column 3, it will be used to scale color on plot
    
    al <- mutnexts
    al[,"V4"]=0
    # V4 column of al contains ""(mut count=0), 0(mut count>=955),1(mut count<955)
    al[which(al[,2]>0),4] <- 2
    al[which(al[,2]<0.1),4] <- 1 
    al[which(al[,2]==0),4] <-3
    
    
    ggplot2::ggplot(al, aes(x = fromfw_torv, y = y)) +
      ggplot2::geom_segment(aes(x = fromfw_torv,xend = fromfw_torv,y = 0,yend = y,linetype="dotted")) +
      ggplot2::scale_linetype_identity() +
      
      ggplot2::geom_point(aes(color = (factor(V4)))) +
      ggplot2::scale_color_manual(values = c("orange","red","white")) +
      
      ggplot2::annotate("text",label = rn,x = ((fendpos+1):(rstartpos-1)) + 0.05,y = -0.001,size = 4,colour = "darkgreen",fontface =2 ) +
      ggplot2::annotate("text",label = refforw,x = (fstartpos:fendpos) + 0.05,y = -0.001,size = 4,colour = "blue" ,fontface =2) +
      ggplot2::annotate("text",label = refreverse,x = (rstartpos:rendpos) + 0.05,y = -0.001,size = 4,colour = "blue",fontface =2 ) +
      
      ggplot2::geom_text(
        data = subset(al, al[, 3]!= ""),
        aes(label = as.character(subset(al, al[, 3]!= "")[,3])),
         colour="black",size=4,
        hjust = -0.5,position=position_dodge(width = .75)) +
      ggplot2::coord_cartesian(xlim = ranges$cx, ylim = ranges$cy,expand = FALSE) +
      ggplot2::xlab("Position") +
      ggplot2::ylab("Percentages of Mutations")+
      ggplot2::ggtitle("Visualizing Primer Binding Region and Amplicon") +
      ggplot2::xlim(fstartpos-0.2, rendpos +0.2)+
      ggplot2:: ylim(-0.2, max(allpos[fromfw_torv, 2])+0.4 )+
      ggplot2::theme_light() +
      ggplot2::theme(plot.title = element_text(color="red", size=14, face="bold",hjust=0.5),
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank()
      )
    
    
  })
  
  output$conservedphylop <- renderPlot({
    
    shiny::validate(shiny::need(input$conserverdprimerTable_rows_selected[1], "Select a primer from the table."))
    
  
    Conservedprimers =Conservedprimers()
    selected = input$conserverdprimerTable_rows_selected[1]
    
    a <- Conservedprimers[selected, 3]
    b <- a:Conservedprimers[selected, 10]
    
    PhyloP <- read.delim("PhyloP", header = FALSE)
    PhyloP[, 2] = as.double(specify_decimal(PhyloP[, 2], 2))
    
    
    new=(PhyloP[b,2])
    new=as.data.frame(cbind(b,new))
    colnames(new) <- c("pos", "p")
    
    
    ggplot2::ggplot(data=new, aes(x=pos, y=p)) + ggplot2::geom_segment(aes(x = pos,xend = pos,y = 0,yend = p)) +
      ggplot2::geom_bar(stat="identity", width=1,fill="steelblue")+ ggplot2::ylim((min(new[,2])-1), (max(new[, 2])+1))+
      ggplot2::geom_text(aes(label=p), vjust=1.6, color="black", size=4, angle=15)+ ggplot2::coord_cartesian(xlim = ranges$cx,
                                                                                           ylim = ranges$cy,
                                                                                           expand = FALSE) +
      ggplot2:: xlab("Position") +
      ggplot2::ylab("PhyloP values")+
      ggplot2::theme_minimal()
    
    
    
  })
  msa2plot <- eventReactive(input$msa2, {
    
    if (input$conserverdprimerTable_rows_selected[1]) {
      Conservedprimers =Conservedprimers()
      selected = input$conserverdprimerTable_rows_selected[1]
      
      reverseprimer = as.character(reverseComplement(DNAString(toString(Conservedprimers[selected, 7])))) #reverse primer
      forwardprimer = as.character(tolower(toString(Conservedprimers[selected, 1])))#forward primer
      
      f = paste(forwardprimer, sep = "", collapse = NULL)#forward primer
      r = paste(reverseprimer, sep = "", collapse = NULL)#reverse primer
      
      stringsets <- DNAStringSet(c(toString(dna$`1`))) # whole Coronavirus genome
      
      
      fwalignment = pairwiseAlignment(stringsets, DNAStringSet(toString(forwardprimer)), type = "local-global")  # alignment
      start_pos = fwalignment@pattern@range@start # start position
      
      rvalignment = pairwiseAlignment(stringsets,DNAStringSet(toString(reverseprimer)) , type = "local-global")  # alignment
      
      end_pos = rvalignment@pattern@range@start + rvalignment@subject@range@width - 1  # end position
      
      ref = dna$`1`[(start_pos):(end_pos)] # reference DNA with certain positions
      
      set = NULL
      set = c(set, paste(ref, collapse = ""))
      set = c(set, paste(f, collapse = ""))
      set = c(set, paste(r, collapse = ""))        
      set = DNAStringSet(set)
      set@ranges@NAMES = c("Reference", "Forward primer", "Reverse primer")
      msaR(AlignSeqs(set), menu = FALSE) # multiple sequence alignment
      
    }
    
  })
  output$conservedMSA <- renderMsaR({
    
  

    msa2plot()
    
    
    
    
  }) 
  mutationinformation_consreactive <- reactive({
    Conservedprimers =Conservedprimers()
    
    ss = input$conserverdprimerTable_rows_selected[1]
    
    fstartpos = Conservedprimers[ss, 3]
    rendpos = Conservedprimers[ss, 10]
    fromfw_torv = fstartpos:rendpos
    
    Mutationdata <- data.frame(Mutationdata)
    
    mutationinformation_consreactive<- data.frame(Mutationdata  %>% dplyr::filter(position%in% fromfw_torv))
  })
  output$mutationinformation_cons <- DT::renderDataTable({
    
    shiny::validate(shiny::need(input$conserverdprimerTable_rows_selected[1], "Select a primer from the table."))
    
  
    DT::datatable(data.frame(mutationinformation_consreactive()),  caption = 'Table: Total sequence number is 6,854,423',extensions = 'Buttons', 
                  options = list(dom = "Blfrtip", buttons = c("copy","csv"), 
                                 paging = TRUE,searching =TRUE,lengthMenu = c(10, 30, 50),
                                                                                     pageLength = 5,
                                                                                     scrollX = TRUE))
    
    
    
  })
  
  ##########  ##########  ##########  ##########  ##########  ##########
  
  
  # ----                                        Degenerate primer table
  
  rangesD <- reactiveValues(st = NULL, en = NULL)
  observeEvent(input$do4, {
    
    range_values <- switch (input$degenerateRegion,
                            "all" = c(266, 29903),
                            "nsp1" = c(266, 805),
                            "nsp2" = c(806, 2719),
                            "nsp3" = c(2720, 8554),
                            "nsp4" = c(8555, 10054),
                            "3C-like proteinase" = c(10055, 10972),
                            "nsp6" = c(10973, 11842),
                            "nsp7" = c(11843, 12091),
                            "nsp8" = c(12092, 12685),
                            "nsp9" = c(12686, 13024),
                            "nsp10" = c(13025, 13441),
                            "RNA-dependent RNA polymerase" = c(13442, 16236),
                            "helicase" = c(16237, 18039),
                            "3'-to-5' exonuclease" = c(18040, 19620),
                            "endoRNAse" = c(19621, 20658),
                            "2'-O-ribose methyltransferase" = c(20659, 21552),
                            "surface" = c(21563, 25384),
                            "ORF3" = c(25393, 26220),
                            "envelope" = c(26245, 26472),
                            "membrane" = c(26523, 27191),
                            "ORF6" = c(27202, 27387),
                            "ORF7a" = c(27394, 27759),
                            "ORF7b" = c(27756, 27887),
                            "ORF8" = c(27894, 28259),
                            "capsid" = c(28274, 29533),
                            "ORF10" = c(29558, 29674)
    )
    
    rangesD$st <- range_values[1]
    rangesD$en <- range_values[2]
    
    
  })
  # degenerate primer
  DegeneratedprimersX <- eventReactive (input$do4, {
    
    # new primer table
    
    degeneratedprimers <-data.frame(read_csv("Degenerateprimers_15Dec.csv"))
    
    degeneratedprimers <- degeneratedprimers[which(degeneratedprimers[, 3] >= rangesD$st), ]
    degeneratedprimers <- degeneratedprimers[which(degeneratedprimers[, 10] <= rangesD$en), ]

    return(as.data.frame(degeneratedprimers))
  })
  
  
  Degeneratedprimers <- reactive({
    req(input$do4)
    degeneratedprimers = DegeneratedprimersX()
    
    colnames(degeneratedprimers)=c("Forward Primer 5'->3'","Fw Length", "Fw Start", "Fw End", "Fw Tm", "Fw GC",
                                   "Reverse Primer  5'->3","Rv Length", "Rv Start", "Rv End", "Rv Tm", "Rv GC",
                                   "Tm Difference", "Amplicon Length","Amplicon Mutations", "Fw SARS-CoV-2 Specificity","Rv SARS-CoV-2 Specificity",
                                   "Degenerate Fw  5'->3", "Fw Mutated Positions", "Fw Lineages","Fw Mutation number", "Fw Mutation Frequency (%)",
                                   "Degenerate Rv  5'->3", "Rv Mutated Positions","Rv Lineages","Rv Mutation number", "Rv Mutation Frequency (%)")
    
    colnames(degeneratedprimers)[c(1, 7)] = c('<span style="color:red">Forward Primer</span>','<span style="color:red">Reverse Primer</span>')
    
    
    if(input$lineage=="all"){
      degeneratedprimers[, 1] = toupper(degeneratedprimers[, 1])
      degeneratedprimers[, 7] = toupper(degeneratedprimers[, 7])
    }
    else {
      degeneratedprimers[, 1] = toupper(degeneratedprimers[, 1])
      degeneratedprimers[, 7] = toupper(degeneratedprimers[, 7])
      
      lin <- input$lineage
      
      degeneratedprimers1 <-  degeneratedprimers[grepl(lin, degeneratedprimers[,20]),]
      degeneratedprimers2 <-  degeneratedprimers[grepl(lin, degeneratedprimers[,25]),]
      
      degeneratedprimers <-as.data.frame(unique(rbind(degeneratedprimers1,degeneratedprimers2)))
    }
    return(degeneratedprimers)
  })
   

  
  output$degenerateprimerTable <- DT::renderDataTable ({
    req(input$do4)
    shiny::validate( shiny::need(input$do4, "Click the 'visualize' button to see the table of primer pairs."))
    
    degeneratedprimers = as.data.frame(Degeneratedprimers())
    colnames(degeneratedprimers)[c(1, 7)] = c('<span style="color:red">Forward Primer</span>','<span style="color:red">Reverse Primer</span>')
    
    
    datatable(
      degeneratedprimers,
      rownames = FALSE,
      escape = FALSE,
      caption = 'Degenerate primer pairs in your gene of interest. You can use the boxes below column names to filter.',
      filter = "top",
      selection = "single",
      options = list(initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#B2B6B9', 'color': 'fff'});",
        "}"), columnDefs=list(list(className='dt-center',targets="_all")),
        headerCallback = DT::JS(
          "function(thead) {",
          "  $(thead).css('font-size', '1.2em');",
          "}"
        ),
        language = list(zeroRecords = "No primer pairs found to display"),
        lengthMenu = c(10, 30, 50),
        pageLength = 3,
        scrollX = TRUE
      )
    ) %>%
      formatStyle(
        '<span style="color:red">Forward Primer</span>',
        fontWeight = 'bold',
        color = 'red'
      ) %>%
      formatStyle(
        '<span style="color:red">Reverse Primer</span>',
        fontWeight = 'bold',
        color = 'red'
      )
    
    
  })
 
  ranges <- reactiveValues(dx = NULL, dy = NULL)
  observeEvent(input$degdblclick, {
    brush <- input$degbrush
    if (!is.null(brush)) {
      ranges$dx <- c(brush$xmin, brush$xmax)
      ranges$dy <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$dx <- NULL
      ranges$dy <- NULL
    }
  })
  
  output$degenerateplot <- renderPlot({
    shiny::validate(shiny::need(input$degenerateprimerTable_rows_selected[1], "Select a primer from the table."))
    
    degeneratedprimers = Degeneratedprimers()
    selected = input$degenerateprimerTable_rows_selected[1]
    
    
    
    fstartpos = degeneratedprimers[selected, 3]
    fendpos=degeneratedprimers[selected, 4]
    rstartpos = degeneratedprimers[selected, 9]
    rendpos = degeneratedprimers[selected, 10]
    
    fromfw_torv = fstartpos:rendpos# from forward start to reverse end
    
    y =as.numeric(allpos[fromfw_torv, 2]) # select the frequencies of amplicon positions
    al = as.data.frame(cbind(data.frame(fromfw_torv), y))  # first column: positions of amplicon, secon: frequencies
    
    referencenucleotides = as.data.frame(read.csv("reference nucleotides", row.names = 1)) # reading reference nucleotides
    rn = referencenucleotides[(fendpos+1):(rstartpos-1),] # ref nt in region between primers
    refforw <- referencenucleotides[fstartpos:fendpos,] # ref nt of forward primer
    refreverse <- referencenucleotides[rstartpos:rendpos,] # ref nt of reverse primer
    
    # for iupac codes on plot
    mutnexts = merge(allpos, data.frame(fromfw_torv), by.x = "X1", by.y = "fromfw_torv") # to find degenerate nt
    
    colnames(mutnexts) <- c("fromfw_torv","y" ,"iupac")
    # isolate count is added to column 3, it will be used to scale color on plot
    
    al <- mutnexts
    al[,"V4"]=0
    # V4 column of al contains ""(mut count=0), 0(mut count>=955),1(mut count<955)
    al[which(al[,2]>0),4] <- 2
    al[which(al[,2]<0.1),4] <- 1 
    al[which(al[,2]==0),4] <-3
    
    ggplot2::ggplot(al, aes(x = fromfw_torv, y = y)) +
      ggplot2::geom_segment(aes(x = fromfw_torv,xend = fromfw_torv,y = 0,yend = y,linetype="dotted")) +
      ggplot2::scale_linetype_identity() +
      
      ggplot2::geom_point(aes(color = (factor(V4)))) +
      ggplot2::scale_color_manual(values = c("orange","red","white")) +
      
      ggplot2::annotate("text",label = rn,x = ((fendpos+1):(rstartpos-1)) + 0.05,y = -0.001,size = 4,colour = "darkgreen",fontface =2 ) +
      ggplot2::annotate("text",label = refforw,x = (fstartpos:fendpos) + 0.05,y = -0.001,size = 4,colour = "blue" ,fontface =2) +
      ggplot2::annotate("text",label = refreverse,x = (rstartpos:rendpos) + 0.05,y = -0.001,size = 4,colour = "blue",fontface =2 ) +
      
      ggplot2::geom_text(
        data = subset(al, al[, 3]!= ""),
        aes(label = as.character(subset(al, al[, 3]!= "")[,3])),
        hjust = 0,
        vjust = 0) +
      ggplot2::coord_cartesian(xlim = ranges$dx, ylim = ranges$dy,expand = FALSE) +
      ggplot2::xlab("Position") +
      ggplot2::ylab("Percentages of Mutations")+
      ggplot2::ggtitle("Visualizing Primer Binding Region and Amplicon") +
      ggplot2::xlim(fstartpos-0.2, rendpos +0.2)+
      ggplot2::ylim(-0.2, max(allpos[fromfw_torv, 2])+0.4 )+
      ggplot2::theme_light() +
      ggplot2::theme(plot.title = element_text(color="red", size=14, face="bold",hjust=0.5),
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank()
      )
    
    
  })
  
  output$degeneratephylop <- renderPlot({
    
    shiny::validate(shiny::need(input$degenerateprimerTable_rows_selected[1], "Select a primer from the table."))
    degeneratedprimers = Degeneratedprimers()
    
    selected = input$degenerateprimerTable_rows_selected[1]
    
    a <- degeneratedprimers[selected, 3]
    b <- a:degeneratedprimers[selected, 10]
    
    PhyloP <- read.delim("PhyloP", header = FALSE)
    PhyloP[, 2] = as.double(specify_decimal(PhyloP[, 2], 2))
    
    new=(PhyloP[b,2])
    new=as.data.frame(cbind(b,new))
    colnames(new) <- c("pos", "p")
    
    ggplot2::ggplot(data=new, aes(x=pos, y=p)) + ggplot2::geom_segment(aes(x = pos,xend = pos,y = 0,yend = p)) +
      ggplot2:: geom_bar(stat="identity", width=1,fill="steelblue")+ ggplot2::ylim((min(new[,2])-1), (max(new[, 2])+1))+
      ggplot2::geom_text(aes(label=p), vjust=1.6, color="black", size=4, angle=15)+ 
      ggplot2::coord_cartesian(xlim = ranges$dx,ylim = ranges$dy, expand = FALSE) +
      ggplot2::xlab("Position") +
      ggplot2:: ylab("PhyloP values")+
      ggplot2::theme_minimal()
    
    
  })
  
  msa3plot <- eventReactive(input$msa3, {
    
    if (input$degenerateprimerTable_rows_selected[1]) {
      
      degeneratedprimers = Degeneratedprimers()
      
      
      selected = input$degenerateprimerTable_rows_selected[1]
      
      reverseprimer = as.character(reverseComplement(DNAString(toString(degeneratedprimers[selected, 7])))) #reverse primer
      forwardprimer = as.character(tolower(toString(degeneratedprimers[selected, 1])))#forward primer
      
      f = paste(forwardprimer, sep = "", collapse = NULL)#forward primer
      r = paste(reverseprimer, sep = "", collapse = NULL)#reverse primer
      
      stringsets <- DNAStringSet(c(toString(dna$`1`))) # whole Coronavirus genome
      
      
      fwalignment = pairwiseAlignment(stringsets, DNAStringSet(toString(forwardprimer)), type = "local-global")  # alignment
      start_pos = fwalignment@pattern@range@start # start position
      
      rvalignment = pairwiseAlignment(stringsets,DNAStringSet(toString(reverseprimer)) , type = "local-global")  # alignment
      
      end_pos = rvalignment@pattern@range@start + rvalignment@subject@range@width - 1  # end position
      
      ref = dna$`1`[(start_pos):(end_pos)] # reference DNA with certain positions
      
      set = NULL
      set = c(set, paste(ref, collapse = ""))
      set = c(set, paste(f, collapse = ""))
      set = c(set, paste(r, collapse = ""))            
      set = DNAStringSet(set)
      set@ranges@NAMES = c("Reference", "Forward primer", "Reverse primer")
      msaR(AlignSeqs(set), menu = FALSE) # multiple sequence alignment
      
    }
    
  })
  output$degenerateMSA <- renderMsaR({
    
   # validate(shiny::need(input$degenerateprimerTable_rows_selected[1], "Select a primer from the table."))
    msa3plot()

    
  }) 
 
   mutationinformation_degreactive <- reactive({
    degeneratedprimers = Degeneratedprimers()
    
    selected = input$degenerateprimerTable_rows_selected[1]
    
    print(degeneratedprimers[selected, 3])
    
    fromfw_torv <- as.integer(degeneratedprimers[selected, 3]):as.integer(degeneratedprimers[selected, 10])
    Mutationdata <- data.frame(Mutationdata)
    
    mutationinformation_degreactive <- data.frame(Mutationdata  %>% dplyr::filter(position%in% fromfw_torv))
  })
  output$mutationinformation_deg <- DT::renderDataTable({
    
    shiny::validate(shiny::need(input$degenerateprimerTable_rows_selected[1], "Select a primer from the table."))
  
    DT::datatable(data.frame(mutationinformation_degreactive()),  caption = 'Table: Total sequence number is 6,854,423',
                  extensions = "Buttons", options = list(dom = "Blfrtip", buttons = c("copy","csv"), 
                                                                        paging = TRUE,searching =TRUE,lengthMenu = c(10, 30, 50),
                                                                                      pageLength = 5,
                                                                                      scrollX = TRUE))
    
  })
  
  
  output$lineagetable <- DT::renderDataTable({
   lineageandmutations <- read_excel("lineageandmutations.xlsx")
      tab <- as.data.frame(lineageandmutations)
    
    return(tab)
    
  }, options = list(pageLength=10, columnDefs = list(list(className = 'dt-center', targets = 3:13)))
  )
  
  output$lineagetable2 <- DT::renderDataTable({
    lineageandmutations <- read_excel("lineageandmutations.xlsx")
    tab <- as.data.frame(lineageandmutations)
    
    return(tab)
    
  }, options = list(pageLength=10 ,columnDefs = list(list(className = 'dt-center', targets = 3:13)))
  )
  output$lineagetable3 <- DT::renderDataTable({
    lineageandmutations <- read_excel("lineageandmutations.xlsx")
    tab <- as.data.frame(lineageandmutations)
    
    return(tab)
    
  }, options = list(pageLength=10 ,columnDefs = list(list(className = 'dt-center', targets = 3:13)))
  )
  
  ##########  ##########  ##########  ##########  ##########  ##########
  
  
  # ----                                      Variable
  
  rangesV <- reactiveValues(st = NULL, en = NULL)
  observeEvent(input$do5, {
    
    range_values <- switch (input$variableregion,
                            "all" = c(266, 29903),
                            "nsp1" = c(266, 805),
                            "nsp2" = c(806, 2719),
                            "nsp3" = c(2720, 8554),
                            "nsp4" = c(8555, 10054),
                            "3C-like proteinase" = c(10055, 10972),
                            "nsp6" = c(10973, 11842),
                            "nsp7" = c(11843, 12091),
                            "nsp8" = c(12092, 12685),
                            "nsp9" = c(12686, 13024),
                            "nsp10" = c(13025, 13441),
                            "RNA-dependent RNA polymerase" = c(13442, 16236),
                            "helicase" = c(16237, 18039),
                            "3'-to-5' exonuclease" = c(18040, 19620),
                            "endoRNAse" = c(19621, 20658),
                            "2'-O-ribose methyltransferase" = c(20659, 21552),
                            "surface" = c(21563, 25384),
                            "ORF3" = c(25393, 26220),
                            "envelope" = c(26245, 26472),
                            "membrane" = c(26523, 27191),
                            "ORF6" = c(27202, 27387),
                            "ORF7a" = c(27394, 27759),
                            "ORF7b" = c(27756, 27887),
                            "ORF8" = c(27894, 28259),
                            "capsid" = c(28274, 29533),
                            "ORF10" = c(29558, 29674)
    )
    
    rangesV$st <- range_values[1]
    rangesV$en <- range_values[2]
    
    
  })
  # VARIABLE primer
  Variableprimers <- eventReactive (input$do5, {
    
    # new primer table
    
    variableprimers <-data.frame(read_csv("ConsprimerforHighlyVariableregions_15Dec.csv"))
   
    variableprimers <- variableprimers[which(variableprimers[, 3] >= rangesV$st), ]
    variableprimers <- variableprimers[which(variableprimers[, 10] <= rangesV$en), ]
    
    return(as.data.frame(variableprimers))
  })
  
  output$variableprimerTable <- DT::renderDataTable ({
    
    shiny::validate( shiny::need(input$do5, "Click the 'visualize' button to see the table of primer pairs."))
    
    variableprimers = Variableprimers()
    colnames(variableprimers)=c("Forward Primer","Fw Length", "Fw Start", "Fw End", "Fw Tm", "Fw GC",
                                "Reverse Primer","Rv Length", "Rv Start", "Rv End", "Rv Tm", "Rv GC",
                                "Tm Difference", "Amplicon Length", "Amplicon Mutations","Fw SARS-CoV-2 Specificity","Rv SARS-CoV-2 Specificity")
    
    colnames(variableprimers)[c(1, 7)] = c('<span style="color:red">Forward Primer</span>','<span style="color:red">Reverse Primer</span>')
    
    
    
    
    datatable(
      variableprimers,
      rownames = FALSE,
      escape = FALSE,
      caption = 'Conserved primer pairs in your gene of interest. You can use the boxes below column names to filter.',
      filter = "top",
      selection = "single",
      options = list(initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': '#B2B6B9', 'color': 'fff'});",
        "}"), columnDefs=list(list(className='dt-center',targets="_all")),
        headerCallback = DT::JS(
          "function(thead) {",
          "  $(thead).css('font-size', '1.2em');",
          "}"
        ),
        language = list(zeroRecords = "No primer pairs found to display"),
        lengthMenu = c(10, 30, 50),
        pageLength = 3,
        scrollX = TRUE
      )
    ) %>%
      formatStyle(
        '<span style="color:red">Forward Primer</span>',
        fontWeight = 'bold',
        color = 'red'
      ) %>%
      formatStyle(
        '<span style="color:red">Reverse Primer</span>',
        fontWeight = 'bold',
        color = 'red'
      )
    
    
  })
  
  ranges <- reactiveValues(vx = NULL, vy = NULL)
  observeEvent(input$variabledblclick, {
    brush <- input$variablebrush
    if (!is.null(brush)) {
      ranges$vx <- c(brush$xmin, brush$xmax)
      ranges$vy <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$vx <- NULL
      ranges$vy <- NULL
    }
  })
  
  output$variableteplot <- renderPlot({
    shiny::validate(shiny::need(input$variableprimerTable_rows_selected[1], "Select a primer from the table."))
    
    variableprimers = Variableprimers()
    selected = input$variableprimerTable_rows_selected[1]
    
    
    
    fstartpos = variableprimers[selected, 3]
    fendpos=variableprimers[selected, 4]
    rstartpos = variableprimers[selected, 9]
    rendpos = variableprimers[selected, 10]
    
    fromfw_torv = fstartpos:rendpos# from forward start to reverse end
    
    y =as.numeric(allpos[fromfw_torv, 2]) # select the frequencies of amplicon positions
    al = as.data.frame(cbind(data.frame(fromfw_torv), y))  # first column: positions of amplicon, secon: frequencies
    
    referencenucleotides = as.data.frame(read.csv("reference nucleotides", row.names = 1)) # reading reference nucleotides
    rn = referencenucleotides[(fendpos+1):(rstartpos-1),] # ref nt in region between primers
    refforw <- referencenucleotides[fstartpos:fendpos,] # ref nt of forward primer
    refreverse <- referencenucleotides[rstartpos:rendpos,] # ref nt of reverse primer
    
   
    
    # for iupac codes on plot
    mutnexts = merge(allpos, data.frame(fromfw_torv), by.x = "X1", by.y = "fromfw_torv") # to find degenerate nt
    
    colnames(mutnexts) <- c("fromfw_torv","y" ,"iupac")
    # isolate count is added to column 3, it will be used to scale color on plot
    
    al <- mutnexts
    al[,"V4"]=0
    # V4 column of al contains ""(mut count=0), 0(mut count>=955),1(mut count<955)
    al[which(al[,2]>0),4] <- 2
    al[which(al[,2]<0.1),4] <- 1 # thresholdun altında kalanları 1 yap
    al[which(al[,2]==0),4] <-3
    
    
    
    ggplot2::ggplot(al, aes(x = fromfw_torv, y = y)) +
      ggplot2::geom_segment(aes(x = fromfw_torv,xend = fromfw_torv,y = 0,yend = y,linetype="dotted")) +
      ggplot2::scale_linetype_identity() +
      
      ggplot2::geom_point(aes(color = (factor(V4)))) +
      ggplot2::scale_color_manual(values = c("orange","red","white")) +
      
      ggplot2::annotate("text",label = rn,x = ((fendpos+1):(rstartpos-1)) + 0.05,y = -0.001,size = 4,colour = "darkgreen",fontface =2 ) +
      ggplot2::annotate("text",label = refforw,x = (fstartpos:fendpos) + 0.05,y = -0.001,size = 4,colour = "blue" ,fontface =2) +
      ggplot2::annotate("text",label = refreverse,x = (rstartpos:rendpos) + 0.05,y = -0.001,size = 4,colour = "blue",fontface =2 ) +
      
      ggplot2::geom_text(
        data = subset(al, al[, 3]!= ""),
        aes(label = as.character(subset(al, al[, 3]!= "")[,3])),colour="black",size=4,
        hjust = -0.5,position=position_dodge(width = .75)) +
      ggplot2::coord_cartesian(xlim = ranges$vx, ylim = ranges$vy,expand = FALSE) +
      ggplot2::xlab("Position") +
      ggplot2::ylab("Percentages of Mutations")+
      ggplot2::ggtitle("Visualizing Primer Binding Region and Amplicon") +
      ggplot2::xlim(fstartpos-0.2, rendpos +0.2)+
      ggplot2::ylim(-0.2, max(allpos[fromfw_torv, 2])+0.4 )+
      ggplot2::theme_light() +
      ggplot2::theme(plot.title = element_text(color="red", size=14, face="bold",hjust=0.5),
            legend.position = "none",
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank()
      )
    
  })
  
  output$variablephylop <- renderPlot({
    
    shiny::validate(shiny::need(input$variableprimerTable_rows_selected[1], "Select a primer from the table."))
    variableprimers = Variableprimers()
    
    selected = input$variableprimerTable_rows_selected[1]
    
    a <- variableprimers[selected, 3]
    b <- a:variableprimers[selected, 10]
    
    PhyloP <- read.delim("PhyloP", header = FALSE)
    PhyloP[, 2] = as.double(specify_decimal(PhyloP[, 2], 2))
    
    new=(PhyloP[b,2])
    new=as.data.frame(cbind(b,new))
    colnames(new) <- c("pos", "p")
    
    ggplot2::ggplot(data=new, aes(x=pos, y=p)) +  
      ggplot2::geom_segment(aes(x = pos,xend = pos,y = 0,yend = p)) +
      ggplot2::geom_bar(stat="identity", width=1,fill="steelblue")+  
      ggplot2::ylim((min(new[,2])-1), (max(new[, 2])+1))+
      ggplot2::geom_text(aes(label=p), vjust=1.6, color="black", size=4, angle=15)+  
      ggplot2::coord_cartesian(xlim = ranges$vx,ylim = ranges$vy, expand = FALSE) +
      ggplot2::xlab("Position") +
      ggplot2::ylab("PhyloP values")+
      ggplot2::theme_minimal()
    
  })
  msaVplot <- eventReactive(input$msaV, {
    
    if (input$variableprimerTable_rows_selected[1]) {
      
      variableprimers =Variableprimers()
      
      
      selected = input$variableprimerTable_rows_selected[1]
      
      reverseprimer = as.character(reverseComplement(DNAString(toString(variableprimers[selected, 7])))) #reverse primer
      forwardprimer = as.character(tolower(toString(variableprimers[selected, 1])))#forward primer
      
      f = paste(forwardprimer, sep = "", collapse = NULL)#forward primer
      r = paste(reverseprimer, sep = "", collapse = NULL)#reverse primer
      
      stringsets <- DNAStringSet(c(toString(dna$`1`))) # whole Coronavirus genome
      
      
      fwalignment = pairwiseAlignment(stringsets, DNAStringSet(toString(forwardprimer)), type = "local-global")  # alignment
      start_pos = fwalignment@pattern@range@start # start position
      
      rvalignment = pairwiseAlignment(stringsets,DNAStringSet(toString(reverseprimer)) , type = "local-global")  # alignment
      
      end_pos = rvalignment@pattern@range@start + rvalignment@subject@range@width - 1  # end position
      
      ref = dna$`1`[(start_pos):(end_pos)] # reference DNA with certain positions
      
      set = NULL
      set = c(set, paste(ref, collapse = ""))
      set = c(set, paste(f, collapse = ""))
      set = c(set, paste(r, collapse = ""))            
      set = DNAStringSet(set)
      set@ranges@NAMES = c("Reference", "Forward primer", "Reverse primer")
      msaR(AlignSeqs(set), menu = FALSE) # multiple sequence alignment
      
    }
    
  })
  output$variableMSA <- renderMsaR({
    
    # shiny::validate(shiny::need(input$degenerateprimerTable_rows_selected[1], "Select a primer from the table."))
    msaVplot()
    
    
  }) 
  
  
  mutationinformation_varreactive <- reactive({
    variableprimers =Variableprimers()
    selected = input$variableprimerTable_rows_selected[1]
    
    fromfw_torv <- as.integer(variableprimers[selected, 3]):as.integer(variableprimers[selected, 10])
    Mutationdata <- data.frame(Mutationdata)
    
    mutationinformation_varreactive <- data.frame(Mutationdata  %>% dplyr::filter(position%in% fromfw_torv))
  })
  output$mutationinformation_var <- DT::renderDataTable({
    
    shiny::validate(shiny::need(input$variableprimerTable_rows_selected[1], "Select a primer from the table."))
   
    DT::datatable(data.frame(mutationinformation_varreactive()),  caption = 'Total sequence number is 6,854,423',
                  extensions= "Buttons", options = list(dom = "Blfrtip", buttons = c("copy","csv"), 
                                                paging = TRUE,searching =TRUE,lengthMenu = c(10, 30, 50),
                                                pageLength = 5,
                                                scrollX = TRUE))
    
  })
  
  
  ##########  ##########  ##########  ##########  ##########  ##########
  
  
  # ----                                       Search tab
  
  
  forward <- reactive({
    input$do3
    codes=c("A","C","G","T","M","R","W","S","Y","K"," V","H","D","B","N")
    noncodes=LETTERS[!LETTERS %in% codes]
    Pattern = paste(noncodes, collapse="|")
    Pattern1 = paste(tolower(noncodes), collapse="|")
    
    shiny::validate(shiny::need( grepl(Pattern, input$forwardtext)==FALSE, "Remove these letters from your sequence: E, F, I, J, L, O, P, Q, U, V, X, Z"))
    shiny::validate(shiny::need( grepl(Pattern1, input$forwardtext)==FALSE, "Remove these letters from your sequence: E, F, I, J, L, O, P, Q, U, V, X, Z"))
    shiny::validate(shiny::need( nchar(input$forwardtext)>10, "Enter longer sequence."))
    stringsets <- DNAStringSet(c(toString(dna$`1`)))
    # remove gaps if there is 
    forwardprimer <- input$forwardtext
    searchString <- ' '  # remove gaps
    replacementString <- ''
    forwardprimer = gsub(searchString, replacementString, forwardprimer)
    
    forwardprimer = as.character(tolower(forwardprimer)) # lower letter
    
    # alignment with reference to find the start/end positions of primer
    fwalignment = pairwiseAlignment(
      stringsets,
      DNAStringSet(toString(forwardprimer)),
      gapOpening = -2,
      gapExtension = -8,
      scoreOnly = FALSE,
      type = "local-global")  # alignment
    
    start_pos = fwalignment@pattern@range@start #start position
    end_pos = fwalignment@pattern@range@start + fwalignment@subject@range@width - 1 # end position
    l = fwalignment@subject@range@width #length
    
    deg_forw = forw_degenerate(forwardprimer, start_pos, end_pos) # find degenerate forward primer
    
    f = data[(start_pos):(end_pos), 2]
    
    mut_forw <- length(which(f != 0)) # number of mutations
    if(mut_forw==0){
      max_Tm_forw <- as.double(calculate_tm(forwardprimer))
      min_Tm_forw <- max_Tm_forw 
      max_Gc_forw <- as.double(calculate_GC(forwardprimer))
      min_Gc_forw <- max_Gc_forw
    }
    else{
      max_Tm_forw <- max(as.double(find_Tm(deg_forw)))
      min_Tm_forw <- min(as.double(find_Tm(deg_forw)))
  
      
      max_Gc_forw <- max(as.double(find_Gc(deg_forw)))
      min_Gc_forw <- min(as.double(find_Gc(deg_forw)))
    }
    # mutation positions/ changes
    forw_mut_pos = paste(unlist(as.character(zt[start_pos:end_pos, 2])), collapse =" ")
    # mutation freq
    mutFreq = paste(unlist(as.character(zt[start_pos:end_pos, 3])), collapse = "/")
    gg <- data.frame(Mutationdata  %>% dplyr::filter(position%in% start_pos:end_pos))[,17:19]
    gg2 <- as.vector(gg[!is.na(gg)])
    lineages <- paste( unique( gg2[gg2!=""]),collapse = " , ")
    forward <- c(
      deg_forw,
      l,
      start_pos,
      end_pos,
      mut_forw,
      round(max_Tm_forw,2),
       round(min_Tm_forw,2),
       round(max_Gc_forw,2),
      round(min_Gc_forw,2),
      forw_mut_pos,
      mutFreq,
      lineages
    )
  })
  
  reversepr <- reactive({
   input$do3
    shiny::validate(shiny::need( nchar(input$reversetext)>10, "Enter longer sequence."))
    codes=c("A","C","G","T","M","R","W","S","Y","K"," V","H","D","B","N")
    noncodes=LETTERS[!LETTERS %in% codes]
    Pattern = paste(noncodes, collapse="|")
    Pattern1 = paste(tolower(noncodes), collapse="|")
    
    
    shiny::validate(shiny::need( grepl(Pattern, input$reversetext)==FALSE, "Remove these letters from your sequence: E, F, I, J, L, O, P, Q, U, V, X, Z"))
    shiny::validate(shiny::need( grepl(Pattern1, input$reversetext)==FALSE, "Remove these letters from your sequence: E, F, I, J, L, O, P, Q, U, V, X, Z"))
    
    shiny::validate(shiny::need(input$reversetext != input$forwardtext, "Please enter different forward or reverse pair."))
    
    reverseprimer <- input$reversetext
    searchString <- ' ' #remove gaps
    replacementString <- ''
    reverseprimer = gsub(searchString, replacementString, reverseprimer)
    reverseprimer = as.character(tolower(reverseprimer))
    # take reverse to make it 3` --- 5`
    
    stringsets <- DNAStringSet(c(toString(dna$`1`)))
    
    rvalignment = pairwiseAlignment(
      stringsets,
      reverseComplement(DNAString(reverseprimer)),
      gapOpening = -2,
      gapExtension = -8,
      scoreOnly = FALSE,
      type = "local-global"
    )
    start_pos = rvalignment@pattern@range@start
    end_pos = rvalignment@pattern@range@start + rvalignment@subject@range@width - 1
    l = rvalignment@subject@range@width
    
    deg_rev = (reverse_degenerate(tolower(as.character(reverseprimer)), start_pos, end_pos))# degenerate reverse primer
    deg_rev <- gsub("A", "a", deg_rev)
    deg_rev <- gsub("T", "t", deg_rev)
    deg_rev <- gsub("C", "c", deg_rev)
    deg_rev <- gsub("G", "g", deg_rev)
    f = data[(start_pos):(end_pos), 2]
    
    mut_rev <- length(which(f != 0)) # number of mutations in reverse
    if(mut_rev==0){
      max_Tm_rev <- as.double(calculate_tm(reverseprimer))
      min_Tm_rev <- max_Tm_rev 
      max_Gc_rev <- as.double(calculate_GC(reverseprimer))
      min_Gc_rev <- max_Gc_rev
    }
    else{
      max_Tm_rev <- max(as.double(find_Tm(deg_rev)))
      min_Tm_rev <- min(as.double(find_Tm(deg_rev)))
      
      max_Gc_rev <- max(as.double(find_Gc(deg_rev)))
      min_Gc_rev <- min(as.double(find_Gc(deg_rev)))
    }
   
    # mutation frequencies
    rev_mut_pos = reverse(paste(unlist(as.character(zt[start_pos:end_pos, 2])), collapse =" "))
    mutFreq = (paste(unlist(as.character(zt[end_pos:start_pos, 3])), collapse ="/"))
    gg <- data.frame(Mutationdata  %>% dplyr::filter(position%in% start_pos:end_pos))[,17:19]
    gg2 <- as.vector(gg[!is.na(gg)])
    lineages <- paste( unique( gg2[gg2!=""]),collapse = " , ")
    reversepr <-
      c(
        deg_rev,
        l,
        start_pos,end_pos,
        mut_rev,
        round(max_Tm_rev,2),
        round(min_Tm_rev,2) ,
        round(max_Gc_rev,2),
        round(min_Gc_rev,2),
        rev_mut_pos,
        mutFreq,
        lineages
      )
    
  })
  
  probe <- reactive({
    input$do3
    if (input$pr=="No") {
      return(NULL)
    }
    else {
      codes=c("A","C","G","T","M","R","W","S","Y","K"," V","H","D","B","N")
      noncodes=LETTERS[!LETTERS %in% codes]
      Pattern = paste(noncodes, collapse="|")
      Pattern1 = paste(tolower(noncodes), collapse="|")
      
      shiny::validate(shiny::need( grepl(Pattern, input$probetext)==FALSE, "Remove these letters from your sequence: E, F, I, J, L, O, P, Q, U, V, X, Z"))
      shiny::validate(shiny::need( grepl(Pattern1, input$probetext)==FALSE, "Remove these letters from your sequence: E, F, I, J, L, O, P, Q, U, V, X, Z"))
      shiny::validate(shiny::need( nchar(input$probetext)>10, "Enter longer sequence."))
      stringsets <- DNAStringSet(c(toString(dna$`1`)))
      # remove gaps if there is 
      probe <- input$probetext
      searchString <- ' '  # remove gaps
      replacementString <- ''
      probe = gsub(searchString, replacementString, probe)
      
      probe = as.character(tolower(probe)) # lower letter
      
      # alignment with reference to find the start/end positions of primer
      probealignment = pairwiseAlignment(
        stringsets,
        DNAStringSet(toString(probe)),
        gapOpening = -2,
        gapExtension = -8,
        scoreOnly = FALSE,
        type = "local-global")  # alignment
      
      start_pos = probealignment@pattern@range@start #start position
      end_pos = probealignment@pattern@range@start + probealignment@subject@range@width - 1 # end position
      l = probealignment@subject@range@width #length
      
      deg_forw = forw_degenerate(probe, start_pos, end_pos) # find degenerate forward primer
      
      f = data[(start_pos):(end_pos), 2]
      
      mut_forw <- length(which(f != 0)) # number of mutations
      if(mut_forw==0){
        max_Tm_forw <- as.double(calculate_tm(probe))
        min_Tm_forw <- max_Tm_forw 
        max_Gc_forw <- as.double(calculate_GC(probe))
        min_Gc_forw <- max_Gc_forw
      }
      else{
        max_Tm_forw <- max(as.double(find_Tm(deg_forw)))
        min_Tm_forw <- min(as.double(find_Tm(deg_forw)))
        
        max_Gc_forw <- max(as.double(find_Gc(deg_forw)))
        min_Gc_forw <- min(as.double(find_Gc(deg_forw)))
      }
      # mutation positions/ changes
      forw_mut_pos = paste(unlist(as.character(zt[start_pos:end_pos, 2])), collapse =" ")
      # mutation freq
      mutFreq = paste(unlist(as.character(zt[start_pos:end_pos, 3])), collapse = "/")
      gg <- data.frame(Mutationdata  %>% dplyr::filter(position%in% start_pos:end_pos))[,17:19]
      gg2 <- as.vector(gg[!is.na(gg)])
      lineages <- paste( unique( gg2[gg2!=""]),collapse = " , ")
      probe <- c(
        deg_forw,
        l,
        start_pos,
        end_pos,
        mut_forw,
        round(max_Tm_forw,2),
        round(min_Tm_forw,2),
        round(max_Gc_forw,2),
        round(min_Gc_forw,2),
        forw_mut_pos,
        mutFreq,
        lineages
      )
    }
  
  })
  output$searchyourprimersTable <- DT::renderDataTable ({
    
    shiny::validate(shiny::need(input$forwardtext, "Enter a valid sequence for your forward primer."))
    
    shiny::validate(shiny::need(input$reversetext, "Enter a valid sequence for your reverse primer."))
    shiny::validate(shiny::need( input$do3, "Click the 'Search' button."))
   
    if (input$pr=="No"){
      tabl <- cbind(forward(), reversepr())
      
      amp_leng <- as.numeric(reversepr()[4]) - as.numeric(forward()[3])+1
     
      tabl = rbind(tabl, amp_leng)
      rownames(tabl) <-
        c(
          "Primers",
          "Primer Length",
          "Start position of Primer (3'end)",
          "End Position of Primer (5'end)",
          "Number of Mutations",
          "Max Tm",
          "Min Tm",
          "Max GC",
          "Min GC",
          "Position of Mutations (5'->3')",
          "Mutation percentages (5'->3')",
          "Lineages having mutations",
          "Amplicon Length"
        )
      colnames(tabl) <- c("Forward Primer", "Reverse Primer")
      
    }
    else{ 
      probe <- probe()
    
    tabl <- cbind(forward(), reversepr(), probe)
    
    amp_leng <- as.numeric(reversepr()[4]) - as.numeric(forward()[3])+1
    probe_reverse <- as.numeric(reversepr()[4]) - as.numeric(probe[3])+1
    a <- c(amp_leng,"",probe_reverse)
    tabl = rbind(tabl, a)
    rownames(tabl) <-
      c(
        "Primers",
        "Primer Length",
        "Start position of Primer (3'end)",
        "End Position of Primer (5'end)",
        "Number of Mutations",
        "Max Tm",
        "Min Tm",
        "Max GC",
        "Min GC",
        "Position of Mutations (5'->3')",
        "Mutation percentages (5'->3')",
        "Lineages having mutations",
        "Amplicon Length"
      )
    colnames(tabl) <- c("Forward Primer", "Reverse Primer","Probe")
    }
    datatable(
      tabl,
      rownames = TRUE,
      options = list(
        initComplete = JS(
          "function(settings, json) {",
          "$(this.api().table().header()).css({'background-color': '#B2B6B9', 'color': 'fff'});",
          "}"
        ), 
        pageLength = 30,
        scrollX = TRUE,
        paging = FALSE
      )
    )
    DT::datatable(data.frame(tabl),  extensions = 'Buttons', options = list(dom = "Blfrtip", buttons = c("copy","csv"
    ),paging = FALSE,searching =TRUE,lengthMenu = FALSE,
    pageLength = 30,
    scrollX = TRUE))
    
  })
  
  ranges <- reactiveValues(a2x = NULL, a2y = NULL)
  
  observeEvent(input$sdblclick, {
    brush <- input$sbrush
    if (!is.null(brush)) {
      ranges$a2x <- c(brush$xmin, brush$xmax)
      ranges$a2y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$a2x <- NULL
      ranges$a2y <- NULL
    }
  })
  
  output$searchplot <- renderPlot({
    shiny::validate(shiny::need(input$forwardtext, "Enter a valid sequence for your primers.") )
    
    shiny::validate(shiny::need(input$reversetext, "Enter a valid sequence for your primers."))
    if (input$align){
      
      fstartpos = as.numeric(forward()[3])
      rendpos = as.numeric(reversepr()[4])
      fendpos = as.numeric(forward()[4])
      rstartpos = as.numeric(reversepr()[3])
      pstartpos = as.numeric(probe()[3])
      pendpos = as.numeric(probe()[4])
      
      pos <- c(fstartpos,rstartpos,pstartpos,fendpos,rendpos,pendpos)
     
      fromfw_torv = min(pos):max(pos)# from forward start to reverse end

      
      y =as.numeric(allpos[fromfw_torv, 2]) # select the frequencies of amplicon positions
      al = as.data.frame(cbind(data.frame(fromfw_torv), y))  # first column: positions of amplicon, second: frequencies
      
      referencenucleotides = as.data.frame(read.csv("reference nucleotides", row.names = 1)) # reading reference nucleotides
      rn = referencenucleotides[(min(pos)):(max(pos)),] # ref nt in region between primers
      refforw <- referencenucleotides[fstartpos:fendpos,] # ref nt of forward primer
      refreverse <- referencenucleotides[rstartpos:rendpos,] # ref nt of reverse primer
      refprobe <- referencenucleotides[pstartpos:pendpos,] # ref nt of probe
      
      mutnexts = merge(allpos, data.frame(fromfw_torv), by.x = "X1", by.y = "fromfw_torv") # to find degenerate nt
      
      colnames(mutnexts) <- c("fromfw_torv","y" ,"iupac")
      # isolate count is added to column 3, it will be used to scale color on plot
      
      al <- mutnexts
      al[,"V4"]=0
      # V4 column of al contains ""(mut count=0), 0(mut count>=955),1(mut count<955)
      al[which(al[,2]>0),4] <- 2
      al[which(al[,2]<0.1),4] <- 1 # thresholdun altında kalanları 1 yap
      al[which(al[,2]==0),4] <-3
      
      
      ggplot2::ggplot(al, aes(x = fromfw_torv, y = y)) +
        ggplot2::geom_segment(aes(x = fromfw_torv,xend = fromfw_torv,y = 0,yend = y,linetype="dotted")) +
        ggplot2::scale_linetype_identity() +
        
        ggplot2::geom_point(aes(color = (factor(V4)))) +
        ggplot2::scale_color_manual(values = c("orange","red","white")) +
       
        ggplot2::annotate("text",label = rn,x = fromfw_torv + 0.05, y = -0.001,size = 4,colour = "black",fontface =2 ) +
        ggplot2::annotate("text",label = refforw,x = (fstartpos:fendpos) + 0.05,y = -0.001,size = 4,colour = "blue" ,fontface =2) +
        ggplot2::annotate("text",label = refreverse,x = (rstartpos:rendpos) + 0.05,y = -0.001,size = 4,colour = "blue",fontface =2 ) +
        ggplot2::annotate("text",label = refprobe,x = (pstartpos:pendpos) + 0.05,y = -0.001,size = 4,colour = "darkgreen",fontface =2 ) +
        
        ggplot2::geom_text(
          data = subset(al, al[, 3]!= ""),
          aes(label = as.character(subset(al, al[, 3]!= "")[,3])),colour="black",size=4,
          hjust = -0.5,position=position_dodge(width = .75)
         ) +
        ggplot2::coord_cartesian(xlim = ranges$a2x, ylim = ranges$a2y,expand = FALSE) +
        ggplot2::xlab("Position") +
        ggplot2::ylab("Percentages of Mutations")+
        ggplot2:: ggtitle("Visualizing Primer Binding Region and Amplicon") +
        ggplot2::xlim(fstartpos-0.2, rendpos +0.2)+
        ggplot2::ylim(-0.2, max(allpos[fromfw_torv, 2])+0.4 )+
        ggplot2::theme_light() +
        ggplot2:: theme(plot.title = element_text(color="red", size=14, face="bold",hjust=0.5),
              legend.position = "none",
              panel.grid.major.x = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank())
      
      
      
      
    }

    
  })
 
  
  output$searchphylop<- renderPlot({
    
    shiny::validate(shiny::need(input$forwardtext, "Enter a valid sequence for your primers.") )
    
    shiny::validate(shiny::need(input$reversetext, "Enter a valid sequence for your primers."))
    
    if (input$align){
      
      fstartpos = as.numeric(forward()[3])
      rendpos = as.numeric(reversepr()[4])
      fendpos = as.numeric(forward()[4])
      rstartpos = as.numeric(reversepr()[3])
      pstartpos = as.numeric(probe()[3])
      pendpos = as.numeric(probe()[4])
      
      pos <- c(fstartpos,rstartpos,pstartpos,fendpos,rendpos,pendpos)
      
      b  = min(pos):max(pos)# from forward start to reverse end
      
     
      
      PhyloP <- read.delim("PhyloP", header = FALSE)
      PhyloP[, 2] = as.double(specify_decimal(PhyloP[, 2], 2))
      
      
      
      new=(PhyloP[b,2])
      new=as.data.frame(cbind(b,new))
      colnames(new) <- c("pos", "p")
      
      ggplot2::ggplot(data=new, aes(x=pos, y=p)) + ggplot2::geom_segment(aes(x = pos,xend = pos,y = 0,yend = p)) +
        ggplot2::geom_bar(stat="identity", width=1,fill="steelblue")+ggplot2::ylim((min(new[,2])-1), (max(new[, 2])+1))+
        ggplot2::geom_text(aes(label=p), vjust=1.6, color="black", size=4,angle=15)+ ggplot2::coord_cartesian(xlim = ranges$a2x,
                                                                                            ylim = ranges$a2y,
                                                                                            expand = FALSE) +
        ggplot2::xlab("Position") +
        ggplot2::ylab("PhyloP values")+
        ggplot2::theme_minimal()
      
      
      
    }
    
    
  })
  
  msa4plot <- eventReactive(input$msa4, {
    
    if (input$align) {
      reverset <- input$reversetext
      forwardt <- input$forwardtext
      probedt <- input$probetext
      searchString <- ' '
      replacementString <- ''
      reverset = gsub(searchString, replacementString,  reverset)  # remove gaps in text
      forwardt = gsub(searchString, replacementString, forwardt)
      probedt = gsub(searchString, replacementString, probedt)
      
      forwardprimer = as.character(tolower(forwardt))
      probe = as.character(tolower(probedt))
      reverseprimer = as.character(reverseComplement(DNAString(reverset))) #reverse primer
      f = paste(forwardprimer,  sep = "", collapse = NULL)
      r = paste(reverseprimer,  sep = "", collapse = NULL)
      p = paste(probe,  sep = "", collapse = NULL)
      
      
      fstartpos = as.numeric(forward()[3])
      rendpos = as.numeric(reversepr()[4])
      fendpos = as.numeric(forward()[4])
      rstartpos = as.numeric(reversepr()[3])
      pstartpos = as.numeric(probe()[3])
      pendpos = as.numeric(probe()[4])
      
      pos <- c(fstartpos,rstartpos,pstartpos,fendpos,rendpos,pendpos)
      
      fromfw_torv = min(pos):max(pos)# from forward start to reverse end
      
      ref = dna$`1`[fromfw_torv] # reference DNA with certain positions
      
      set = NULL
      set = c(set, paste(ref, collapse = ""))
      set = c(set, paste(f, collapse = ""))
      set = c(set, paste(r, collapse = ""))    
      set = c(set, paste(p, collapse = "")) 
      set = DNAStringSet(set)
      set@ranges@NAMES = c("Reference", "Forward Primer", "Reverse Primer","Probe")
      msaR::msaR(AlignSeqs(set), menu = F)
      
      
    
      
      
      
      
      
     # AlignSeqs(set)
      #set1 <- AlignSeqs(set)
      #msaR(as.sequences(set1))
      
     # msaR(DNAMultipleAlignment(set1))
    }
    
  })
  output$searchMSA <- renderMsaR({
    
    shiny::validate(shiny::need(input$forwardtext, "Enter a valid sequence for your primers.") )
    
    shiny::validate(shiny::need(input$reversetext, "Enter a valid sequence for your primers."))
    msa4plot()
  })
  
 mutationinfor_search <- reactive({
   if (input$align){
     if (input$mutinfo=="f2r"){
       fstartpos = as.numeric(forward()[3])
     }
     else {
       fstartpos = as.numeric(probe()[3])
     }
     
     rendpos = as.numeric(reversepr()[4])
     fromfw_torv = fstartpos:rendpos
     Mutationdata <- data.frame(Mutationdata)
     
     m <- data.frame(Mutationdata  %>% dplyr::filter(position%in% fromfw_torv))
     
     mutationinfor_search <-  m
     
   }
 })
  output$mutationinformation_search <- DT::renderDataTable({
    
    shiny::validate(shiny::need(input$forwardtext, "Enter a valid sequence for your primers.") )
    
    shiny::validate(shiny::need(input$reversetext, "Enter a valid sequence for your primers."))
    
      DT::datatable(data.frame(mutationinfor_search()),
                    caption = 'Table: Total sequence number is 6,854,423',
                    extensions = 'Buttons', options = list(dom = "Blfrtip", buttons = c("copy","csv")
                                                                           ,paging = TRUE,searching =TRUE,lengthMenu = c(10, 30, 50),
                                                                                               pageLength = 5,
                                                                                               scrollX = TRUE))
    
  })

  
  output$down_mut_cons <- downloadHandler(
    filename = function() {
      paste0("Detailed Variation Information For Conserved", ".tsv")
    },
    content = function(file) {
      vroom::vroom_write(mutationinformation_consreactive(), file)
    }
  )
  
  output$down_mut_deg <- downloadHandler(
    filename = function() {
      paste0("Detailed Variation Information For Degenerate", ".tsv")
    },
    content = function(file) {
      vroom::vroom_write(mutationinformation_degreactive(), file)
    }
  )
  
  
  output$down_mut_var <- downloadHandler(
    filename = function() {
      paste0("Detailed Variation Information", ".tsv")
    },
    content = function(file) {
      vroom::vroom_write(mutationinformation_varreactive(), file)
    }
  )
  
  
   output$down_mut_search <- downloadHandler(
    filename = function() {
      paste0(input$mutinfo, ".tsv")
    },
    content = function(file) {
      vroom::vroom_write(mutationinfor_search(), file)
    }
  )
  # conserved table
  output$conserverdprimerTabledownload = downloadHandler(
    filename = "selected_conserved_primers.tsv",
    content = function(file) {
      vroom::vroom_write(Conservedprimers(), file)
    }
  )
  # degenerate table
  output$degenerateprimerTabledownload = downloadHandler(
    filename = "selected_degenerate_primers.tsv",
    content = function(file) {
      vroom::vroom_write(Degeneratedprimers(), file)
    }
  )
  
  # variable table
  output$variableprimerTabledownload = downloadHandler(
    filename = "selected_conserved_primers_for_variableregions.tsv",
    content = function(file) {
      vroom::vroom_write(Variableprimers(), file)
    }
  )
  

  
  
  ##########  ##########  ##########  ##########  ##########  ##########
  

  Mutationdata <- read.csv("MutationDatawithlineages15Dec.csv",header = T)
  Mutationdata <-  Mutationdata %>% mutate_if(is.numeric, ~round(., 3))
  # freq is count/total,  not percentage.
 # Mutationdata <- Mutationdata[,-1]
  colnames(  Mutationdata)[1] <- "position"
  
  ### Most used codes in app----------------------------------------------------------
  dbConn <- dbConnect(SQLite(), ":memory:") # create database
  Seqs2DB(
    "reference_genome.fasta",
    type = "FASTA",
    dbFile = dbConn,
    identifier = "reference"
  ) 
  
  dna <- SearchDB(dbConn)
  dna1 = strsplit(as.character(dna), split = "")
  xx = as.integer(1:29903)
  data2 = as.data.frame(cbind(xx, dna1)) # positions and reference nt
  
  # all positions
  allp <-  data.frame(read_csv("allpositions-iupac-change-percentDec15.csv"))

  #allp[,6] <- round((allp[,6]),5)
  allpos = cbind(as.integer(1:29903), as.integer(rep(0, 29903)))
  allpos[allp[, 1], 2] = allp[, 6]
  allpos <- data.frame(allpos)

  allpos[, 3] = ""
  allpos[allp[, 1], 3] = unlist(allp[, 4])
  
  
  #  only positions >0.1
  nexts <-  data.frame(read_csv("Dec15_0.1threshold-iupac-change-percent.csv"))
  
  nexts[,6] <- round((nexts[,6]),3)
  x = as.integer(1:29903)
  y = as.integer(rep(0, 29903))
  data = cbind(x, y)
  data[nexts[, 1], 2] = nexts[, 6]
  

  
  zt <- data.frame(1:29903,0)
  zt[nexts[, 1], 2] = nexts[, 5] # number of mutation changes in nt
  zt = as.data.frame(zt)
  zt[, 3] = 0
  zt[nexts[, 1], 3] = nexts[, 6] # mutation frequencies
  
}

