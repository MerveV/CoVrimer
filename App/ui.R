
library(shiny)
library(BiocManager)
library(dashboardthemes)
library(shinydashboard)
library(DT)
library(msaR)
library(msa)
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

library(shinycssloaders)
library(plotly)
source("calculate_GC.R")
source("calculate_tm.R")
source("find_Tm.R") 
source("find_Gc.R") 
source("forw_degenerate.R")
source("reverse_degenerate.R") 
source("specify_decimal.R")

library(shinyjs)
ui <- dashboardPage(
  skin='red',
  #-----
  dashboardHeader(titleWidth = 300,title = "CoVrimer",
                  dropdownMenu(type = "notifications",
                               notificationItem( text = "Last update of data: 2022-09-14",icon = icon("calendar"))),
                  dropdownMenu(type = "messages",icon = icon("at"),
                               messageItem(from = "Contact information",message = tags$div("merve.vural@bilkent.edu.tr" , tags$br(),"konu@fen.bilkent.edu.tr",style = "display: inline-block; vertical-align: middle;")))
  ),
  dashboardSidebar( width = 300, 
                    sidebarMenu(id = "sider",
                                menuItem(h4("Documentation"),tabName = "introduction",icon = icon("book")),
                                menuItem(h4("Align Published Primer Sets"), tabName = "publishedprimers",icon = icon("random")),
                                menuItem(h4("Select New Primer Pairs"),tabName = "design",icon = icon("flask"),
                                         menuSubItem('Select Conserved/Degenerate Primers', tabName = "select_tab", icon = icon('hand-pointer')),
                                         menuSubItem('Get Degenerate Primer Pairs', tabName = "get_tab", icon = icon("search") ))
                    )
  ),
  
  #-----
  
  dashboardBody(  useShinyjs(debug = TRUE),
  tags$head(tags$style(HTML('
      .main-header .logo {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 44px;
      }
    '))),tags$head(tags$style(HTML('
                               /* body */
                               .content-wrapper, .right-side {
                               background-color: #fff;
                               }

                               '))),
                ###----
                tabItems(
                  
                  tabItem(tabName = "introduction",  img(src = "CoVrimer.svg", align = "center"),
                          p(strong(h2("How to use the CoVrimer app:"))),
                          includeMarkdown("mdfiles/intro.md"),
                          
                          actionButton("go", h4("Click here for more details")),
                          uiOutput("details"),
                  ),
                  tabItem(tabName = "publishedprimers",
                          fluidRow(column(8, align="center", offset = 2,h3(strong("Table of published primer sets")), 
                                          h4("Please select a primer pair to display mutation and conservation plots for the amplicon and alignment to the reference genome")),
                                   style="text-align:justify;color:white;background-color:#535557;padding:1px;border-radius:10px" ),
                          
                              fluidRow(actionButton("description1", "Click for Description",icon("info")),
                                       bsModal("modal", "", "description1", size = "large",includeMarkdown("mdfiles/descriptionforpublishedprimers.md"))
                                       ,DT::dataTableOutput("publishedprimerTable"),
                                       downloadButton("publishedprimerdownload", "Download")) ,
                       
                          fluidRow(box(title="Detailed Variation Information",width = 12,collapsible = TRUE,solidHeader = TRUE,
                                       DT::dataTableOutput("click_info"))),
                          
                          fluidRow( includeMarkdown("mdfiles/plotexplain.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"),
                          
                          fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                               color = getOption("spinner.color", default = "#0275D8"),
                                               div(plotOutput("publishedprimerplot",dblclick = "A1dblclick",
                                                          brush = brushOpts(id = "A1brush",resetOnNew = TRUE)),
                                   
                          
                          
                        plotOutput("publishedphylop",dblclick = "A1dblclick",brush = brushOpts(id = "A1brush",resetOnNew = TRUE)))))
                          ,
                          fluidRow(actionButton("msa1", "Click for MSA result"),
                           box(width = 12,msaROutput("publishedMSA", width = "100%", height = "100%")))
                  ),
                  # Design primers Conserved/Degenerate ----
                  tabItem(tabName = "select_tab",
                          #####
                          tabsetPanel(
                            tabPanel(h4("Conserved Primer pairs"),
                                     fluidRow(tags$style(HTML("
     .tabbable > .nav > li > a                  {background-color: #B2B6B9;  color:black}
    
    .tabbable > .nav > li[class=active]    > a {background-color: #404040; color:white}
  ")),
                                              tags$style(HTML("
.box.box-solid.box-danger>.box-header {
  color:#fff;
  background:#606060 }
.box.box-solid.box-danger{
border-bottom-color:#606060;
border-left-color:#606060;
border-right-color:#606060;
border-top-color:#606060;
}

.box.box-danger>.box-header {
  color:#000000;
  background:#fff
  
                    }")),
                                              
                                              
                                              box(title=h4("Select the primer pair to display mutation and Phylop plots for the amplicon  and alignment to the reference genome"),
                                                  status = "danger", solidHeader = TRUE,width = 12,collapsible = TRUE,height="2em",
                                                  
                                                  column(width = 6, selectInput("conservedRegion",h5("Select the region that you want to amplify"),
                                                                                c( "All"="all", "nsp1" = "nsp1",
                                                                                   "nsp2" = "nsp2",
                                                                                   "nsp3" = "nsp3",
                                                                                   "nsp4" = "nsp4",
                                                                                   "3C-like proteinase" = "3C-like proteinase",
                                                                                   "nsp6" = "nsp6",
                                                                                   "nsp7" ="nsp7",
                                                                                   "nsp8" = "nsp8",
                                                                                   "nsp9" = "nsp9",
                                                                                   "nsp10" = "nsp10",
                                                                                   "RNA-dependent RNA polymerase" = "RNA-dependent RNA polymerase",
                                                                                   "Helicase" = "helicase",
                                                                                   "3'-to-5' exonuclease" = "3'-to-5' exonuclease",
                                                                                   "endoRNAse" ="endoRNAse",
                                                                                   "2'-O-ribose methyltransferase" = "2'-O-ribose methyltransferase",
                                                                                   "Surface glycoprotein" ="surface",
                                                                                   "ORF3"="ORF3",
                                                                                   "Envelope"="envelope",
                                                                                   "Membrane glycoprotein" = "membrane",
                                                                                   "ORF6" ="ORF6",
                                                                                   "ORF7a"="ORF7a","ORF7b"="ORF7b","ORF8"="ORF8",
                                                                                   "Nucleocapsid"="capsid"  , "ORF10"= "ORF10"
                                                                                ) 
                                                  )),
                                                  column(width = 2, actionButton("do2", h4(strong("Visualize")),  icon("eye"), width = "90px"))
                                              )),
                                     fluidRow(  downloadButton(outputId = "conserverdprimerTabledownload", label = "Download Table"),
                                                actionButton("iu12", "IUPAC Codes"),
                                                bsModal("modalExample12", "IUPAC Codes", "iu12", size = "small",includeMarkdown("mdfiles/iupac.md")),
                                       actionButton("explanationbutton1", "Click for Description", icon("info")),
                                               bsModal("modal1", "", "explanationbutton1", size = "large",includeMarkdown("mdfiles/explanation1.md")),
                                       
                                             actionButton("tabBut3", "View Lineage Mutations"), 
                                             bsModal("modallineage3", "Data Table", "tabBut3", size = "large",
                                                     dataTableOutput("lineagetable3")),
                                               DT::dataTableOutput("conserverdprimerTable")),br(),
                                   
                                     br(),
                                     fluidRow(box(title=("Detailed Variation Information"),width = 12,collapsible = TRUE,solidHeader = TRUE,
                                                  DT::dataTableOutput("mutationinformation_cons"), br(),
                                                  downloadButton(outputId = "down_mut_cons", label = "Download"))),
                                     fluidRow( includeMarkdown("mdfiles/plotexplain.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"
                                               ),
                                     
                                     tags$style(type = 'text/css',".modal-lg{ max-width: 1500px !important; width: 1500px !important; }") ,
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                                          plotOutput("conservedplot",
                                                                     dblclick = "condblclick",brush = brushOpts(id = "conbrush",resetOnNew = TRUE)))),
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                                          plotOutput("conservedphylop",
                                                                     dblclick = "condblclick",brush = brushOpts(id = "conbrush",resetOnNew = TRUE)))),
                                     fluidRow(actionButton("msa2", "Click for MSA result (works on only Safari)"),
                                              msaROutput("conservedMSA", width = "100%", height = "100%") ) ),
                            # tabset panel
                            #####
                            tabPanel(h4("Degenerate Primer pairs"),
                                     fluidRow(
                                       box(
                                         title = h4("Select the primer pair to display mutation and PhyloP plots for the amplicon and alignment to the reference genome"),
                                         status = "danger",collapsible = TRUE,solidHeader = TRUE,width = 12,
                                         fluidRow(strong("*Degenerate primer pairs that exhibit mutations found in different SARS-CoV-2 Lineages")),
                                         column(width = 6, selectInput("degenerateRegion",h5("Select the region that you want to amplify"),
                                                                       c( "All"="all",
                                                                         "nsp1" = "nsp1",
                                                                          "nsp2" = "nsp2",
                                                                          "nsp3" = "nsp3",
                                                                          "nsp4" = "nsp4",
                                                                          "3C-like proteinase" = "3C-like proteinase",
                                                                          "nsp6" = "nsp6","nsp7" ="nsp7", "nsp8" = "nsp8",
                                                                          "nsp9" = "nsp9",
                                                                           "nsp10" = "nsp10",
                                                                          "RNA-dependent RNA polymerase" = "RNA-dependent RNA polymerase",
                                                                          "Helicase" = "helicase",
                                                                           "3'-to-5' exonuclease" = "3'-to-5' exonuclease",
                                                                          "endoRNAse" ="endoRNAse", "2'-O-ribose methyltransferase" = "2'-O-ribose methyltransferase",
                                                                          "Surface glycoprotein" ="surface",
                                                                          "ORF3"="ORF3",
                                                                          "Envelope"="envelope",
                                                                          "Membrane glycoprotein" = "membrane","ORF6" ="ORF6",
                                                                          "ORF7a"="ORF7a","ORF7b"="ORF7b",
                                                                         "ORF8"="ORF8",
                                                                          "Nucleocapsid"="capsid" , "ORF10"= "ORF10" 
                                                                       )
                                         )),
                                         
                                         
                                         column(width = 2, actionButton("do4", h4(strong("Visualize")),  icon("eye"), width = "90px"))
                                       ) #box
                                     ), #fluidrow
                                     
                                
                                     selectInput("lineage", h4("Select lineage to display the primers"),
                                                 c("All"="all", "B.1.1.529 (omicron)"="B.1.529","P.1 (gamma)"="P.1","P.2"="P.2","B.1.1.7  (alpha)"="B.1.1.7","B.1.351 (beta)"="B.1.351",
                                                   "B.1.617.1"="B.1.617.1","B.1.617.2 (Delta)"="B.1.617.2","C.37 (Lambda)"="C.37","B.1.621 (Mu)"="B.1.621","C.1.2"="C.1.2")),
                                     fluidRow( downloadButton(outputId = "degenerateprimerTabledownload", label = "Download Table"),actionButton("iu", "IUPAC Codes"),
                                               bsModal("modalExample", "IUPAC Codes", "iu", size = "small",includeMarkdown("mdfiles/iupac.md")),
                                               actionButton("explanationbutton", "Click for Description", icon("info")),
                                               bsModal("modal2", "", "explanationbutton", size = "large",includeMarkdown("mdfiles/explanation2.md")),
                                               actionButton("tabBut", "View Lineage Mutations"), 
                                               bsModal("modallineage", "Data Table", "tabBut", size = "large",
                                                       dataTableOutput("lineagetable")), 
                                               DT::dataTableOutput("degenerateprimerTable")
                                               
                                     ),
                                     
                                     br(),
                                     fluidRow(box(title=("Detailed Variation Information"),width = 12,collapsible = TRUE,solidHeader = TRUE,
                                                  DT::dataTableOutput("mutationinformation_deg"), 
                                                  br(),downloadButton(outputId = "down_mut_deg", label = "Download"))),
                                     fluidRow( includeMarkdown("mdfiles/plotexplain.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"),
                                     
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                       plotOutput("degenerateplot",dblclick = "degdblclick",brush = brushOpts(id = "degbrush",resetOnNew = TRUE)))),
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                                          plotOutput("degeneratephylop",dblclick = "degdblclick",brush = brushOpts(id = "degbrush",resetOnNew = TRUE)))),
                                     fluidRow(box(actionButton("msa3", "Click for MSA result (works on only Safari)"),
                                                  msaROutput("degenerateMSA", width = "100%", height = "100%") )  )
                                     
                            ), #degenerate tab panel
                            #####
                            tabPanel(h4("Conserved Primer pairs for highly variable regions"),
                                     fluidRow(
                                       box(
                                         title = h4("Select the primer pair to display mutation and PhyloP plots for the amplicon and alignment to the reference genome"),
                                         status = "danger",collapsible = TRUE,solidHeader = TRUE,width = 12,
                                         fluidRow(strong("*Conserved primers flanking highly variable regions")),
                                         column(width = 6, selectInput("variableregion",h5("Select the region that you want to amplify"),
                                                                       c( "All"="all",
                                                                          "nsp1" = "nsp1",
                                                                          "nsp2" = "nsp2",
                                                                          "nsp3" = "nsp3",
                                                                          "nsp4" = "nsp4",
                                                                          "3C-like proteinase" = "3C-like proteinase",
                                                                           "nsp6" = "nsp6","nsp7" ="nsp7", "nsp8" = "nsp8",
                                                                          "nsp9" = "nsp9",
                                                                           "nsp10" = "nsp10",
                                                                          "RNA-dependent RNA polymerase" = "RNA-dependent RNA polymerase",
                                                                          "Helicase" = "helicase",
                                                                          "3'-to-5' exonuclease" = "3'-to-5' exonuclease",
                                                                          "endoRNAse" ="endoRNAse", "2'-O-ribose methyltransferase" = "2'-O-ribose methyltransferase",
                                                                          "Surface glycoprotein" ="surface",
                                                                          "ORF3"="ORF3",
                                                                          "Envelope"="envelope",
                                                                           "Membrane glycoprotein" = "membrane",  "ORF6" ="ORF6",
                                                                         "ORF7a"="ORF7a","ORF7b"="ORF7b",
                                                                             "ORF8"="ORF8",
                                                                          "Nucleocapsid"="capsid"    , "ORF10"= "ORF10" 
                                                                       )
                                         )),
                                         
                                         
                                         column(width = 2, actionButton("do5", h4(strong("Visualize")),  icon("eye"), width = "90px"))
                                       ) #box
                                     ), #fluidrow
                                     
                              
                                     fluidRow( downloadButton(outputId = "variableprimerTabledownload", label = "Download Table"),
                                               actionButton("iu2", "IUPAC Codes"),
                                               bsModal("modalex", "IUPAC Codes", "iu2", size = "small",includeMarkdown("mdfiles/iupac.md")),
                                               actionButton("explanationbutton3", "Click for Description", icon("info")),
                                               bsModal("modal3", "", "explanationbutton3", size = "large",includeMarkdown("mdfiles/explanation2.md")),
                                               br(), DT::dataTableOutput("variableprimerTable")
                                               
                                     ),
                                     
                                     br(),
                                     fluidRow(box(title=("Detailed Variation Information"),width = 12,collapsible = TRUE,solidHeader = TRUE,
                                                  DT::dataTableOutput("mutationinformation_var"),
                                                  br(),
                                                  downloadButton(outputId = "down_mut_var", label = "Download"))),
                                     fluidRow( includeMarkdown("mdfiles/plotexplain.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"),
                                     
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                                          plotOutput("variableteplot",
                                                                     dblclick = "variabledblclick",brush = brushOpts(id = "variablebrush",resetOnNew = TRUE)))),
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                                          plotOutput("variablephylop",
                                                                     dblclick = "variabledblclick",brush = brushOpts(id = "variablebrush",resetOnNew = TRUE)))),
                                     fluidRow(box(actionButton("msaV", "Click for MSA result (works on only Safari)"),
                                                  msaROutput("variableMSA", width = "100%", height = "100%") )  )
                                     
                            ) #highly variable tab panel
                            #####
                          )# tabset panel
                  ),# tabItem
                  # Search primer
                  tabItem(tabName = "get_tab",
                          fluidRow(
                            box(
                              title = "Search your primers in reference genome",
                              status = "danger",
                              solidHeader = TRUE,
                              width = 12,
                              collapsible = TRUE,
                              fluidRow(
                              column(
                                  width = 3,
                                  textInput(
                                    "forwardtext",
                                    label = h5("Your Forward Primer (5'->3')"),
                                    value = "GCTAACAATGCTGCAATCGTGCT"
                                  ),
                                  fluidRow(includeMarkdown("mdfiles/words.md"))
                                ),
                                column(
                                  width = 3,
                                  offset = 1,
                                  textInput(
                                    "reversetext",
                                    label = h5("Your Reverse Primer (5'->3')"),
                                    value = "CCTACTGCTGCCTGGAGTTGAA")
                                ),
                                column(
                                  width = 3,
                                  offset = 1,
                                  fluidRow( textInput(
                                    "probetext",
                                    label = h5("Your Probe (5'->3')"),
                                    value = "CAGTCAAGCCTCTTCTCGTTCCT")),
                                  fluidRow(radioButtons("pr","Include probe",c("Yes","No")))
                                 
                                ),
                                column(2, actionButton("do3", h4(strong("Search")), icon("search-plus")))),
                              fluidRow(box(width = 12, DT::dataTableOutput("searchyourprimersTable")))
                            )
                          ),
                          # fluid row ends
                          fluidRow(
                            box(
                              width = 12,column(2,actionButton("tabBut2", "View Lineage Mutations"), 
                                                bsModal("modallineage2", "", "tabBut2", size = "large",
                                                        dataTableOutput("lineagetable2"))),
                             column(4,offset=3,
                              actionButton("align", h4(strong("Align your primers")), icon("align-justify"))),
                              fluidRow(box(title="Detailed Variation Information",width = 12,collapsible = TRUE,solidHeader = T,
                                           column(4,conditionalPanel(condition = "input.pr=='Yes'",
                                                                     wellPanel(radioButtons("mutinfo","Nucleotide positions;",
                                                                                            c("From Forward to Reverse"="f2r","From Probe to Reverse"="p2r"), inline = T)))),
                                           column(2, offset = 6,downloadButton(outputId = "down_mut_search", label = "Download")),
                                      DT::dataTableOutput("mutationinformation_search")
                             
                                         )),
                            
                              fluidRow( includeMarkdown("mdfiles/plotexplainsearch.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"),
                              
                              fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                   color = getOption("spinner.color", default = "#0275D8"),
                                                   plotOutput("searchplot", click = "plot_click", hover = "plot_hover",
                                                              dblclick = "sdblclick",brush = brushOpts(id = "sbrush",resetOnNew = TRUE)))),
                                       #fluidRow(bsModal("modalExample", "Your plot", "go", size = "large",dataTableOutput ("interact"),downloadButton('downloadmut', 'Download')),
                                       #verbatimTextOutput("plot_clickinfo"),
                              fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                   color = getOption("spinner.color", default = "#0275D8"),
                                                   plotOutput("searchphylop" ,dblclick = "sdblclick",brush = brushOpts(id = "sbrush",resetOnNew = TRUE)))),
                              actionButton("msa4", "Click for MSA result (works on only Safari)"),
                              msaROutput("searchMSA", width = "100%", height = "100%")
                            )
                          ))
                  
                  
                  
                )#tabitems
                
  )#dahsboard body
)# fluidpage
######