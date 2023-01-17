
library(shiny)
library(BiocManager)
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
library(shinydashboardPlus)
library(shinycssloaders)
source("calculate_GC.R")
source("calculate_tm.R")
source("find_Tm.R") 
source("find_Gc.R") 
source("forw_degenerate.R")
source("reverse_degenerate.R") 
source("specify_decimal.R")

ui <- dashboardPage(
  skin='red',
  #-----
  dashboardHeader(titleWidth = 300,title = "CoVrimer",
                  dropdownMenu(type = "notifications",
                               notificationItem( text = "Last update of data: 2021-05-12",icon = icon("calendar"),status = "warning")),
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
  
  dashboardBody(tags$head(tags$style(HTML('
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
                  
                  tabItem(tabName = "introduction",  img(src = "cov.svg", align = "center"),
                          p(strong(h2("How to use the CoVrimer app:"))),
                          includeMarkdown("mdfiles/intro.md"),
                          
                          actionButton("go", h4("Click here for more details")),
                          uiOutput("details"),
                  ),
                  tabItem(tabName = "publishedprimers",
                          fluidRow(column(8, align="center", offset = 2,h3(strong("Table of published primer sets")), 
                                          h4("Please select a primer pair to display mutation and conservation plots for the amplicon and alignment to the reference genome")),
                                   style="text-align:justify;color:white;background-color:#535557;padding:1px;border-radius:10px" ),
                          box(solidHeader = FALSE,width=12,collapsible = TRUE,
                              style = 'display:block;width:100%;overflow-x: scroll;',
                              fluidRow(actionButton("description1", "Click for Description",icon("info")),
                                       bsModal("modal", "", "description1", size = "large",includeMarkdown("mdfiles/descriptionforpublishedprimers.md"))
                                       ,DT::dataTableOutput("publishedprimerTable")) ),
                       
                          fluidRow(box(title="Detailed Variation Information",width = 12,collapsible = TRUE,solidHeader = TRUE,
                                       DT::dataTableOutput("click_info"))),
                          
                          fluidRow( includeMarkdown("mdfiles/plotexplain.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"),
                          
                          fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                               color = getOption("spinner.color", default = "#0275D8"),
                                               plotOutput("publishedprimerplot",dblclick = "A1dblclick",
                                                          brush = brushOpts(id = "A1brush",resetOnNew = TRUE)))
                                   
                             ),
                          
                          fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                               color = getOption("spinner.color", default = "#0275D8"), plotOutput("publishedphylop",dblclick = "A1dblclick",
                                               brush = brushOpts(id = "A1brush",resetOnNew = TRUE))))
                          ,
                          fluidRow(actionButton("msa1", "Click for MSA result (works on only Safari)"),
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
                                                                                   #"nsp7" ="nsp7",
                                                                                   "nsp8" = "nsp8",
                                                                                   "nsp9" = "nsp9",
                                                                                   "nsp10" = "nsp10",
                                                                                   "RNA-dependent RNA polymerase" = "RNA-dependent RNA polymerase",
                                                                                   "Helicase" = "helicase",
                                                                                   "3'-to-5' exonuclease" = "3'-to-5' exonuclease",
                                                                                   "endoRNAse" ="endoRNAse",
                                                                                   #"2'-O-ribose methyltransferase" = "2'-O-ribose methyltransferase",
                                                                                   "Surface glycoprotein" ="surface",
                                                                                   "ORF3"="ORF3",
                                                                                  # "Envelope"="envelope",
                                                                                   "Membrane glycoprotein" = "membrane",
                                                                                   #"ORF6" ="ORF6",
                                                                                   "ORF7a"="ORF7a",#"ORF7b"="ORF7b","ORF8"="ORF8",
                                                                                   "Nucleocapsid"="capsid" # , "ORF10"= "ORF10"
                                                                                ) 
                                                  )),
                                                  column(width = 2, actionButton("do2", h4(strong("Visualize")),  icon("eye"), width = "90px"))
                                              )),
                                     fluidRow( actionButton("explanationbutton1", "Click for Description", icon("info")),
                                               bsModal("modal1", "", "explanationbutton1", size = "large",includeMarkdown("mdfiles/explanation1.md"))
                                               ,downloadButton(outputId = "conserverdprimerTabledownload", label = "Download Table"),DT::dataTableOutput("conserverdprimerTable")),br(),
                                   
                                     
                                     fluidRow(box(title=("Detailed Variation Information"),width = 12,collapsible = TRUE,solidHeader = TRUE,
                                                  DT::dataTableOutput("mutationinformation_cons"))),
                                     fluidRow( includeMarkdown("mdfiles/plotexplain.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"),
                                     
                                     
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                                          plotOutput("conservedplot",dblclick = "condblclick",brush = brushOpts(id = "conbrush",resetOnNew = TRUE)))),
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                                          plotOutput("conservedphylop" ,dblclick = "condblclick",brush = brushOpts(id = "conbrush",resetOnNew = TRUE)))),
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
                                                                         # "nsp1" = "nsp1",
                                                                          "nsp2" = "nsp2",
                                                                          "nsp3" = "nsp3",
                                                                          "nsp4" = "nsp4",
                                                                          "3C-like proteinase" = "3C-like proteinase",
                                                                          "nsp6" = "nsp6",#"nsp7" ="nsp7", "nsp8" = "nsp8",
                                                                          "nsp9" = "nsp9",
                                                                          # "nsp10" = "nsp10",
                                                                          #"RNA-dependent RNA polymerase" = "RNA-dependent RNA polymerase",
                                                                          "Helicase" = "helicase",
                                                                          # "3'-to-5' exonuclease" = "3'-to-5' exonuclease",
                                                                          #"endoRNAse" ="endoRNAse", "2'-O-ribose methyltransferase" = "2'-O-ribose methyltransferase",
                                                                          "Surface glycoprotein" ="surface",
                                                                          "ORF3"="ORF3",
                                                                          "Envelope"="envelope",
                                                                         # "Membrane glycoprotein" = "membrane","ORF6" ="ORF6",
                                                                          "ORF7a"="ORF7a",#"ORF7b"="ORF7b",
                                                                         "ORF8"="ORF8",
                                                                          "Nucleocapsid"="capsid" #, "ORF10"= "ORF10" 
                                                                       )
                                         )),
                                         
                                         
                                         column(width = 2, actionButton("do4", h4(strong("Visualize")),  icon("eye"), width = "90px"))
                                       ) #box
                                     ), #fluidrow
                                     
                                     br(),
                                     fluidRow( downloadButton(outputId = "degenerateprimerTabledownload", label = "Download Table"),actionButton("iu", "IUPAC Codes"),
                                               bsModal("modalExample", "IUPAC Codes", "iu", size = "small",includeMarkdown("mdfiles/iupac.md")),
                                               actionButton("explanationbutton", "Click for Description", icon("info")),
                                               bsModal("modal2", "", "explanationbutton", size = "large",includeMarkdown("mdfiles/explanation2.md")),
                                               DT::dataTableOutput("degenerateprimerTable")
                                               
                                     ),
                                     
                                     
                                     fluidRow(box(title=("Detailed Variation Information"),width = 12,collapsible = TRUE,solidHeader = TRUE,
                                                  DT::dataTableOutput("mutationinformation_deg"))),
                                     fluidRow( includeMarkdown("mdfiles/plotexplain.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"),
                                     
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                       plotOutput("degenerateplot",dblclick = "degdblclick",brush = brushOpts(id = "degbrush",resetOnNew = TRUE)))),
                                     fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                          color = getOption("spinner.color", default = "#0275D8"),
                                                          plotOutput("degeneratephylop",dblclick = "degdblclick",brush = brushOpts(id = "degbrush",resetOnNew = TRUE) ))),
                                     fluidRow(box(actionButton("msa3", "Click for MSA result (works on only Safari)"),
                                                  msaROutput("degenerateMSA", width = "100%", height = "100%") )  )
                                     
                            ) #degenerate tab panel
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
                                  width = 2,
                                  textInput(
                                    "forwardtext",
                                    label = h5("Your Forward Primer (5'->3')"),
                                    value = "ACAGGTACGTTAATAGTTAATAGCGT"
                                  )
                                ),
                                column(
                                  width = 2,
                                  offset = 1,
                                  textInput(
                                    "reversetext",
                                    label = h5("Your Reverse Primer (5'->3')"),
                                    value = "ATATTGCAGCAGTACGCACACA")
                                ),
                                column(2, actionButton("do3", h4(strong("Search")), icon("search-plus")))),
                              fluidRow(box(width = 12, DT::dataTableOutput("searchyourprimersTable")))
                            )
                          ),
                          # fluid row ends
                          fluidRow(
                            box(
                              width = 12,
                              title = "See the Alignment",
                              actionButton("align", h4(strong("Align")), icon("align-justify")),
                              fluidRow(box(h4("Detailed Variation Information"),width = 12,collapsible = TRUE,solidHeader = TRUE,
                                           DT::dataTableOutput("mutationinformation_search"))),
                              fluidRow( includeMarkdown("mdfiles/plotexplain.md"),style="text-align:justify;color:black;background-color:#DEE3E6;padding:1px;border-radius:10px"),
                              
                              fluidRow(withSpinner(type = getOption("spinner.type", default = 6),
                                                   color = getOption("spinner.color", default = "#0275D8"),
                                                   plotOutput("searchplot",dblclick = "sdblclick",brush = brushOpts(id = "sbrush",resetOnNew = TRUE)))),
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
