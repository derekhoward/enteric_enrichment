library(shiny)
library(ggplot2)
library(magrittr)
library(here)
library(readr)
library(plotly)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(shinyjs)
source("./AUCFunction.R")
source("./string_processing.R")

apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  #tags$head(includeHTML("google-analytics.html")),
  # App title ----
  titlePanel("Polygenic tester for human enteric nervous system"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Selector for choosing dataset ----
      selectInput(
        inputId = 'dataset',
        label = 'Dataset:',
        choices = c(
          'Enteric neurons and glia',
          'Enteric cells'
        )
      ),
      textAreaInput(
        inputId = "genelist",
        label = "Input your gene list:",
        value = 'ITPR3\nSPRY2\nUBE2M\nCKB\nHSPA1A\nSAMD5\nMANEA\nZC3H7B\nTMEM106B\nASXL3\nEXT1\nRERE\nPPP6C',
        rows = 15
      ),
      selectInput(
        inputId = 'species',
        label = 'Species of input genes:',
        choices = c('Human', 'Mouse', 'Rhesus Macaque')
      ),
      actionButton(inputId = "submit",
                   label = "Submit"),
      br(),
      br(),
      downloadButton(outputId = "download_data", label = "Download results as .csv"),
      hr(),
      tags$b("Data made available by the Broad institute's Single Cell Portal:"),
      br(),
      tags$a(href="https://singlecell.broadinstitute.org/single_cell/study/SCP1038/the-human-and-mouse-enteric-nervous-system-at-single-cell-resolution", "The human and mouse enteric nervous system at single cell resolution"),
      br()
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(
        id = "main",
        # Output: Verbatim text for data summary ----
        verbatimTextOutput("summary"),
        br(),
        plotlyOutput("dotplot"),
        br(),
        dataTableOutput("view"),
        verbatimTextOutput("info")
      )
    )
  )
)


# Define server logic process and output top cortical layers/zones ----
server <- function(input, output) {
  output$summary <- renderPrint({
    cat("\nResults will load here when complete")
    cat("\n")
    #print(gc())
    #print(Sys.info()['nodename'])
  })
  
  observeEvent(input$submit, {
    start <- Sys.time()
    
    cleaned_gene_list <-
      isolate(process_input_genes(input$genelist))
    
    # load reference data
    if (input$dataset == 'Enteric neurons and glia') {
      ranks_matrix <- read_csv(here('data/processed/ent_neurons_glia_ranks.csv'))
      unique_genes <- ranks_matrix$gene_symbol
      paste('loaded data')
      
      # could be useful to load in raw expression data instead of just ranks matrix in order to plot expression values
      
      #load('./data/processed/developing_cortical_zones_expression_matrix.Rdata', verbose = TRUE)
      enteric_expression <- read_csv(here('data/processed/ent_neurons_glia_logExp.csv'))
      tidy_expression <- pivot_longer(enteric_expression, -gene_symbol, names_to='cell_id', values_to = 'logExp')
      #tidy_expression <- cortical_zones_expression_matrix %>%
      #  gather(key = zones, value = expression,-gene_symbol)
      
      #target_gene_symbols <- 'Human'
      cell_types <- tribble(
        ~ cell_id, ~ cell_type,
        'PEMN_1', 'putative excitatory motor neurons 1',
        'PEMN_2', 'putative excitatory motor neurons 2',
        'PEMN_3', 'putative excitatory motor neurons 3',
        'PEMN_4', 'putative excitatory motor neurons 4',
        'PIMN_1', 'putative inhibitory motor neurons 1',
        'PIMN_2', 'putative inhibitory motor neurons 2',
        'PIMN_3', 'putative inhibitory motor neurons 3',
        'PIMN_4', 'putative inhibitory motor neurons 4',
        'PIMN_5', 'putative inhibitory motor neurons 5',
        'PIN_1', 'putative inhibitory neurons 1',
        'PIN_2', 'putative inhibitory neurons 2',
        'PSN', 'putative sensory neurons',
        'PSVN_1', 'putative secretomotor/vasodilator neurons 1',
        'PSVN_2', 'putative secretomotor/vasodilator neurons 2',
        'Glia_1', 'glia 1',
        'Glia_2', 'glia 2',
        'Glia_3', 'glia 3',
        'Glia_4', 'glia 4',
        'Glia_5', 'glia 5',
        'Glia_6', 'glia 6'
      )
      
      #cleaned_gene_list <- convert2human(input_genes = cleaned_gene_list, in_species = input$species)
      
    } else {
      ranks_matrix <- read_csv(here('./data/processed/all_ent_ranks.csv'))
      unique_genes <- ranks_matrix$gene_symbol
      paste('loaded data')
      
      enteric_expression <- read_csv(here('data/processed/all_ent_logExp.csv'))
      tidy_expression <- pivot_longer(enteric_expression, -gene_symbol, names_to='cell_id', values_to = 'logExp')
      
      cell_types <- tribble(
        ~ cell_id, ~ cell_type,
        'Adipose', 'Adipose',
        'Epithelial', 'Epithelial',
        'Fibroblast_1', 'Fibroblast 1',
        'Fibroblast_2', 'Fibroblast 2',
        'Glia_1', 'Glia 1',
        'Glia_2', 'Glia 2',
        'Glia_3', 'Glia 3',
        'H3', 'Patient specific - H3',
        'ICCs', 'ICCs',
        'Lymphatic', 'Lymphatic',
        'Macrophage', 'Macrophage',
        'Mesothelial', 'Mesothelial',
        'MHC.I_H1', 'Patient specific - MHC.I H1',
        'MHC.I_H9', 'Patient specific - MHC.I H9',
        'Myocyte_1', 'Myocyte 1',
        'Myocyte_2', 'Myocyte 2',
        'Myocyte_3', 'Myocyte 3',
        'Myocyte_4', 'Myocyte 4',
        'Myocyte_5', 'Myocyte 5',
        'Neuron', 'Neuron',
        'OXPHOS_H3', 'Patient specific - OXPHOS H3',
        'Pericytes', 'Pericytes',
        'T_cells', 'T cells',
        'Vascular_1', 'Vascular 1',
        'Vascular_2', 'Vascular 2'
      )
    }
    print(paste0("Before time taken:", Sys.time() - start))
    
    cleaned_gene_list <- convert2human(input_genes = cleaned_gene_list, in_species = input$species)
    
    #for indices - use dplyr for ease
    forIndices <- as_tibble(ranks_matrix$gene_symbol)
    names(forIndices) <- 'gene_symbol'
    forIndices %<>% mutate(isTargetGene = gene_symbol %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    # only columns from cortical zones remain in df
    df <- ranks_matrix %>%
      select(-gene_symbol)
    
    AUROC <- map_df(df, auroc_analytic, targetIndices)
    wilcox_tests <- map_df(df, apply_MWU, targetIndices)
    
    # group results together in a single table
    table <- bind_cols(gather(AUROC, key = cell_id, value = AUROC), 
                       gather(wilcox_tests, value = pValue)) %>%
      select(-key)
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    

    
    # these are the values for the results table
    table %<>% arrange(-AUROC)
    table %<>% mutate(pValue = signif(pValue, digits = 3), 
                      AUROC = signif(AUROC, digits = 3),
                      adjusted_P = signif(p.adjust(pValue), digits = 3))
    table %<>% left_join(cell_types, by='cell_id') %>%
      select(cell_id, cell_type, AUROC, pValue, adjusted_P)
    
    
    #table$order
#    if (input$dataset == 'Cortical layers from Developing Human Brain Atlas') {
#      zone_order <-c("ventricular zone", "subventricular zone", "intermediate zone", "subplate zone", "cortical plate",
#                     "marginal zone","subpial granular zone")
#      x <- tibble(laminar_order = 1:7, zone = zone_order)
#    } else {
#      #from Blueprint documentation: http://download.alleninstitute.org/nhp/Prenatal_Macaque_LMD_Microarray/neuroanatomical_guides_for_LMD_sampling/V1/
#      #The layers annotated for LMD are indicated in the right panel and include marginal zone (mz), layer 2 (2), layer 2/3 (2/3), layer 3 (3), 
#      zone_order <- c("white matter", "marginal zone", "layer I", "layer II", "layer II/III", "layer III", 
#                      #outer cortical plate (cpo), layer 4 (4), layer 4A (4A), layer 4B (4B), layer 4Ca (4Ca), layer 4Cb (4Cb), layer 5 (5), layer 6 (6), cortical plate (cp), 
#                      "cortical plate (outer)",  "layer IV", "layer V", "layer VI", "cortical plate", 
#                      #inner cortical plate (cpi), subplate (sp), intermediate zone (iz), intermediate cell dense zone (icd), transitory migratory zone (tmz), 
#                      "cortical plate (inner)","subplate","intermediate zone","transitory migratory zone", 
#                      #outer fiber (plexiform) zone (ofz), subventricular zone (sz), outer subventricular zone (szo), inner fiber (plexiform) zone (ifz), 
#                      "outer fiber zone","subventricular zone","subventricular zone (outer)", "inner fiber zone", 
                      #inner subventricular zone (szi), outer ventricular zone (vzo), inner ventricular zone (vzi), and ventricular zone (vz). 
#                      "subventricular zone (inner)","ventricular zone (outer)","ventricular zone (inner)", "ventricular zone")
#      zone_order <- rev(zone_order) #reverse so it lines up with human direction
      
#      x <- tibble(laminar_order = 1:length(zone_order), zone = zone_order)
#    }
#    table <- inner_join(table, x, by = 'zone')

    #      zone_order <-c("ventricular zone", "subventricular zone", "intermediate zone", "subplate zone", "cortical plate",
    #                     "marginal zone","subpial granular zone")
    #      x <- tibble(laminar_order = 1:7, zone = zone_order)
    

    
    selected_values <- reactive({
      req(cleaned_gene_list)
      selected_values <- tidy_expression %>% filter(gene_symbol %in% cleaned_gene_list)
      #selected_values <- selected_values %<>% left_join(cell_types, by='cell_id')
      #selected_values %<>% mutate(zones = factor(zones, levels = zone_order))
    })
    
    # these are the values used for the plot
    #selected_values <- tidy_expression %>% filter(gene_symbol %in% cleaned_gene_list)
    #selected_values %<>% mutate(zones = factor(zones, levels=c("ventricular zone", "subventricular zone", "intermediate zone", "subplate zone", "cortical plate", "marginal zone", "subpial granular zone")))
    
    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste(
        "\nGenes found in data:",
        sum(cleaned_gene_list %in% unique_genes),
        "of",
        length(cleaned_gene_list)
      ))
    })
    
    output$view <- renderDataTable({
      table
    }, escape = FALSE)
    
    output$dotplot <- renderPlotly({
      p <- ggplot(selected_values(), aes(x = cell_id, y = logExp, names=gene_symbol))
      #figure depends on size
      if (length(cleaned_gene_list) > 20) {
        p <-  p + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.25, width = 0.1)
      } else {
        p <- p + geom_jitter(width = 0.1)
      }
      p <- p + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 0))
      p <- p + geom_hline(yintercept = 0, color='darkgrey', size=0.4)
      ggplotly(p) #%>% layout(dragmode = "select")
    })
    
    output$download_data <-
      downloadHandler(
        filename = "polygenic_enteric_AUC_results.csv",
        content = function(file) {
          write_csv(table, file)
        }
      )
    
  })
}

shinyApp(ui, server)