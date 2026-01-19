# =============================================================================
# RNA-seq Results Explorer - Multi-Species Dashboard
# =============================================================================
# Interactive Shiny dashboard for exploring DE results across species
# Compares weep mutant vs wild-type in Arabidopsis, Pennycress, and Soybean

library(shiny)
library(shinydashboard)
library(shinycssloaders)  # BUG_007 FIX: For loading spinners
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(ggVennDiagram)  # For Venn diagram visualization

# =============================================================================
# Data Paths
# =============================================================================

base_path <- "/Users/benficaa/Desktop/RNAseq/RNASeq_Results_Explorer/Data"

# =============================================================================
# Load All Data
# =============================================================================

# HPI/Overall DE data (WEEP vs WT) - All genes with DE results
# BUG_005 FIX: Use *_DE_all_genes.csv for all species (includes all genes, not just DEGs)
arab_hpi <- tryCatch(read.csv(file.path(base_path, "Arabidopsis/Arabidopsis_DE_all_genes.csv")), error = function(e) NULL)
penny_hpi <- tryCatch(read.csv(file.path(base_path, "Pennycress/Pennycress_DE_all_genes.csv")), error = function(e) NULL)
soy_hpi <- tryCatch(read.csv(file.path(base_path, "Soybean/Soybean_DE_all_genes.csv")), error = function(e) NULL)

# Segment-specific data
arab_seg1 <- tryCatch(read.csv(file.path(base_path, "Arabidopsis/Arabidopsis_Segment1_WT_vs_weep.csv")), error = function(e) NULL)
arab_seg2 <- tryCatch(read.csv(file.path(base_path, "Arabidopsis/Arabidopsis_Segment2_WT_vs_weep.csv")), error = function(e) NULL)
arab_seg3 <- tryCatch(read.csv(file.path(base_path, "Arabidopsis/Arabidopsis_Segment3_WT_vs_weep.csv")), error = function(e) NULL)

penny_seg1 <- tryCatch(read.csv(file.path(base_path, "Pennycress/Pennycress_Segment1_WT_vs_weep.csv")), error = function(e) NULL)
penny_seg2 <- tryCatch(read.csv(file.path(base_path, "Pennycress/Pennycress_Segment2_WT_vs_weep.csv")), error = function(e) NULL)
penny_seg3 <- tryCatch(read.csv(file.path(base_path, "Pennycress/Pennycress_Segment3_WT_vs_weep.csv")), error = function(e) NULL)

soy_seg1 <- tryCatch(read.csv(file.path(base_path, "Soybean/Soybean_Segment1_WT_vs_weep.csv")), error = function(e) NULL)
soy_seg2 <- tryCatch(read.csv(file.path(base_path, "Soybean/Soybean_Segment2_WT_vs_weep.csv")), error = function(e) NULL)
soy_seg3 <- tryCatch(read.csv(file.path(base_path, "Soybean/Soybean_Segment3_WT_vs_weep.csv")), error = function(e) NULL)

# FEATURE_003: Load combined cross-species DEG data
all_species_deg <- tryCatch(read.csv(file.path(base_path, "all_species_deg.csv")), error = function(e) NULL)

# Species colors
species_colors <- c("Arabidopsis" = "#27ae60", "Pennycress" = "#9b59b6", "Soybean" = "#f39c12")
deg_colors <- c("UP in weep" = "#e74c3c", "DOWN in weep" = "#3498db", "NOT Significant" = "#95a5a6")

# BUG_009 FIX: Segment tissue type descriptions
segment_labels <- c(
  "seg1" = "Seg1 (young laterals)",
  "seg2" = "Seg2 (no LR/primordia)",
  "seg3" = "Seg3 (root tip)"
)

segment_labels_full <- c(
  "seg1" = "Segment 1 - Primary root with young laterals",
  "seg2" = "Segment 2 - Primary root (no visible LR/primordia)",
  "seg3" = "Segment 3 - Primary root tip"
)

# Helper function to get segment label
get_seg_label <- function(seg_code, full = FALSE) {
  if (full) {
    return(segment_labels_full[seg_code])
  }
  return(segment_labels[seg_code])
}

# =============================================================================
# UI
# =============================================================================

ui <- dashboardPage(
  skin = "green",
  
  dashboardHeader(title = "RNA-seq DE Explorer"),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Explore WEEP vs WT", tabName = "wt_vs_weep", icon = icon("dna"),
               badgeLabel = "HPI", badgeColor = "green"),
      menuItem("Explore Within Segments", tabName = "within_segments", icon = icon("layer-group")),
      menuItem("Compare Between Segments", tabName = "between_segments", icon = icon("exchange-alt")),
      menuItem("Explore Between Species", tabName = "between_species", icon = icon("globe"),
               badgeLabel = "OG", badgeColor = "purple")
    ),
    hr(),
    h5("  Quick Stats", style = "padding-left: 15px; color: #888;"),
    verbatimTextOutput("quick_stats", placeholder = TRUE)
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #ecf0f5; }
        .box { box-shadow: 0 1px 3px rgba(0,0,0,0.12); }
        .info-box { min-height: 80px; }
        .nav-tabs-custom > .tab-content { padding: 10px; }
        .species-arabidopsis { color: #27ae60; font-weight: bold; }
        .species-pennycress { color: #9b59b6; font-weight: bold; }
        .species-soybean { color: #f39c12; font-weight: bold; }
      "))
    ),
    
    tabItems(
      # =========================================================================
      # Tab 1: Explore WEEP vs WT (HPI genes)
      # =========================================================================
      tabItem(
        tabName = "wt_vs_weep",
        fluidRow(
          box(
            title = "Settings", status = "success", solidHeader = TRUE, width = 3,
            selectInput("hpi_species", "Select Species:",
                        choices = c("Arabidopsis", "Pennycress", "Soybean"),
                        selected = "Arabidopsis"),
            hr(),
            h5("Significance Filters"),
            sliderInput("hpi_padj", "Max adjusted p-value:",
                        min = 0.001, max = 0.1, value = 0.05, step = 0.001),
            # BUG_008/BUG_010 FIX: Add text boxes alongside slider for log2FC range
            # Extended range to -10 to 10 to accommodate all data values
            h5("log2 Fold Change Range:"),
            fluidRow(
              column(6, numericInput("hpi_lfc_min", "Min:", value = -10, min = -15, max = 15, step = 0.25)),
              column(6, numericInput("hpi_lfc_max", "Max:", value = 10, min = -15, max = 15, step = 0.25))
            ),
            sliderInput("hpi_lfc", NULL,
                        min = -10, max = 10, value = c(-10, 10), step = 0.25),
            # BUG_008 FIX: Show UP/DOWN gene counts
            wellPanel(
              style = "padding: 8px; margin-bottom: 10px; background-color: #f9f9f9;",
              uiOutput("hpi_gene_counts")
            ),
            hr(),
            # BUG_013 FIX: Default to "All" to show all genes
            selectInput("hpi_status", "DE Status:",
                        choices = c("All", "UP in weep", "DOWN in weep", "NOT Significant"),
                        selected = "All"),
            hr(),
            # HPI Filter checkbox with tooltip
            tags$div(
              style = "display: flex; align-items: center;",
              checkboxInput("hpi_filter", "Show HPI Genes Only", value = FALSE),
              tags$span(
                title = "This will show only the HPI gene candidates that are differentially expressed",
                style = "margin-left: 5px; cursor: help; color: #3498db;",
                icon("info-circle")
              )
            ),
            hr(),
            downloadButton("download_hpi", "Download Filtered Data")
          ),
          box(
            title = "Volcano Plot: weep vs WT", status = "primary", solidHeader = TRUE, width = 9,
            # BUG_007 FIX: Add loading spinner
            withSpinner(plotlyOutput("hpi_volcano", height = "450px"), type = 6, color = "#27ae60")
          )
        ),
        fluidRow(
          box(
            title = "DE Genes Table", status = "info", solidHeader = TRUE, width = 12,
            # BUG_007 FIX: Add loading spinner
            withSpinner(DT::dataTableOutput("hpi_table"), type = 6, color = "#27ae60"),
            hr(),
            verbatimTextOutput("hpi_summary")
          )
        )
      ),
      
      # =========================================================================
      # Tab 2: Explore DE Within Segments
      # =========================================================================
      tabItem(
        tabName = "within_segments",
        fluidRow(
          box(
            title = "Settings", status = "warning", solidHeader = TRUE, width = 3,
            selectInput("seg_species", "Select Species:",
                        choices = c("Arabidopsis", "Pennycress", "Soybean"),
                        selected = "Arabidopsis"),
            selectInput("seg_segment", "Select Segment:",
                        # BUG_009 FIX: Add tissue type descriptions
                        choices = c("Segment 1 - Primary root with young laterals" = "seg1",
                                    "Segment 2 - Primary root (no visible LR/primordia)" = "seg2",
                                    "Segment 3 - Primary root tip" = "seg3"),
                        selected = "seg1"),
            hr(),
            h5("Significance Filters"),
            sliderInput("seg_padj", "Max adjusted p-value:",
                        min = 0.001, max = 0.1, value = 0.05, step = 0.001),
            # BUG_008/BUG_010 FIX: Add text boxes alongside slider for log2FC range
            # Extended range to -10 to 10 to accommodate all data values
            h5("log2 Fold Change Range:"),
            fluidRow(
              column(6, numericInput("seg_lfc_min", "Min:", value = -10, min = -15, max = 15, step = 0.25)),
              column(6, numericInput("seg_lfc_max", "Max:", value = 10, min = -15, max = 15, step = 0.25))
            ),
            sliderInput("seg_lfc", NULL,
                        min = -10, max = 10, value = c(-10, 10), step = 0.25),
            # BUG_008 FIX: Show UP/DOWN gene counts
            wellPanel(
              style = "padding: 8px; margin-bottom: 10px; background-color: #f9f9f9;",
              uiOutput("seg_gene_counts")
            ),
            hr(),
            # BUG_013 FIX: Default to "All" to show all genes
            selectInput("seg_status", "DE Status:",
                        choices = c("All", "UP in weep", "DOWN in weep", "NOT Significant"),
                        selected = "All"),
            hr(),
            # HPI Filter checkbox with tooltip
            tags$div(
              style = "display: flex; align-items: center;",
              checkboxInput("seg_hpi_filter", "Show HPI Genes Only", value = FALSE),
              tags$span(
                title = "This will show only the HPI gene candidates that are differentially expressed",
                style = "margin-left: 5px; cursor: help; color: #3498db;",
                icon("info-circle")
              )
            ),
            hr(),
            downloadButton("download_seg", "Download Filtered Data")
          ),
          box(
            title = "Volcano Plot: Segment-Specific weep vs WT", status = "primary", solidHeader = TRUE, width = 9,
            # BUG_007 FIX: Add loading spinner
            withSpinner(plotlyOutput("seg_volcano", height = "450px"), type = 6, color = "#f39c12")
          )
        ),
        fluidRow(
          box(
            title = "Segment DE Genes Table", status = "info", solidHeader = TRUE, width = 12,
            # BUG_007 FIX: Add loading spinner
            withSpinner(DT::dataTableOutput("seg_table"), type = 6, color = "#f39c12"),
            hr(),
            verbatimTextOutput("seg_summary")
          )
        )
      ),
      
      # =========================================================================
      # Tab 3: Compare Between Segments
      # =========================================================================
      tabItem(
        tabName = "between_segments",
        fluidRow(
          box(
            title = "Settings", status = "danger", solidHeader = TRUE, width = 3,
            selectInput("comp_species", "Select Species:",
                        choices = c("Arabidopsis", "Pennycress", "Soybean"),
                        selected = "Arabidopsis"),
            selectInput("comp_segments", "Compare Segments:",
                        # BUG_009 FIX: Add tissue type descriptions
                        choices = c("Seg1 (young laterals) vs Seg2 (no LR)" = "seg1_vs_seg2",
                                    "Seg1 (young laterals) vs Seg3 (root tip)" = "seg1_vs_seg3",
                                    "Seg2 (no LR) vs Seg3 (root tip)" = "seg2_vs_seg3"),
                        selected = "seg1_vs_seg2"),
            hr(),
            h5("Significance Filters"),
            sliderInput("comp_padj", "Max adjusted p-value:",
                        min = 0.001, max = 0.1, value = 0.05, step = 0.001),
            # BUG_010 FIX: Extended range to -10 to 10
            sliderInput("comp_lfc", "log2 Fold Change Range:",
                        min = -10, max = 10, value = c(-10, 10), step = 0.25),
            hr(),
            checkboxInput("comp_show_shared", "Show shared DE genes only", value = FALSE),
            # HPI Filter checkbox with tooltip
            tags$div(
              style = "display: flex; align-items: center;",
              checkboxInput("comp_hpi_filter", "Show HPI Genes Only", value = FALSE),
              tags$span(
                title = "This will show only the HPI gene candidates that are differentially expressed",
                style = "margin-left: 5px; cursor: help; color: #3498db;",
                icon("info-circle")
              )
            ),
            hr(),
            h5("Downloads"),
            downloadButton("download_comp", "Download All DEGs"),
            br(), br(),
            downloadButton("download_common", "Download Common Genes")
          ),
          box(
            title = "Segment Comparison Overview", status = "primary", solidHeader = TRUE, width = 9,
            fluidRow(
              column(6, plotlyOutput("comp_bar", height = "350px")),
              column(6, plotOutput("comp_venn", height = "350px"))
            )
          )
        ),
        fluidRow(
          box(
            title = "Common Genes Between Segments", status = "success", solidHeader = TRUE, width = 12,
            collapsible = TRUE,
            DT::dataTableOutput("comp_common_table"),
            hr(),
            verbatimTextOutput("comp_common_summary")
          )
        ),
        fluidRow(
          box(
            title = "Segment A DE Genes", status = "info", solidHeader = TRUE, width = 6,
            collapsible = TRUE, collapsed = TRUE,
            DT::dataTableOutput("comp_table_a")
          ),
          box(
            title = "Segment B DE Genes", status = "warning", solidHeader = TRUE, width = 6,
            collapsible = TRUE, collapsed = TRUE,
            DT::dataTableOutput("comp_table_b")
          )
        ),
        fluidRow(
          box(
            title = "Comparison Summary", status = "primary", solidHeader = TRUE, width = 12,
            collapsible = TRUE,
            verbatimTextOutput("comp_summary")
          )
        )
      ),

      # =========================================================================
      # Tab 4: Explore Between Species (Cross-Species Comparison)
      # =========================================================================
      tabItem(
        tabName = "between_species",
        fluidRow(
          box(
            title = "Downloads", status = "info", solidHeader = TRUE, width = 4,
            downloadButton("download_three_species", "Download Conserved OGs"),
            br(), br(),
            downloadButton("download_two_species", "Download Shared OGs")
          ),
          valueBoxOutput("total_three_species_box", width = 4),
          valueBoxOutput("total_two_species_box", width = 4)
        ),
        fluidRow(
          box(
            title = "Conserved Orthogroups (Arabidopsis, Pennycress, Soybean)",
            status = "primary", solidHeader = TRUE, width = 12,
            p("DEGs with orthologs in all three species. Direction: UP/DOWN = consistent across species, MIXED = different direction."),
            DT::dataTableOutput("three_species_table")
          )
        ),
        fluidRow(
          box(
            title = "Shared Orthogroups (2 Species)",
            status = "warning", solidHeader = TRUE, width = 12,
            p("DEGs with orthologs in exactly two species."),
            DT::dataTableOutput("two_species_table")
          )
        )
      )
    )
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {

  # ===========================================================================
  # BUG_008 FIX: Sync text boxes with sliders for Tab 1 (HPI)
  # ===========================================================================

  # When slider changes, update text boxes
  observeEvent(input$hpi_lfc, {
    updateNumericInput(session, "hpi_lfc_min", value = input$hpi_lfc[1])
    updateNumericInput(session, "hpi_lfc_max", value = input$hpi_lfc[2])
  }, ignoreInit = TRUE)

  # When min text box changes, update slider
  observeEvent(input$hpi_lfc_min, {
    req(input$hpi_lfc_min, input$hpi_lfc_max)
    if (input$hpi_lfc_min != input$hpi_lfc[1]) {
      updateSliderInput(session, "hpi_lfc", value = c(input$hpi_lfc_min, input$hpi_lfc[2]))
    }
  }, ignoreInit = TRUE)

  # When max text box changes, update slider
  observeEvent(input$hpi_lfc_max, {
    req(input$hpi_lfc_min, input$hpi_lfc_max)
    if (input$hpi_lfc_max != input$hpi_lfc[2]) {
      updateSliderInput(session, "hpi_lfc", value = c(input$hpi_lfc[1], input$hpi_lfc_max))
    }
  }, ignoreInit = TRUE)

  # ===========================================================================
  # BUG_008 FIX: Sync text boxes with sliders for Tab 2 (Segment)
  # ===========================================================================

  # When slider changes, update text boxes
  observeEvent(input$seg_lfc, {
    updateNumericInput(session, "seg_lfc_min", value = input$seg_lfc[1])
    updateNumericInput(session, "seg_lfc_max", value = input$seg_lfc[2])
  }, ignoreInit = TRUE)

  # When min text box changes, update slider
  observeEvent(input$seg_lfc_min, {
    req(input$seg_lfc_min, input$seg_lfc_max)
    if (input$seg_lfc_min != input$seg_lfc[1]) {
      updateSliderInput(session, "seg_lfc", value = c(input$seg_lfc_min, input$seg_lfc[2]))
    }
  }, ignoreInit = TRUE)

  # When max text box changes, update slider
  observeEvent(input$seg_lfc_max, {
    req(input$seg_lfc_min, input$seg_lfc_max)
    if (input$seg_lfc_max != input$seg_lfc[2]) {
      updateSliderInput(session, "seg_lfc", value = c(input$seg_lfc[1], input$seg_lfc_max))
    }
  }, ignoreInit = TRUE)

  # ===========================================================================
  # Reactive: Get HPI data by species
  # ===========================================================================
  hpi_data <- reactive({
    df <- switch(input$hpi_species,
                 "Arabidopsis" = arab_hpi,
                 "Pennycress" = penny_hpi,
                 "Soybean" = soy_hpi)

    if (is.null(df)) return(data.frame())

    # Filter by status
    if (input$hpi_status != "All") {
      df <- df[df$DEG_status == input$hpi_status, ]
    }

    # Filter by significance (only for significant genes)
    if (input$hpi_status != "All" || "padj" %in% colnames(df)) {
      df <- df[!is.na(df$padj), ]
    }

    # Filter by HPI genes only if checkbox is checked
    if (isTRUE(input$hpi_filter) && "HPI" %in% colnames(df)) {
      df <- df[df$HPI == "Yes", ]
    }

    df
  })
  
  # ===========================================================================
  # Reactive: Get segment data by species and segment
  # ===========================================================================
  segment_data <- reactive({
    species <- input$seg_species
    segment <- input$seg_segment

    df <- switch(paste(species, segment, sep = "_"),
                 "Arabidopsis_seg1" = arab_seg1,
                 "Arabidopsis_seg2" = arab_seg2,
                 "Arabidopsis_seg3" = arab_seg3,
                 "Pennycress_seg1" = penny_seg1,
                 "Pennycress_seg2" = penny_seg2,
                 "Pennycress_seg3" = penny_seg3,
                 "Soybean_seg1" = soy_seg1,
                 "Soybean_seg2" = soy_seg2,
                 "Soybean_seg3" = soy_seg3)

    if (is.null(df)) return(data.frame())

    # Filter by status
    if (input$seg_status != "All") {
      df <- df[df$DEG_status == input$seg_status, ]
    }

    # Filter by HPI genes only if checkbox is checked
    if (isTRUE(input$seg_hpi_filter) && "HPI" %in% colnames(df)) {
      df <- df[df$HPI == "Yes", ]
    }

    df
  })
  
  # ===========================================================================
  # Quick Stats (Hard-coded based on actual file contents)
  # ===========================================================================
  # Data verification (2025-01-18) - BUG_005 FIX:
  # - Using *_DE_all_genes.csv files for all species
  # - Arabidopsis_DE_all_genes.csv: 64,332 total genes, 794 DEGs
  # - Pennycress_DE_all_genes.csv: 67,733 total genes, 498 DEGs
  # - Soybean_DE_all_genes.csv: 117,810 total genes, 853 DEGs
  output$quick_stats <- renderText({
    paste0("=== Data Summary ===\n\n",
           "Total Genes:\n",
           "  Arabidopsis: 64,332\n",
           "  Pennycress: 67,733\n",
           "  Soybean: 117,810\n\n",
           "DEGs (padj < 0.05):\n",
           "  Arabidopsis: 794\n",
           "  Pennycress: 498\n",
           "  Soybean: 853")
  })

  # ===========================================================================
  # BUG_008 FIX: Gene count displays for Tab 1 (HPI)
  # ===========================================================================
  output$hpi_gene_counts <- renderUI({
    df <- switch(input$hpi_species,
                 "Arabidopsis" = arab_hpi,
                 "Pennycress" = penny_hpi,
                 "Soybean" = soy_hpi)

    if (is.null(df) || nrow(df) == 0) {
      return(tags$div(style = "color: #999;", "No data"))
    }

    # Filter by current log2FC range and padj
    filtered <- df[!is.na(df$padj) & df$padj < input$hpi_padj &
                     df$log2FoldChange >= input$hpi_lfc[1] &
                     df$log2FoldChange <= input$hpi_lfc[2], ]

    up_count <- sum(filtered$DEG_status == "UP in weep", na.rm = TRUE)
    down_count <- sum(filtered$DEG_status == "DOWN in weep", na.rm = TRUE)

    tags$div(
      tags$span(style = "color: #e74c3c; font-weight: bold;",
                icon("arrow-up"), paste0(" UP: ", up_count)),
      tags$span(" | "),
      tags$span(style = "color: #3498db; font-weight: bold;",
                icon("arrow-down"), paste0(" DOWN: ", down_count))
    )
  })

  # ===========================================================================
  # BUG_008 FIX: Gene count displays for Tab 2 (Segment)
  # ===========================================================================
  output$seg_gene_counts <- renderUI({
    species <- input$seg_species
    segment <- input$seg_segment

    df <- switch(paste(species, segment, sep = "_"),
                 "Arabidopsis_seg1" = arab_seg1,
                 "Arabidopsis_seg2" = arab_seg2,
                 "Arabidopsis_seg3" = arab_seg3,
                 "Pennycress_seg1" = penny_seg1,
                 "Pennycress_seg2" = penny_seg2,
                 "Pennycress_seg3" = penny_seg3,
                 "Soybean_seg1" = soy_seg1,
                 "Soybean_seg2" = soy_seg2,
                 "Soybean_seg3" = soy_seg3)

    if (is.null(df) || nrow(df) == 0) {
      return(tags$div(style = "color: #999;", "No data"))
    }

    # Filter by current log2FC range and padj
    filtered <- df[!is.na(df$padj) & df$padj < input$seg_padj &
                     df$log2FoldChange >= input$seg_lfc[1] &
                     df$log2FoldChange <= input$seg_lfc[2], ]

    up_count <- sum(filtered$DEG_status == "UP in weep", na.rm = TRUE)
    down_count <- sum(filtered$DEG_status == "DOWN in weep", na.rm = TRUE)

    tags$div(
      tags$span(style = "color: #e74c3c; font-weight: bold;",
                icon("arrow-up"), paste0(" UP: ", up_count)),
      tags$span(" | "),
      tags$span(style = "color: #3498db; font-weight: bold;",
                icon("arrow-down"), paste0(" DOWN: ", down_count))
    )
  })

  # ===========================================================================
  # Tab 1: HPI Volcano Plot
  # ===========================================================================
  output$hpi_volcano <- renderPlotly({
    # BUG_007 FIX: Add req() for proper validation
    req(input$hpi_species, input$hpi_status)

    df <- hpi_data()
    req(nrow(df) > 0)  # Ensures we have data before proceeding

    # Ensure required columns exist
    if (!"log2FoldChange" %in% colnames(df) || !"pvalue" %in% colnames(df)) {
      return(plotly_empty())
    }

    # BUG_012 FIX: Filter data to only show genes within log2FC range
    # This ensures slider from 0 to 2 only shows upregulated genes
    df <- df[df$log2FoldChange >= input$hpi_lfc[1] & df$log2FoldChange <= input$hpi_lfc[2], ]

    if (nrow(df) == 0) {
      return(plotly_empty() %>% layout(title = "No genes in selected log2FC range"))
    }

    # Color based on padj significance threshold
    df$significant <- ifelse(!is.na(df$padj) & df$padj < input$hpi_padj,
                             df$DEG_status, "NOT Significant")

    # Create tooltip - include Gene_Name and OG if available
    df$tooltip <- paste0("Gene: ", df$gene_id,
                         "\nlog2FC: ", round(df$log2FoldChange, 2),
                         "\npadj: ", signif(df$padj, 3),
                         "\nStatus: ", df$DEG_status)

    if ("Gene_Name" %in% colnames(df)) {
      df$tooltip <- paste0(df$tooltip, "\nName: ", df$Gene_Name)
    }
    if ("Segment" %in% colnames(df)) {
      df$tooltip <- paste0(df$tooltip, "\nSegment: ", df$Segment)
    }
    if ("X.OG" %in% colnames(df)) {
      df$tooltip <- paste0(df$tooltip, "\nOrthogroup: ", df$X.OG)
    }

    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant, text = tooltip)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(values = deg_colors, name = "Status") +
      geom_vline(xintercept = c(input$hpi_lfc[1], input$hpi_lfc[2]), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(input$hpi_padj), linetype = "dashed", color = "grey50") +
      labs(title = paste(input$hpi_species, "- HPI Candidates (weep vs WT)"),
           x = "log2 Fold Change (weep/WT)",
           y = "-log10(p-value)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "right")

    ggplotly(p, tooltip = "text") %>% layout(legend = list(orientation = "v"))
  })
  
  # ===========================================================================
  # Tab 1: HPI Table
  # ===========================================================================
  output$hpi_table <- DT::renderDataTable({
    df <- hpi_data()
    if (nrow(df) == 0) return(DT::datatable(data.frame(Message = "No data available")))

    # BUG_011 FIX: Filter by log2FC range and padj thresholds
    df <- df[!is.na(df$padj) & df$padj < input$hpi_padj &
               df$log2FoldChange >= input$hpi_lfc[1] &
               df$log2FoldChange <= input$hpi_lfc[2], ]

    if (nrow(df) == 0) return(DT::datatable(data.frame(Message = "No genes match current filters")))

    # Select display columns
    display_cols <- c("gene_id", "Gene_Name", "baseMean", "log2FoldChange", "padj", "DEG_status")
    if ("Segment" %in% colnames(df)) display_cols <- c(display_cols, "Segment")
    if ("HPI" %in% colnames(df)) display_cols <- c(display_cols, "HPI")
    if ("OG" %in% colnames(df)) display_cols <- c(display_cols, "OG")

    display_cols <- display_cols[display_cols %in% colnames(df)]
    df <- df[, display_cols]

    # Format numeric columns
    if ("baseMean" %in% colnames(df)) df$baseMean <- round(df$baseMean, 1)
    if ("log2FoldChange" %in% colnames(df)) df$log2FoldChange <- round(df$log2FoldChange, 3)
    if ("padj" %in% colnames(df)) df$padj <- signif(df$padj, 3)

    DT::datatable(df, options = list(pageLength = 15, scrollX = TRUE),
                  filter = "top", rownames = FALSE)
  })

  output$hpi_summary <- renderText({
    df <- hpi_data()
    if (nrow(df) == 0) return("No data")

    # BUG_011 FIX: Filter by log2FC range and padj thresholds
    filtered_df <- df[!is.na(df$padj) & df$padj < input$hpi_padj &
                        df$log2FoldChange >= input$hpi_lfc[1] &
                        df$log2FoldChange <= input$hpi_lfc[2], ]

    total <- nrow(filtered_df)
    up <- sum(filtered_df$DEG_status == "UP in weep", na.rm = TRUE)
    down <- sum(filtered_df$DEG_status == "DOWN in weep", na.rm = TRUE)

    # Check available annotation columns
    has_gene_name <- "Gene_Name" %in% colnames(df)
    has_og <- "X.OG" %in% colnames(df)

    annotations <- c()
    if (has_gene_name) annotations <- c(annotations, "Gene Names")
    if (has_og) annotations <- c(annotations, "Orthogroups")
    annot_text <- if (length(annotations) > 0) paste0("\nAnnotations: ", paste(annotations, collapse = ", ")) else ""

    paste0("Species: ", input$hpi_species, "\n",
           "Filters: padj < ", input$hpi_padj, ", log2FC: [", input$hpi_lfc[1], ", ", input$hpi_lfc[2], "]\n",
           "Filtered genes: ", total, " (UP: ", up, ", DOWN: ", down, ")",
           annot_text)
  })
  
  output$download_hpi <- downloadHandler(
    filename = function() paste0(input$hpi_species, "_HPI_filtered_", Sys.Date(), ".csv"),
    content = function(file) {
      # BUG_011 FIX: Download filtered data based on current slider/filter settings
      df <- hpi_data()
      df <- df[!is.na(df$padj) & df$padj < input$hpi_padj &
                 df$log2FoldChange >= input$hpi_lfc[1] &
                 df$log2FoldChange <= input$hpi_lfc[2], ]
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # ===========================================================================
  # Tab 2: Segment Volcano Plot
  # ===========================================================================
  output$seg_volcano <- renderPlotly({
    # BUG_007 FIX: Add req() for proper validation
    req(input$seg_species, input$seg_segment, input$seg_status)

    df <- segment_data()
    req(nrow(df) > 0)  # Ensures we have data before proceeding

    # BUG_012 FIX: Filter data to only show genes within log2FC range
    df <- df[df$log2FoldChange >= input$seg_lfc[1] & df$log2FoldChange <= input$seg_lfc[2], ]

    if (nrow(df) == 0) {
      return(plotly_empty() %>% layout(title = "No genes in selected log2FC range"))
    }

    # Color based on padj significance threshold
    df$significant <- ifelse(!is.na(df$padj) & df$padj < input$seg_padj,
                             df$DEG_status, "NOT Significant")

    df$tooltip <- paste0("Gene: ", df$gene_id,
                         "\nlog2FC: ", round(df$log2FoldChange, 2),
                         "\npadj: ", signif(df$padj, 3),
                         "\nStatus: ", df$DEG_status)

    # BUG_009 FIX: Use descriptive segment label with tissue type
    seg_label <- get_seg_label(input$seg_segment)

    p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant, text = tooltip)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = deg_colors, name = "Status") +
      geom_vline(xintercept = c(input$seg_lfc[1], input$seg_lfc[2]), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(input$seg_padj), linetype = "dashed", color = "grey50") +
      labs(title = paste(input$seg_species, "-", seg_label, "(weep vs WT)"),
           x = "log2 Fold Change (weep/WT)",
           y = "-log10(p-value)") +
      theme_bw(base_size = 12)

    ggplotly(p, tooltip = "text")
  })
  
  # ===========================================================================
  # Tab 2: Segment Table
  # ===========================================================================
  output$seg_table <- DT::renderDataTable({
    df <- segment_data()
    if (nrow(df) == 0) return(DT::datatable(data.frame(Message = "No data available")))

    # BUG_008 FIX: Use range slider values
    # Filter for significant genes if requested
    if (input$seg_status == "All") {
      df <- df[!is.na(df$padj) & df$padj < input$seg_padj &
                 df$log2FoldChange >= input$seg_lfc[1] &
                 df$log2FoldChange <= input$seg_lfc[2], ]
    }

    display_cols <- c("gene_id", "Gene_Name", "baseMean", "log2FoldChange", "padj", "DEG_status")
    if ("HPI" %in% colnames(df)) display_cols <- c(display_cols, "HPI")
    if ("OG" %in% colnames(df)) display_cols <- c(display_cols, "OG")
    display_cols <- display_cols[display_cols %in% colnames(df)]
    df <- df[, display_cols]

    if ("baseMean" %in% colnames(df)) df$baseMean <- round(df$baseMean, 1)
    if ("log2FoldChange" %in% colnames(df)) df$log2FoldChange <- round(df$log2FoldChange, 3)
    if ("padj" %in% colnames(df)) df$padj <- signif(df$padj, 3)

    DT::datatable(df, options = list(pageLength = 15, scrollX = TRUE),
                  filter = "top", rownames = FALSE)
  })

  output$seg_summary <- renderText({
    df <- segment_data()
    if (nrow(df) == 0) return("No data")

    # BUG_008 FIX: Use range slider values
    sig_df <- df[!is.na(df$padj) & df$padj < input$seg_padj &
                   df$log2FoldChange >= input$seg_lfc[1] &
                   df$log2FoldChange <= input$seg_lfc[2], ]
    up <- sum(sig_df$DEG_status == "UP in weep", na.rm = TRUE)
    down <- sum(sig_df$DEG_status == "DOWN in weep", na.rm = TRUE)

    # BUG_009 FIX: Use descriptive segment label
    paste0("Species: ", input$seg_species, " | ", get_seg_label(input$seg_segment, full = TRUE), "\n",
           "Thresholds: padj < ", input$seg_padj, ", log2FC range: [", input$seg_lfc[1], ", ", input$seg_lfc[2], "]\n",
           "Significant DEGs: ", nrow(sig_df), " (UP: ", up, ", DOWN: ", down, ")")
  })
  
  output$download_seg <- downloadHandler(
    filename = function() paste0(input$seg_species, "_", input$seg_segment, "_filtered_", Sys.Date(), ".csv"),
    content = function(file) write.csv(segment_data(), file, row.names = FALSE)
  )
  
  # ===========================================================================
  # Tab 3: Comparison Data
  # ===========================================================================
  comp_data <- reactive({
    species <- input$comp_species
    comparison <- input$comp_segments

    # Get segment data for comparison
    seg_a <- switch(paste(species, substr(comparison, 1, 4), sep = "_"),
                    "Arabidopsis_seg1" = arab_seg1,
                    "Arabidopsis_seg2" = arab_seg2,
                    "Pennycress_seg1" = penny_seg1,
                    "Pennycress_seg2" = penny_seg2,
                    "Soybean_seg1" = soy_seg1,
                    "Soybean_seg2" = soy_seg2)

    seg_b <- switch(paste(species, substr(comparison, 9, 12), sep = "_"),
                    "Arabidopsis_seg2" = arab_seg2,
                    "Arabidopsis_seg3" = arab_seg3,
                    "Pennycress_seg2" = penny_seg2,
                    "Pennycress_seg3" = penny_seg3,
                    "Soybean_seg2" = soy_seg2,
                    "Soybean_seg3" = soy_seg3)

    # Filter by HPI genes only if checkbox is checked
    if (isTRUE(input$comp_hpi_filter)) {
      if (!is.null(seg_a) && "HPI" %in% colnames(seg_a)) {
        seg_a <- seg_a[seg_a$HPI == "Yes", ]
      }
      if (!is.null(seg_b) && "HPI" %in% colnames(seg_b)) {
        seg_b <- seg_b[seg_b$HPI == "Yes", ]
      }
    }

    list(seg_a = seg_a, seg_b = seg_b)
  })
  
  # ===========================================================================
  # Tab 3: Comparison Bar Plot
  # ===========================================================================
  output$comp_bar <- renderPlotly({
    data <- comp_data()
    if (is.null(data$seg_a) || is.null(data$seg_b)) return(plotly_empty())
    
    # Get significant genes
    sig_a <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj & 
                          data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]
    sig_b <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj & 
                          data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]
    
    # BUG_009 FIX: Use descriptive segment labels
    seg_codes <- strsplit(input$comp_segments, "_vs_")[[1]]
    seg_labels <- sapply(seg_codes, get_seg_label)

    df <- data.frame(
      Segment = c(seg_labels[1], seg_labels[1], seg_labels[2], seg_labels[2]),
      Direction = c("UP in weep", "DOWN in weep", "UP in weep", "DOWN in weep"),
      Count = c(sum(sig_a$DEG_status == "UP in weep", na.rm = TRUE),
                sum(sig_a$DEG_status == "DOWN in weep", na.rm = TRUE),
                sum(sig_b$DEG_status == "UP in weep", na.rm = TRUE),
                sum(sig_b$DEG_status == "DOWN in weep", na.rm = TRUE))
    )
    
    plot_ly(df, x = ~Segment, y = ~Count, color = ~Direction, type = "bar",
            colors = c("UP in weep" = "#e74c3c", "DOWN in weep" = "#3498db")) %>%
      layout(title = "DEGs per Segment",
             barmode = "group",
             xaxis = list(title = ""),
             yaxis = list(title = "Number of DEGs"))
  })
  
  # ===========================================================================
  # Tab 3: Venn Diagram for DEG Overlap
  # ===========================================================================
  output$comp_venn <- renderPlot({
    data <- comp_data()
    if (is.null(data$seg_a) || is.null(data$seg_b)) return(NULL)

    # Get significant gene IDs
    sig_a <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj &
                          data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]$gene_id
    sig_b <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj &
                          data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]$gene_id

    # BUG_009 FIX: Use descriptive segment labels
    seg_codes <- strsplit(input$comp_segments, "_vs_")[[1]]
    seg_labels <- sapply(seg_codes, get_seg_label)

    # Create named list for Venn diagram
    gene_lists <- list(sig_a, sig_b)
    names(gene_lists) <- seg_labels

    # Create Venn diagram
    ggVennDiagram(gene_lists,
                  label_alpha = 0,
                  category.names = seg_labels,
                  set_size = 5) +
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      scale_color_manual(values = c("#27ae60", "#e74c3c")) +
      labs(title = "DEG Overlap Between Segments") +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "none")
  })
  
  # ===========================================================================
  # Tab 3: Common Genes Table
  # ===========================================================================
  output$comp_common_table <- DT::renderDataTable({
    data <- comp_data()
    if (is.null(data$seg_a) || is.null(data$seg_b)) {
      return(DT::datatable(data.frame(Message = "No data available")))
    }

    # Get significant genes from both segments
    sig_a_df <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj &
                             data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]
    sig_b_df <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj &
                             data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]

    # Find common gene IDs
    common_ids <- intersect(sig_a_df$gene_id, sig_b_df$gene_id)

    if (length(common_ids) == 0) {
      return(DT::datatable(data.frame(Message = "No common genes found between segments")))
    }

    # BUG_009 FIX: Use descriptive segment labels for column naming
    seg_codes <- strsplit(input$comp_segments, "_vs_")[[1]]
    # Use short codes for column names (Seg1, Seg2, Seg3)
    seg_short <- gsub("seg", "Seg", seg_codes)

    # Create merged dataframe with info from both segments
    cols_to_get <- c("gene_id", "log2FoldChange", "padj", "DEG_status")
    if ("Gene_Name" %in% colnames(sig_a_df)) cols_to_get <- c(cols_to_get, "Gene_Name")
    if ("HPI" %in% colnames(sig_a_df)) cols_to_get <- c(cols_to_get, "HPI")
    if ("OG" %in% colnames(sig_a_df)) cols_to_get <- c(cols_to_get, "OG")
    cols_to_get <- cols_to_get[cols_to_get %in% colnames(sig_a_df)]

    common_a <- sig_a_df[sig_a_df$gene_id %in% common_ids, c("gene_id", "log2FoldChange", "padj", "DEG_status")]
    common_b <- sig_b_df[sig_b_df$gene_id %in% common_ids, c("gene_id", "log2FoldChange", "padj", "DEG_status")]

    # Get extra columns from seg_a for merged table
    extra_cols <- sig_a_df[sig_a_df$gene_id %in% common_ids, ]
    extra_info <- data.frame(gene_id = extra_cols$gene_id)
    if ("Gene_Name" %in% colnames(extra_cols)) extra_info$Gene_Name <- extra_cols$Gene_Name
    if ("HPI" %in% colnames(extra_cols)) extra_info$HPI <- extra_cols$HPI
    if ("OG" %in% colnames(extra_cols)) extra_info$OG <- extra_cols$OG

    colnames(common_a) <- c("gene_id", paste0("log2FC_", seg_short[1]), paste0("padj_", seg_short[1]), paste0("Status_", seg_short[1]))
    colnames(common_b) <- c("gene_id", paste0("log2FC_", seg_short[2]), paste0("padj_", seg_short[2]), paste0("Status_", seg_short[2]))

    # Merge by gene_id
    common_merged <- merge(common_a, common_b, by = "gene_id")

    # Add extra info columns (Gene_Name, HPI, OG)
    if (ncol(extra_info) > 1) {
      common_merged <- merge(common_merged, extra_info, by = "gene_id")
    }

    # Round numeric columns
    num_cols <- grep("log2FC|padj", colnames(common_merged), value = TRUE)
    for (col in num_cols) {
      if (grepl("log2FC", col)) {
        common_merged[[col]] <- round(common_merged[[col]], 3)
      } else {
        common_merged[[col]] <- signif(common_merged[[col]], 3)
      }
    }

    # BUG_009 FIX: Use descriptive labels in caption
    seg_labels_full <- sapply(seg_codes, function(x) get_seg_label(x, full = TRUE))

    DT::datatable(common_merged,
                  options = list(pageLength = 10, scrollX = TRUE),
                  filter = "top",
                  rownames = FALSE,
                  caption = paste("Common DE genes between", seg_labels_full[1], "and", seg_labels_full[2]))
  })

  output$comp_common_summary <- renderText({
    data <- comp_data()
    if (is.null(data$seg_a) || is.null(data$seg_b)) return("No data")

    sig_a_df <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj &
                             data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]
    sig_b_df <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj &
                             data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]

    common_ids <- intersect(sig_a_df$gene_id, sig_b_df$gene_id)

    # Count concordant vs discordant
    if (length(common_ids) > 0) {
      common_a <- sig_a_df[sig_a_df$gene_id %in% common_ids, c("gene_id", "DEG_status")]
      common_b <- sig_b_df[sig_b_df$gene_id %in% common_ids, c("gene_id", "DEG_status")]
      merged <- merge(common_a, common_b, by = "gene_id")
      concordant <- sum(merged$DEG_status.x == merged$DEG_status.y)
      discordant <- nrow(merged) - concordant
    } else {
      concordant <- 0
      discordant <- 0
    }

    paste0("Common genes: ", length(common_ids), "\n",
           "Concordant direction (same UP/DOWN): ", concordant, "\n",
           "Discordant direction (opposite): ", discordant)
  })

  # ===========================================================================
  # Tab 3: Comparison Tables
  # ===========================================================================
  output$comp_table_a <- DT::renderDataTable({
    data <- comp_data()
    if (is.null(data$seg_a)) return(DT::datatable(data.frame(Message = "No data")))
    
    df <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj & 
                       data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]
    
    if (input$comp_show_shared && !is.null(data$seg_b)) {
      sig_b <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj & 
                            data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]$gene_id
      df <- df[df$gene_id %in% sig_b, ]
    }
    
    display_cols <- c("gene_id", "Gene_Name", "log2FoldChange", "padj", "DEG_status", "HPI", "OG")
    display_cols <- display_cols[display_cols %in% colnames(df)]
    df <- df[, display_cols]
    df$log2FoldChange <- round(df$log2FoldChange, 3)
    df$padj <- signif(df$padj, 3)

    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  output$comp_table_b <- DT::renderDataTable({
    data <- comp_data()
    if (is.null(data$seg_b)) return(DT::datatable(data.frame(Message = "No data")))

    df <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj &
                       data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]

    if (input$comp_show_shared && !is.null(data$seg_a)) {
      sig_a <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj &
                            data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]$gene_id
      df <- df[df$gene_id %in% sig_a, ]
    }

    display_cols <- c("gene_id", "Gene_Name", "log2FoldChange", "padj", "DEG_status", "HPI", "OG")
    display_cols <- display_cols[display_cols %in% colnames(df)]
    df <- df[, display_cols]
    df$log2FoldChange <- round(df$log2FoldChange, 3)
    df$padj <- signif(df$padj, 3)

    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })
  
  output$comp_summary <- renderText({
    data <- comp_data()
    if (is.null(data$seg_a) || is.null(data$seg_b)) return("No data available")
    
    sig_a <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj & 
                          data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]$gene_id
    sig_b <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj & 
                          data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]$gene_id
    
    shared <- length(intersect(sig_a, sig_b))
    only_a <- length(setdiff(sig_a, sig_b))
    only_b <- length(setdiff(sig_b, sig_a))
    
    # BUG_009 FIX: Use descriptive segment labels
    seg_codes <- strsplit(input$comp_segments, "_vs_")[[1]]
    seg_labels <- sapply(seg_codes, get_seg_label)

    paste0("=== Segment Comparison: ", input$comp_species, " ===\n\n",
           seg_labels[1], ": ", length(sig_a), " DEGs\n",
           seg_labels[2], ": ", length(sig_b), " DEGs\n\n",
           "Shared between segments: ", shared, " DEGs\n",
           seg_labels[1], " only: ", only_a, " DEGs\n",
           seg_labels[2], " only: ", only_b, " DEGs\n\n",
           "Overlap percentage: ", round(shared / max(1, length(union(sig_a, sig_b))) * 100, 1), "%")
  })
  
  output$download_comp <- downloadHandler(
    filename = function() paste0(input$comp_species, "_", input$comp_segments, "_comparison_", Sys.Date(), ".csv"),
    content = function(file) {
      data <- comp_data()
      sig_a <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj &
                            data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]
      sig_b <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj &
                            data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]
      sig_a$Comparison_Segment <- "Segment_A"
      sig_b$Comparison_Segment <- "Segment_B"
      combined <- rbind(sig_a[, intersect(colnames(sig_a), colnames(sig_b))],
                        sig_b[, intersect(colnames(sig_a), colnames(sig_b))])
      write.csv(combined, file, row.names = FALSE)
    }
  )

  # ===========================================================================
  # Tab 3: Download Common Genes
  # ===========================================================================
  output$download_common <- downloadHandler(
    filename = function() paste0(input$comp_species, "_", input$comp_segments, "_common_genes_", Sys.Date(), ".csv"),
    content = function(file) {
      data <- comp_data()

      # Get significant genes from both segments
      sig_a_df <- data$seg_a[!is.na(data$seg_a$padj) & data$seg_a$padj < input$comp_padj &
                               data$seg_a$log2FoldChange >= input$comp_lfc[1] & data$seg_a$log2FoldChange <= input$comp_lfc[2], ]
      sig_b_df <- data$seg_b[!is.na(data$seg_b$padj) & data$seg_b$padj < input$comp_padj &
                               data$seg_b$log2FoldChange >= input$comp_lfc[1] & data$seg_b$log2FoldChange <= input$comp_lfc[2], ]

      # Find common gene IDs
      common_ids <- intersect(sig_a_df$gene_id, sig_b_df$gene_id)

      # BUG_009 FIX: Use short codes for column naming
      seg_codes <- strsplit(input$comp_segments, "_vs_")[[1]]
      seg_short <- gsub("seg", "Seg", seg_codes)

      # Create merged dataframe with info from both segments
      common_a <- sig_a_df[sig_a_df$gene_id %in% common_ids, c("gene_id", "log2FoldChange", "padj", "DEG_status")]
      common_b <- sig_b_df[sig_b_df$gene_id %in% common_ids, c("gene_id", "log2FoldChange", "padj", "DEG_status")]

      colnames(common_a) <- c("gene_id", paste0("log2FC_", seg_short[1]), paste0("padj_", seg_short[1]), paste0("Status_", seg_short[1]))
      colnames(common_b) <- c("gene_id", paste0("log2FC_", seg_short[2]), paste0("padj_", seg_short[2]), paste0("Status_", seg_short[2]))

      # Merge by gene_id
      common_merged <- merge(common_a, common_b, by = "gene_id")

      write.csv(common_merged, file, row.names = FALSE)
    }
  )

  # ===========================================================================
  # Tab 4: Cross-Species Comparison - Reactive Data
  # ===========================================================================

  # Reactive: Cross-species data (all segments combined)
  cross_species_data <- reactive({
    if (is.null(all_species_deg)) return(data.frame())

    data <- all_species_deg

    # Remove rows without OG
    data <- data[!is.na(data$OG) & data$OG != "", ]

    # Filter by HPI genes only if checkbox is checked
    if (isTRUE(input$cross_hpi_filter) && "HPI" %in% colnames(data)) {
      data <- data[data$HPI == "Yes", ]
    }

    return(data)
  })

  # Reactive: Group by Orthogroup
  og_grouped <- reactive({
    data <- cross_species_data()
    if (nrow(data) == 0) return(data.frame())

    # Group by OG and summarize
    og_summary <- data %>%
      group_by(OG) %>%
      summarise(
        # Count species
        Num_Species = n_distinct(Species),
        Species_List = paste(sort(unique(Species)), collapse = ", "),

        # Arabidopsis info
        Arab_gene = paste(unique(gene_id[Species == "Arabidopsis"]), collapse = "; "),
        Arab_LFC = ifelse(any(Species == "Arabidopsis"),
                          paste(round(log2FoldChange[Species == "Arabidopsis"], 2), collapse = "; "),
                          "-"),
        Arab_status = ifelse(any(Species == "Arabidopsis"),
                             paste(unique(DEG_status[Species == "Arabidopsis"]), collapse = "; "),
                             "-"),

        # Pennycress info
        Penny_gene = paste(unique(gene_id[Species == "Pennycress"]), collapse = "; "),
        Penny_LFC = ifelse(any(Species == "Pennycress"),
                           paste(round(log2FoldChange[Species == "Pennycress"], 2), collapse = "; "),
                           "-"),
        Penny_status = ifelse(any(Species == "Pennycress"),
                              paste(unique(DEG_status[Species == "Pennycress"]), collapse = "; "),
                              "-"),

        # Soybean info
        Soy_gene = paste(unique(gene_id[Species == "Soybean"]), collapse = "; "),
        Soy_LFC = ifelse(any(Species == "Soybean"),
                         paste(round(log2FoldChange[Species == "Soybean"], 2), collapse = "; "),
                         "-"),
        Soy_status = ifelse(any(Species == "Soybean"),
                            paste(unique(DEG_status[Species == "Soybean"]), collapse = "; "),
                            "-"),

        # Direction info
        Num_UP = sum(DEG_status == "UP in weep"),
        Num_DOWN = sum(DEG_status == "DOWN in weep"),

        # Gene name (from Arabidopsis HPI)
        Gene_Name = paste(unique(Gene_Name[!is.na(Gene_Name)]), collapse = "; "),

        .groups = "drop"
      ) %>%
      mutate(
        # Determine direction consistency
        Direction = case_when(
          Num_UP > 0 & Num_DOWN == 0 ~ "ALL UP in weep",
          Num_DOWN > 0 & Num_UP == 0 ~ "ALL DOWN in weep",
          TRUE ~ "MIXED"
        ),
        Is_Consistent = Direction != "MIXED"
      )

    return(og_summary)
  })

  # ===========================================================================
  # Tab 4: Value Boxes
  # ===========================================================================
  output$total_three_species_box <- renderValueBox({
    og_data <- og_grouped()
    count <- sum(og_data$Num_Species == 3)
    valueBox(
      count,
      "Conserved OGs (All 3 Species)",
      icon = icon("layer-group"),
      color = "green"
    )
  })

  output$total_two_species_box <- renderValueBox({
    og_data <- og_grouped()
    count <- sum(og_data$Num_Species == 2)
    valueBox(
      count,
      "Shared OGs (2 Species)",
      icon = icon("code-branch"),
      color = "yellow"
    )
  })

  # ===========================================================================
  # Tab 4: Three Species Table
  # ===========================================================================
  output$three_species_table <- DT::renderDataTable({
    og_data <- og_grouped()
    if (nrow(og_data) == 0) return(DT::datatable(data.frame(Message = "No data available")))

    # Filter to 3 species only
    three_sp <- og_data[og_data$Num_Species == 3, ]

    if (nrow(three_sp) == 0) return(DT::datatable(data.frame(Message = "No OGs found in all 3 species")))

    # Select columns to display
    display_data <- three_sp %>%
      select(OG, Gene_Name, Direction,
             Arab_gene, Arab_LFC, Arab_status,
             Penny_gene, Penny_LFC, Penny_status,
             Soy_gene, Soy_LFC, Soy_status)

    DT::datatable(display_data,
                  filter = "top",
                  options = list(
                    pageLength = 25,
                    scrollX = TRUE
                  ),
                  rownames = FALSE) %>%
      formatStyle(
        'Direction',
        backgroundColor = styleEqual(
          c("ALL UP in weep", "ALL DOWN in weep", "MIXED"),
          c("#d4edda", "#f8d7da", "#fff3cd")
        )
      ) %>%
      formatStyle(
        'Arab_status',
        color = styleEqual(
          c("UP in weep", "DOWN in weep", "-"),
          c("#e74c3c", "#3498db", "#6c757d")
        )
      ) %>%
      formatStyle(
        'Penny_status',
        color = styleEqual(
          c("UP in weep", "DOWN in weep", "-"),
          c("#e74c3c", "#3498db", "#6c757d")
        )
      ) %>%
      formatStyle(
        'Soy_status',
        color = styleEqual(
          c("UP in weep", "DOWN in weep", "-"),
          c("#e74c3c", "#3498db", "#6c757d")
        )
      )
  })

  # ===========================================================================
  # Tab 4: Two Species Table
  # ===========================================================================
  output$two_species_table <- DT::renderDataTable({
    og_data <- og_grouped()
    if (nrow(og_data) == 0) return(DT::datatable(data.frame(Message = "No data available")))

    # Filter to 2 species only
    two_sp <- og_data[og_data$Num_Species == 2, ]

    if (nrow(two_sp) == 0) return(DT::datatable(data.frame(Message = "No OGs found in exactly 2 species")))

    # Select columns to display
    display_data <- two_sp %>%
      select(OG, Gene_Name, Species_List, Direction,
             Arab_gene, Arab_LFC, Arab_status,
             Penny_gene, Penny_LFC, Penny_status,
             Soy_gene, Soy_LFC, Soy_status)

    DT::datatable(display_data,
                  filter = "top",
                  options = list(
                    pageLength = 25,
                    scrollX = TRUE
                  ),
                  rownames = FALSE) %>%
      formatStyle(
        'Direction',
        backgroundColor = styleEqual(
          c("ALL UP in weep", "ALL DOWN in weep", "MIXED"),
          c("#d4edda", "#f8d7da", "#fff3cd")
        )
      ) %>%
      formatStyle(
        'Arab_status',
        color = styleEqual(
          c("UP in weep", "DOWN in weep", "-"),
          c("#e74c3c", "#3498db", "#6c757d")
        )
      ) %>%
      formatStyle(
        'Penny_status',
        color = styleEqual(
          c("UP in weep", "DOWN in weep", "-"),
          c("#e74c3c", "#3498db", "#6c757d")
        )
      ) %>%
      formatStyle(
        'Soy_status',
        color = styleEqual(
          c("UP in weep", "DOWN in weep", "-"),
          c("#e74c3c", "#3498db", "#6c757d")
        )
      )
  })

  # ===========================================================================
  # Tab 4: Download Handlers
  # ===========================================================================
  output$download_three_species <- downloadHandler(
    filename = function() {
      paste0("OGs_three_species_", Sys.Date(), ".csv")
    },
    content = function(file) {
      og_data <- og_grouped()
      three_sp <- og_data[og_data$Num_Species == 3, ]
      write.csv(three_sp, file, row.names = FALSE)
    }
  )

  output$download_two_species <- downloadHandler(
    filename = function() {
      paste0("OGs_two_species_", Sys.Date(), ".csv")
    },
    content = function(file) {
      og_data <- og_grouped()
      two_sp <- og_data[og_data$Num_Species == 2, ]
      write.csv(two_sp, file, row.names = FALSE)
    }
  )
}

# =============================================================================
# Run App
# =============================================================================

shinyApp(ui = ui, server = server)
