################################################################################
#@@@@@@@ REQUIRED LIBRARIES
################################################################################

library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(shinyWidgets) 
library(DT)
library(readr)
library(stringr)
library(tidyverse)

# Load pre-processed table for introgression and the gene info
results_list <- readRDS("results_list.rds")
gene_info <- readRDS("gene_info.rds")
#results_list[1:2]
results_list_rename <- readRDS("results_list_rename.rds")

# Function to filter samples based on gene name for introgressed regions
find_introgressed_samples_by_gene <- function(gene_info, results_list_rename, gene_name) {
  gene_coords <- gene_info[gene_info$gene_id == gene_name, ]
  if(nrow(gene_coords) == 0) return("Gene not found")
  
  chromosome <- gene_coords$chromosome[1]
  start_position <- gene_coords$start[1]
  end_position <- gene_coords$end[1]
  
  introgressed_samples <- names(results_list_rename)[sapply(results_list_rename, function(sample_data) {
    any(sample_data$V1 == paste0("chr", chromosome) & sample_data$V4 == "Introgression" &
          ((sample_data$V2 >= start_position & sample_data$V2 <= end_position) |
             (sample_data$V3 >= start_position & sample_data$V3 <= end_position) |
             (sample_data$V2 <= start_position & sample_data$V3 >= end_position)))
  })]
  
  return(introgressed_samples)
}

# Function to filter samples based on chromosome, start, and end position for introgressed regions
find_introgressed_samples_by_region <- function(results_list_rename, chromosome, start_position, end_position) {
  introgressed_samples <- names(results_list_rename)[sapply(results_list_rename, function(sample_data) {
    any(sample_data$V1 == paste0("chr", chromosome) & sample_data$V4 == "Introgression" &
          ((sample_data$V2 >= start_position & sample_data$V2 <= end_position) |
             (sample_data$V3 >= start_position & sample_data$V3 <= end_position)))
  })]
  
  return(introgressed_samples)
}

# Function to find the samples by Gene Name
check_introgression_in_gene <- function(sample_data, chromosome, gene_start, gene_end) {
  # Ensure chromosome names match (e.g., 'chr1' vs '1')
  chromosome <- paste0("chr", chromosome)
  
  introgressed <- sample_data %>%
    filter(V1 == chromosome, V4 == "Introgression") %>%
    filter((V2 <= gene_end & V2 >= gene_start) | (V3 >= gene_start & V3 <= gene_end) | (V2 <= gene_start & V3 >= gene_end))
  
  return(nrow(introgressed) > 0)
}


ui <- dashboardPage(
  dashboardHeader(title = "BZea Introgression Sample Finder and Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Gene Search", tabName = "gene_search", icon = icon("search")),
      menuItem("Introgression Analysis", tabName = "introgression_analysis", icon = icon("chart-bar"))
    )
  ),
  dashboardBody(
    tabItems(
      # Gene Search Tab
      tabItem(tabName = "gene_search",
              fluidRow(
                box(title = "Gene Search",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    width = 12,
                    textInput("gene_name", "Enter Gene Name:", value = "Zm00001eb121780"),
                    textInput("chromosome", "Enter Chromosome:"),
                    textInput("start_position", "Enter Start Position:"),
                    textInput("end_position", "Enter End Position:"),
                    actionButton("submit", "Submit"),
                    h3("Samples with Introgressions in the Region or Gene:"),
                    verbatimTextOutput("sample_output"),
                    selectInput("gene_search_sample_select", "Select Sample:", choices = NULL),
                    plotOutput("gene_search_introgression_plot")
                )
              )
      ),
      # Introgression Analysis Tab
      tabItem(tabName = "introgression_analysis",
              fluidRow(
                box(title = "Select PN Group:",
                    selectInput("pn_group", "Select PN Group:", choices = NULL),
                    actionButton("update", "Update"),
                    selectizeInput("sample_select", "Select Sample:", choices = NULL),
                    width = 12
                )
              ),
              fluidRow(
                column(12, dataTableOutput("summary_stats"))
              ),
              fluidRow(
                column(6, plotOutput("introgression_plot")),
                column(6, dataTableOutput("introgression_details"))
              ),
              fluidRow(
                column(12, box(textOutput("introgression_percentage"), width = 12))
              )
      )
    )
  )
)


server <- function(input, output, session) {
  
  # Defining the reactive values to store the selected samples and groups
  selected_samples <- reactiveVal()
  selected_group <- reactiveVal()
  
  # Reactivity for gene-based search
  observeEvent(input$submit, {
    gene_name <- input$gene_name
    chromosome <- input$chromosome
    start_pos <- as.numeric(input$start_position)
    end_pos <- as.numeric(input$end_position)
    
    if (!is.na(start_pos) && !is.na(end_pos) && chromosome != "") {
      introgressed_samples <- find_introgressed_samples_by_region(results_list_rename, chromosome, start_pos, end_pos)
    } else if (gene_name != "") {
      introgressed_samples <- find_introgressed_samples_by_gene(gene_info, results_list_rename, gene_name)
    } else {
      introgressed_samples <- "Please enter valid search criteria."
    }
    
    # Update the reactive values based on the search
    selected_samples(introgressed_samples)
    selected_group(names(introgressed_samples)) # Assuming introgressed_samples is a named list
    
    # Update the UI elements
    updateSelectInput(session, "pn_group", choices = names(grouped_results_list))
    updateSelectizeInput(session, "sample_select", choices = names(introgressed_samples))
    
    # Make sure to trigger re-render of the output after search
    output$sample_output <- renderPrint({
      selected_samples()
    })
  })
  
  # Group the results by the PN group
  patterns <- c("Zlux", "Bals", "Chal", "Hueh", "Mesa", "Nabo", "Zdip", "Dura")
  # Filter and split the list based on the presence of specific substrings in the names
  grouped_results_list <- lapply(patterns, function(pat) {
    filtered_names <- names(results_list_rename)[grep(pat, names(results_list_rename))]
    if (length(filtered_names) > 0) {
      results_list_rename[filtered_names]
    } else {
      NULL
    }
  })
  names(grouped_results_list) <- patterns
  
  # Remove NULL elements if any patterns were not found
  grouped_results_list <- grouped_results_list[sapply(grouped_results_list, length) > 0]
  
  observe({
    updateSelectInput(session, "pn_group", choices = names(grouped_results_list))
  })
  
  observeEvent(input$update, {
    if (!is.null(input$pn_group)) {
      updateSelectizeInput(session, "sample_select", choices = names(grouped_results_list[[input$pn_group]]))
    }
  })
  
  output$summary_stats <- renderDataTable({
    req(input$pn_group)
    group <- input$pn_group
    if (is.null(group) || is.null(grouped_results_list[[group]])) {
      return(data.frame())  # Return an empty dataframe if the group is NULL or not found
    }
    summary_stats <- lapply(grouped_results_list[[selected_group]], function(sample_df) {
      if (is.data.frame(sample_df) && "V1" %in% names(sample_df)) {
        introgressions <- sample_df %>%
          filter(V4 == "Introgression") %>%
          mutate(size = V3 - V2) %>%
          group_by(V1) %>%
          summarise(total_number_of_introgressions = n(),
                    average_size_introgression_mb = mean(size / 1000000))
        return(introgressions)
      } else {
        return(data.frame(V1 = character(), total_number_of_introgressions = integer(),
                          average_size_introgression_mb = numeric()))
      }
    }) %>%
      bind_rows()
    
    return(summary_stats)
  })
  
  output$introgression_plot <- renderPlot({
    req(input$sample_select, input$pn_group)
    selected_group <- input$pn_group
    selected_sample <- input$sample_select
    
    if (selected_sample %in% names(grouped_results_list[[selected_group]])) {
      sample_data <- grouped_results_list[[selected_group]][[selected_sample]] %||% data.frame()
      
      # Assuming sample_data is a data frame with the expected structure
      if (is.data.frame(sample_data) && all(c("V1", "V2", "V3", "V4") %in% names(sample_data))) {
        sample_data$V1 <- factor(sample_data$V1, levels = unique(sample_data$V1))
        sample_data$V4 <- factor(sample_data$V4, levels = c("B73", "Introgression"))
        
        # Plot
        ggplot(sample_data, aes(x = V2, xend = V3, y = V1, yend = V1, color = V4)) +
          geom_segment(size = 10) +
          scale_color_manual(values = c("B73" = "blue", "Introgression" = "red")) +
          theme_minimal(base_size = 18) +
          theme(
            legend.position = "bottom",
            panel.grid.major = element_line(color = "gray", linetype = "dotted"),
            axis.text.x = element_text(angle = 45, hjust = 1)
          ) +
          labs(x = "Position along chromosome", y = "Chromosome", color = "Region") +
          ggtitle(paste("Introgression Locations of", selected_sample))
      }
    }
    
  })
  
  output$introgression_details <- renderDataTable({
    req(input$sample_select, input$pn_group)
    selected_group <- input$pn_group
    selected_sample <- input$sample_select
    
    sample_data <- grouped_results_list[[selected_group]][[selected_sample]] %||% data.frame()
    # Prepare and display introgression details
    introgression_details <- sample_data %>%
      filter(V4 == "Introgression") %>%
      mutate(Size_Mb = (V3 - V2) / 1000000) %>%
      select(Chromosome = V1, Introgression_start = V2, Introgression_end = V3, Size_Mb)
    
    datatable(introgression_details, options = list(pageLength = 5))
  })
  
  output$introgression_percentage <- renderText({
    req(input$sample_select, input$pn_group)
    selected_group <- input$pn_group
    selected_sample <- input$sample_select
    
    sample_data <- grouped_results_list[[selected_group]][[selected_sample]] %||% data.frame()
    
    # Calculate and display the percentage of introgression
    total_chromosome_length_mb <- 2131.847
    
    # Filter for introgressions only and calculate the total size
    introgression_total_size_mb <- sum(sample_data %>%
                                         dplyr::filter(V4 == "Introgression") %>%
                                         dplyr::mutate(Size_Mb = (V3 - V2) / 1000000) %>%
                                         .$Size_Mb)
    
    # Calculate the percentage of introgression
    introgression_percentage <- introgression_total_size_mb / total_chromosome_length_mb * 100
    
    # Output the % introgression
    paste("Total introgression percentage is:", round(introgression_percentage, 4), "%")
  })
  
  # For gene search - summary stats
  observe({
    updateSelectInput(session, "gene_search_sample_select",
                      choices = selected_samples(),
                      selected = selected_samples()[1])
  })
  
  # For gene search - introgression plot
  output$gene_search_introgression_plot <- renderPlot({
    req(selected_samples())  # Make sure there are selected samples
    selected_sample <- input$gene_search_sample_select
    
    # Extract the data for the selected sample
    sample_data <- results_list_rename[[selected_sample]] %||% data.frame()
    
    # Assuming sample_data is a data frame with the expected structure
    if (is.data.frame(sample_data) && all(c("V1", "V2", "V3", "V4") %in% names(sample_data))) {
      sample_data$V1 <- factor(sample_data$V1, levels = unique(sample_data$V1))
      sample_data$V4 <- factor(sample_data$V4, levels = c("B73", "Introgression"))
      
      # Plot
      ggplot(sample_data, aes(x = V2, xend = V3, y = V1, yend = V1, color = V4)) +
        geom_segment(size = 10) +
        scale_color_manual(values = c("B73" = "blue", "Introgression" = "red")) +
        theme_minimal(base_size = 18) +
        theme(
          legend.position = "bottom",
          panel.grid.major = element_line(color = "gray", linetype = "dotted"),
          axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(x = "Position along chromosome", y = "Chromosome", color = "Region") +
        ggtitle(paste("Introgression Locations of", selected_sample))
    }
  })
}


shinyApp(ui = ui, server = server)