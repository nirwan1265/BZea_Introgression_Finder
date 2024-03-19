library(shiny)
library(dplyr)
library(ggplot2)
library(DT)
library(readr)
library(stringr)
library(shinydashboard)

# Load pre-processed data outside of server function so it's only done once
#results_list_rename <- readRDS("/Users/nirwantandukar/Documents/Research/results/introgression_map/PN3/results_list_rename.rds")
gene_info <- readRDS("gene_info.rds")
results_list_rename <- readRDS("results_list_rename.rds")

# Create a new list that groups the results by the PN group
grouped_results_list_rename <- split(results_list_rename, gsub("_SID.*", "", names(results_list_rename)))

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


# UI
ui <- fluidPage(
  titlePanel("BZea Introgression Sample Finder"),
  sidebarLayout(
    sidebarPanel(
      textInput("gene_name", "Enter Gene Name:", value = "Zm00001eb121780"),
      textInput("chromosome", "Enter Chromosome:"),
      textInput("start_position", "Enter Start Position:"),
      textInput("end_position", "Enter End Position:"),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      h3("Samples with Introgressions in the Region or Gene:"),
      verbatimTextOutput("sample_output")
    )
  )
)

# Server logic
server <- function(input, output) {
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
    
    output$sample_output <- renderPrint({
      introgressed_samples
    })
  })
}



# Run the application
shinyApp(ui = ui, server = server)


