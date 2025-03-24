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

# Get the gene information
gene_info <- readRDS("gene_info.rds")

results_list_rename <- readRDS("results_list_new_name.rds")
head(results_list_rename)


# Step 1: Filter the results_list_rename to only include "Introgression" regions
results_list_rename_filtered <- lapply(results_list_rename, function(sample_data) {
  sample_data %>% filter(V4 == "Introgression")
})

# Step 2: Separate the filtered introgressions by chromosome
results_list_by_chr <- list()

# Loop through chromosome numbers (1 to 10 in your case, modify if needed)
for (chr_num in 1:10) {
  chr_label <- paste0("chr", chr_num)
  results_list_by_chr[[chr_label]] <- lapply(results_list_rename_filtered, function(sample_data) {
    sample_data %>% filter(V1 == chr_label)  # Filter by chromosome
  })
  
  # Remove any empty entries
  results_list_by_chr[[chr_label]] <- Filter(function(x) nrow(x) > 0, results_list_by_chr[[chr_label]])
}



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
# find_introgressed_samples_by_region <- function(results_list_rename, chromosome, start_position, end_position) {
#   introgressed_samples <- names(results_list_rename)[sapply(results_list_rename, function(sample_data) {
#     any(sample_data$V1 == paste0("chr", chromosome) & sample_data$V4 == "Introgression" &
#           ((sample_data$V2 >= start_position & sample_data$V2 <= end_position) |
#              (sample_data$V3 >= start_position & sample_data$V3 <= end_position)))
#   })]
# 
#   return(introgressed_samples)
# }
find_introgressed_samples_by_region <- function(results_list_rename, chromosome, start_position, end_position) {
  introgressed_samples <- names(results_list_rename)[sapply(results_list_rename, function(sample_data) {
    # Check if the introgression fully encompasses the specified region
    any(sample_data$V1 == paste0("chr", chromosome) & sample_data$V4 == "Introgression" &
          sample_data$V2 <= start_position & sample_data$V3 >= end_position)  # Encompasses the region
  })]
  
  return(introgressed_samples)
}



# Function to filter samples based on chromosome, start, and end position for introgressed regions
# Function to find continuous introgressions within a specified region
# Function to search introgressions by chromosome and position range

## Extended for selection
find_introgressed_samples_by_extended_gene <- function(gene_info, results_list_rename, gene_name, extension = 5000) {
  gene_coords <- gene_info[gene_info$gene_id == gene_name, ]
  if(nrow(gene_coords) == 0) return("Gene not found")
  
  chromosome <- gene_coords$chromosome[1]
  start_position <- max(gene_coords$start[1] - extension, 0)  # Ensure the start is not less than 0
  end_position <- gene_coords$end[1] + extension
  
  introgressed_samples <- names(results_list_rename)[sapply(results_list_rename, function(sample_data) {
    any(sample_data$V1 == paste0("chr", chromosome) & sample_data$V4 == "Introgression" &
          ((sample_data$V2 >= start_position & sample_data$V2 <= end_position) |
             (sample_data$V3 >= start_position & sample_data$V3 <= end_position) |
             (sample_data$V2 <= start_position & sample_data$V3 >= end_position)))
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


### PHENOTYPE
bzea_merged <- read.csv("Bzea_merged.csv", stringsAsFactors = FALSE, na.strings = c("#N/A")) %>%
  filter(!is.na(Lines))                # remove rows where Lines is NA

str(bzea_merged)
bzea_merged <- bzea_merged %>%
  mutate(Class = if_else(
    Lines == "B73",
    "B73",
    sub("\\..*", "", Lines)
  )) %>%
  mutate(Class = recode(
    Class,
    "Zh" = "Huehuetenango",
    "Zl" = "Luxurians",
    "Zx" = "Mexicana",
    "Zv" = "Parviglumis",
    "Zd" = "Diploperennis"
  ))

# Convert date columns
bzea_merged$DTS <- as.Date(bzea_merged$DTS, format = "%m/%d/%y")
bzea_merged$DTA <- as.Date(bzea_merged$DTA, format = "%m/%d/%y")

# Define variable categories
continuous_vars <- c("PH", "EH", "BW", "BL", "SL", "SPAD1", "SPAD2", "elevation", "latitude", "longitude")
date_vars <- c("DTS", "DTA")
categorical_vars <- c("Class","EN","ST", "StPi", "StPu", "Kinki", "Prolif", "NBR",
                      "location", "race", "country", "state", "county")

ui <- dashboardPage(
  dashboardHeader(title = "BZea Introgression Sample Finder and Analysis"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Gene Search", tabName = "gene_search", icon = icon("search")),
      menuItem("Introgression Analysis", tabName = "introgression_analysis", icon = icon("chart-bar")),
      # Main parent item (BZea Phenotypes)
      menuItem("BZea Phenotypes", tabName = "main", icon = icon("leaf"),
               # Sub-items (Table, Statistics, Figures)
               menuSubItem("Table", tabName = "table", icon = icon("table")),
               menuSubItem("Statistics", tabName = "stats", icon = icon("chart-line")),
               # Use menuItem (not menuSubItem) for Figures to allow nesting plots
               menuItem("Figures", tabName = "figures", icon = icon("bar-chart"),
                        menuSubItem("Histogram", tabName = "histogram"),
                        menuSubItem("Violin Plot", tabName = "violin"),
                        menuSubItem("Boxplot", tabName = "boxplot"),
                        menuSubItem("Density Plot", tabName = "density"),
                        menuSubItem("Scatter Plot", tabName = "scatter")
               )
      )
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
              ),
              fluidRow(
                box(title = "Extended Gene Region Samples",
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    width = 12,
                    dataTableOutput("extended_gene_region_samples")
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
                    pickerInput("sample_select", "Select Sample:", choices = NULL,
                                options = list(`actions-box` = TRUE, `live-search` = TRUE)), # This replaces selectizeInput
                    
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
      ),
      tabItem(tabName = "table",
              fluidRow(
                column(6, uiOutput("categorical_filters")),
                column(6, uiOutput("numeric_date_filters")),
                column(12, actionButton("reset_filters", "Reset Filters"), downloadButton('downloadData', 'Download CSV'))
              ),
              fluidRow(
                column(12, DT::dataTableOutput("filteredTable"))
              )
      ),
      tabItem(tabName = "stats",
              selectInput("stat_var", "Select a variable", choices = NULL),
              verbatimTextOutput("stats")
      ),
      tabItem(tabName = "histogram",
              selectInput("hist_var", "Select a variable", choices = NULL),
              plotOutput("hist")
      ),
      tabItem(tabName = "violin",
              selectInput("violin_var", "Select a variable", choices = NULL),
              plotOutput("violin")
      ),
      tabItem(tabName = "boxplot",
              selectInput("box_var", "Select a variable", choices = NULL),
              plotOutput("boxplot")
      ),
      tabItem(tabName = "density",
              selectInput("density_var", "Select a variable", choices = NULL),
              plotOutput("density")
      ),
      tabItem(tabName = "scatter",
              selectInput("scatter_var1", "Select variable for x-axis", choices = NULL),
              selectInput("scatter_var2", "Select variable for y-axis", choices = NULL),
              plotOutput("scatter")
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
  
  # 1. Define how each short code matches itself (because your new names *already* contain Zx, Zl, etc.)
  pattern_groups <- list(
    "Zl" = "Zl",  # will match names containing "Zl"
    "Zv" = "Zv",  # will match names containing "Zv"
    "Zx" = "Zx",  # will match names containing "Zx"
    "Zh" = "Zh",  # will match names containing "Zh"
    "Zd" = "Zd"   # will match names containing "Zd"
  )
  
  # 2. Define final descriptive labels for each group
  group_labels <- c(
    "Zl" = "Luxurians",
    "Zv" = "Parviglumis",
    "Zx" = "Mexicana",
    "Zh" = "Huehuetenango",
    "Zd" = "Diploperennis"
  )
  
  # 3. Create an empty list to store grouped results
  grouped_results_list <- list()
  
  # 4. Filter and split the list based on these short codes in the names
  for (group_code in names(pattern_groups)) {
    
    # Build a regex that matches the short code (e.g., "Zl", "Zx")
    search_pattern <- paste(pattern_groups[[group_code]], collapse = "|")
    
    # Find all sample names in 'results_list_rename' that contain these codes
    filtered_names <- names(results_list_rename)[grep(search_pattern, names(results_list_rename))]
    
    # If any matches are found, store them under the final descriptive label
    if (length(filtered_names) > 0) {
      group_name <- group_labels[group_code]  # e.g., "Luxurians", "Mexicana", etc.
      grouped_results_list[[group_name]] <- results_list_rename[filtered_names]
    }
  }
  
  # 5. Remove any empty groups
  grouped_results_list <- grouped_results_list[sapply(grouped_results_list, length) > 0]
  
  observe({
    updateSelectInput(session, "pn_group", choices = names(grouped_results_list))
  })
  
  observeEvent(input$update, {
    if (!is.null(input$pn_group)) {
      updateSelectizeInput(session, "sample_select", choices = names(grouped_results_list[[input$pn_group]]))
    }
  })
  
  observeEvent(input$pn_group, {
    req(input$pn_group)
    # This will update the 'sample_select' input with search functionality
    updatePickerInput(session, "sample_select",
                      choices = names(grouped_results_list[[input$pn_group]]),
                      options = list(`actions-box` = TRUE, `live-search` = TRUE)
    )
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
    
    # Extract the gene information
    gene_name <- input$gene_name
    gene_coords <- gene_info[gene_info$gene_id == gene_name, ]
    if (nrow(gene_coords) == 0) {
      stop("Gene not found")
    }
    
    # Get the chromosome and start position of the gene
    chromosome <- paste0("chr", gene_coords$chromosome)
    gene_start <- gene_coords$start
    
    # Extract the data for the selected sample
    sample_data <- results_list_rename[[selected_sample]] %||% data.frame()
    
    # Plot introgressions
    p <- ggplot(sample_data, aes(x = V2, xend = V3, y = V1, yend = V1, color = V4)) +
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
    
    # Subset for the chromosome of the gene
    chromosome_data <- sample_data[sample_data$V1 == chromosome, ]
    
    # Check if there is data for the chromosome
    if (nrow(chromosome_data) > 0) {
      # Add a thick vertical line for the gene start position only on the gene's chromosome
      p <- p + geom_vline(aes(xintercept = gene_start),
                          color = "green", size = 1.5, linetype = "longdash")
      
      # Optionally, add a label for the gene start line
      p <- p + geom_text(aes(x = gene_start, label = "Gene Start", y = chromosome),
                         vjust = -0.5, color = "green", size = 5)
    }
    
    # Return the plot
    p
  })
  # Add this inside your server function
  # output$extended_gene_region_samples <- renderDataTable({
  #   req(input$gene_name)  # Ensure a gene name has been inputted
  #   gene_name <- input$gene_name
  #   extended_samples <- find_introgressed_samples_by_extended_gene(gene_info, results_list_rename, gene_name)
  #   data.frame(Sample = extended_samples)
  # })
  
  data <- reactive({
    df <- read.csv("Bzea_merged.csv", stringsAsFactors = FALSE, na.strings = c("#N/A")) %>%
      filter(!is.na(Lines)) %>%
      mutate(Class = if_else(
        Lines == "B73",
        "B73",
        sub("\\..*", "", Lines)
      )) %>%
      mutate(Class = recode(
        Class,
        "Zh" = "Huehuetenango",
        "Zl" = "Luxurians",
        "Zx" = "Mexicana",
        "Zv" = "Parviglumis",
        "Zd" = "Diploperennis"
      ))
    df <- df %>%
      mutate(across(all_of(continuous_vars), ~ suppressWarnings(as.numeric(.))))
    
    df$DTS <- as.Date(df$DTS, format = "%m/%d/%y")
    df$DTA <- as.Date(df$DTA, format = "%m/%d/%y")
    
    updateSelectInput(session, "hist_var", choices = continuous_vars)
    updateSelectInput(session, "violin_var", choices = continuous_vars)
    updateSelectInput(session, "box_var", choices = continuous_vars)
    updateSelectInput(session, "density_var", choices = continuous_vars)
    updateSelectInput(session, "scatter_var1", choices = continuous_vars)
    updateSelectInput(session, "scatter_var2", choices = continuous_vars)
    updateSelectInput(session, "stat_var", choices = continuous_vars)
    
    df
  })
  
  renderFilter <- function(df, column) {
    if (column %in% names(df)) {
      col_data <- df[[column]]
      if (column %in% categorical_vars) {
        selectInput(
          inputId = column,
          label = paste("Select", column),
          choices = c("All" = "", sort(unique(col_data))),
          selected = ""
        )
      } else if (column %in% continuous_vars) {
        sliderInput(
          inputId = column,
          label = paste("Range of", column),
          min = min(col_data, na.rm = TRUE),
          max = max(col_data, na.rm = TRUE),
          value = c(min(col_data, na.rm = TRUE), max(col_data, na.rm = TRUE))
        )
      } else if (column %in% date_vars) {
        dateRangeInput(
          inputId = column,
          label = paste("Date range for", column),
          start = min(col_data, na.rm = TRUE),
          end = max(col_data, na.rm = TRUE)
        )
      }
    }
  }
  
  
  output$categorical_filters <- renderUI({
    df <- data()
    tagList(lapply(categorical_vars, function(col) renderFilter(df, col)))
  })
  
  output$numeric_date_filters <- renderUI({
    df <- data()
    tagList(
      lapply(continuous_vars, function(col) renderFilter(df, col)),
      lapply(date_vars, function(col) renderFilter(df, col))
    )
  })
  
  filteredData <- reactive({
    df <- data()
    for (column in c(categorical_vars, continuous_vars, date_vars)) {
      if (column %in% names(df) && !is.null(input[[column]])) {
        if (is.character(df[[column]]) || is.factor(df[[column]])) {
          if (input[[column]] != "") df <- df %>% filter(df[[column]] == input[[column]] | is.na(df[[column]]))
        } else if (is.numeric(df[[column]])) {
          df <- df %>% 
            filter(
              (df[[column]] >= input[[column]][1] & df[[column]] <= input[[column]][2]) |
                is.na(df[[column]]))
        } else if (inherits(df[[column]], "Date") || inherits(df[[column]], "POSIXct")) {
          df <- df %>% 
            filter(
              (df[[column]] >= input[[column]][1] & df[[column]] <= input[[column]][2]) |
                is.na(df[[column]]))
        }
      }
    }
    df
  })
  
  output$filteredTable <- DT::renderDataTable({
    DT::datatable(filteredData())
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(filteredData(), file, row.names = FALSE)
    }
  )
  
  observeEvent(input$reset_filters, {
    df <- data()
    for (var in categorical_vars) updateSelectInput(session, var, selected = "")
    for (var in continuous_vars) updateSliderInput(session, var, value = range(df[[var]], na.rm = TRUE))
    for (var in date_vars) updateDateRangeInput(session, var, start = min(df[[var]], na.rm = TRUE), end = max(df[[var]], na.rm = TRUE))
  })
  
  output$stats <- renderPrint({
    df <- data()
    if (!is.null(input$stat_var) & input$stat_var %in% names(df)) {
      variable <- df[[input$stat_var]]
      na_count <- sum(is.na(variable))
      list(
        Mean = mean(variable, na.rm = TRUE),
        Median = median(variable, na.rm = TRUE),
        SD = sd(variable, na.rm = TRUE),
        Min = min(variable, na.rm = TRUE),
        Max = max(variable, na.rm = TRUE),
        "Count of NA" = na_count,
        "Percentage of NA" = round(na_count / nrow(df) * 100, 2)
      )
    }
  })
  
  output$hist <- renderPlot({
    df <- data()
    if (!is.null(input$hist_var) & input$hist_var %in% names(df)) {
      hist(df[[input$hist_var]], main = paste(input$hist_var, "Histogram"), xlab = NULL, col = "skyblue")
    }
  })
  
  output$violin <- renderPlot({
    df <- data()
    if (!is.null(input$violin_var) & input$violin_var %in% names(df)) {
      ggplot(df, aes(x = "", y = .data[[input$violin_var]])) +
        geom_violin(fill = "skyblue") +
        labs(x = NULL, y = NULL, title = paste(input$violin_var, "Violin Plot")) +
        theme_minimal()
    }
  })
  
  output$boxplot <- renderPlot({
    df <- data()
    if (!is.null(input$box_var) & input$box_var %in% names(df)) {
      boxplot(df[[input$box_var]], main = paste(input$box_var, "Boxplot"), xlab = NULL, col = "skyblue")
    }
  })
  
  output$density <- renderPlot({
    df <- data()
    if (!is.null(input$density_var) & input$density_var %in% names(df)) {
      plot(density(df[[input$density_var]], na.rm = TRUE), main = paste(input$density_var, "Density Plot"))
    }
  })
  
  output$scatter <- renderPlot({
    df <- data()
    if (!is.null(input$scatter_var1) & !is.null(input$scatter_var2) &
        input$scatter_var1 %in% names(df) & input$scatter_var2 %in% names(df)) {
      ggplot(df, aes_string(x = input$scatter_var1, y = input$scatter_var2)) +
        geom_point() +
        labs(title = "Scatter Plot", x = input$scatter_var1, y = input$scatter_var2) +
        theme_minimal()
    }
  })
  
}


shinyApp(ui = ui, server = server)

