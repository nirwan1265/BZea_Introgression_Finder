# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)

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
  dashboardHeader(title = "Bzea FieldBook Analysis"),
  dashboardSidebar(
    sidebarMenu(
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

shinyApp(ui, server)
