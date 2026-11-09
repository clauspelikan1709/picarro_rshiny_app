# R Shiny app to process and analyze raw Picarro CRDS analyzer data
# Modified version with corrected CH4 solubility units, water-vapour correction,
# and parallel 12C/13C isotope concentration calculations

# Packages
library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)


# Allow uploads larger than Shiny's 5 MB default (e.g. large Picarro .dat files)
options(shiny.maxRequestSize = 500 * 1024^2)


# Solubility and concentration conversion helpers
# Concentrations are given in ppm, volumes in mL, and temperature in K.

## solubility of CO2 in mol/L/atm at a given temperature and salinity
solubility_co2 <- function(temp_k, salin) {
    exp(-58.0931 + 90.5069 * (100 / temp_k) + 22.294 * log(temp_k / 100) +
            salin * (0.027766 - 0.025888 * (temp_k / 100) + 0.0050578 * (temp_k / 100)^2))
}

## solubility of CH4 as L gas (STP) per L water per atm at a given temperature and salinity
solubility_ch4 <- function(temp_k, salin) {
    exp(-68.8862 + 101.4956 * (100 / temp_k) + 28.7314 * log(temp_k / 100) +
            salin * (-0.076146 + 0.043970 * (temp_k / 100) - 0.0068672 * (temp_k / 100)^2))
}

## concentration of a gas in the headspace in micromole/L
headspace_uM <- function(ppm, temp_k, p_dry_atm) {
    (ppm * 1e-06 * p_dry_atm * 101325 / (8.314 * temp_k)) * 1000
}

## final dissolved CO2 concentration in the original water sample in micromole/L
dissolved_co2_uM <- function(ppm, temp_k, p_dry_atm, v_hs, v_liq, sol_co2) {
    hs_uM <- headspace_uM(ppm, temp_k, p_dry_atm)
    hs_tot <- hs_uM * v_hs / 1000
    aq_uM <- sol_co2 * (ppm * 1e-06 * p_dry_atm) * 1e6
    aq_tot <- aq_uM * v_liq / 1000
    ((hs_tot + aq_tot) / v_liq) * 1000
}

## final dissolved CH4 concentration in the original water sample in micromole/L
dissolved_ch4_uM <- function(ppm, temp_k, p_dry_atm, v_hs, v_liq, beta_ch4) {
    hs_uM <- headspace_uM(ppm, temp_k, p_dry_atm)
    hs_tot <- hs_uM * v_hs / 1000
    aq_uM <- beta_ch4 * (ppm * 1e-06 * p_dry_atm) / 22.414 * 1e6
    aq_tot <- aq_uM * v_liq / 1000
    ((hs_tot + aq_tot) / v_liq) * 1000
}


# Define UI
ui <- fluidPage( # begin ui
    # Application title
    titlePanel("Processing and analysis tool for Picarro CRDS analyzer data for headspace gas measurements"),
    
    
    # Split layout into side and main panel
    sidebarLayout( # begin sidebar layout
        
        # Sidebar for data import
        sidebarPanel( # begin sidebar panel
            width = 3,
            
            # Upload button for importing data files
            h3("Data import"),
            fileInput("dat_files", "Import data file(s)",
                      multiple = TRUE,
                      accept = c(".dat", "text/plain")),
            
            
            ## reactive dropdown menu to select individual file for analysis
            uiOutput("data_file"),
            
            
            ## enter sample name (optional)
            textInput("sample_name", "Enter sample name (optional)",
                      value = NA),
            
            # Reactive dropdown menus to select data columns for analysis
            h3("Select raw data columns"),
            
            p("12C columns must be dry mole fractions; 13C columns are measured wet and are converted to dry values in the app."),
            
            ## 12C and 13C CO2 columns
            uiOutput("co2_12_vals"),
            uiOutput("co2_13_vals"),
            
            ## 12C and 13C CH4 columns
            uiOutput("ch4_12_vals"),
            uiOutput("ch4_13_vals"),
            p(em("Use HR columns above 12 ppm CH4 and HP columns below 12 ppm.")),
            
            
            # Options for automatic plateau detection
            h3("Settings for automatic plateau detection"),
            
            ## threshold concentration for CO2
            numericInput("co2_thresh", "CO2 threshold concentration (ppm)",
                         value = 200, step = 10),
            
            ## variation tolerance
            numericInput("co2_tol", "Variation tolerance",
                         value = 0.001, step = 0.001),
            
            
            # Options for default measurement conditions
            h3("Default values for measurement conditions"),
            
            ## liquid sample volume
            numericInput("v_liq", "Liquid volume (mL)",
                         value = 60, step = 1),
            
            ## headspace volume
            numericInput("v_hs", "Headspace volume (mL)",
                         value = 60, step = 1),
            
            ## volume injected
            numericInput("v_inj", "Volume of injected gas sample (mL)",
                         value = 120, step = 1),
            
            ## salinity
            numericInput("sal", "Salinity (ppt)",
                         value = 0, step = 0.1),
            
            ## pressure
            numericInput("pressure", "Pressure (atm)",
                         value = 1, step = 0.1),


            # Options for table export
            h3("Table export options"),
            
            ## field separator
            selectInput("field_sep", "Field separator",
                        choices = c("semicolon (;)" = ";",
                                    "tab space (\\t)" = "\t"),
                        selected = "semicolon (;)"
            ),
            
            ## decimal sign
            selectInput("dec_sep", "Decimal separator",
                        choices = c("period (.)" = ".",
                                    "comma (,)" = ","),
                        selected = "period (.)"
            )
            
            
        ), # end sidebar panel
        
        
        # Main panel for plots
        mainPanel( # begin main panel
            
            # Tab layout
            tabsetPanel(type = "tabs", # begin tabs
                        
                        # Figure tab
                        tabPanel("Figures", # begin figure tab
                                 
                                 br(),
                                 
                                 # download buttons to save figures
                                 downloadButton("plots_pdf", label = "Export figure (.pdf)"),
                                 downloadButton("plots_png", label = "Export figure (.png)"),
                                 downloadButton("plots_svg", label = "Export figure (.svg)"),
                                 
                                 br(),
                                 
                                 # display figures
                                 plotOutput("plots", width = 800, height = 500),
                                 
                                 
                                 # slider widget to adjust time window for plateau
                                 uiOutput("plateau_window"),
                                 
                                 
                                 # calculation results for currently selected data file
                                 DTOutput("results_i"),
                                 
                                 br(),
                                 
                                 # action button to write results to collective table
                                 actionButton("add_results", "Submit data to Results")
                                 
                        ), # end figure tab
                        
                        
                        # Results tab
                        tabPanel("Results", # being results tab
                                 
                                 fluidRow(
                                     
                                     # table options and buttons 
                                     br(),
                                     column(3,     
                                            ## reactive dropdown menu for excluding rows
                                            uiOutput("excl_rows")
                                     ),
                                     
                                     column(4,
                                            ## download button to save table
                                            p(),
                                            downloadButton("table1_exp", "Export current table (.csv)")
                                     ),
                                     
                                     column(4, offset = 0.2,
                                            p(),
                                            ## action button to clear table
                                            actionButton("clear_results_tab", "Clear table")
                                     ),
                                     
                                 ), 
                                 
                                 # table with results for all data files
                                 DTOutput("results_all"),
                                 
                                 br(),
                                 
                                 # action button to submit samples to user inputs for calculating final concentrations
                                 actionButton("proceed", "Continue with selected samples")
                                 
                                 
                        ), # end results tab
                        
                        
                        # Measurement conditions tab for reactive user input
                        tabPanel("Measurement conditions", # begin measurement conditions tab
                                 
                                 # brief explanation
                                 br(),
                                 
                                 p(strong("Double-click on a field inside the table to enter or change values for:")),
                                 p("- Sample name (optional)"),
                                 p("- Sample temperature in", strong("°C.")),
                                 p("- Volumes of liquid, headspace, and the injected gas sample, respectively, in", strong("mL.")),
                                 p("- Salinity in", strong("ppt"), "(for freshwater samples, leave at 0)."),
                                 p("- Pressure in", strong("atm.")),
                                 
                                 br(),
                                 
                                 
                                 # buttons
                                 column(4,
                                        
                                        ## download button to save table
                                        downloadButton("table2_exp", "Export table (.csv)")
                                 ),
                                 
                                 column(4, offset = 1,
                                        
                                        ## action button to rest values in table to default
                                        actionButton("reset_tab", "Reset values to default")
                                 ),
                                 
                                 
                                 # editable table to enter measurement conditions
                                 DTOutput("meas_cond"),
                                 
                                 br(),
                                 
                                 # action button to pass data on for calculating final concentrations
                                 actionButton("calc_conc", "Calculate final concentrations"),
                                 
                                 
                        ), # end measurement conditions tab
                        
                        
                        
                        # Tab with calculated final concentrations of dissolved gases in original samples
                        tabPanel("Final concentrations",
                                 
                                 br(),   
                                 
                                 # download button to save table
                                 downloadButton("table3_exp", "Export table (.csv)"),
                                 
                                 br(),
                                 
                                 p("Final concentrations are given in µmol/L (µM) and nmol/L. delta13C and the atom fraction of 13C (AF_13C) are calculated from the final 12C and 13C concentrations."),
                                 
                                 # table with calculated results
                                 DTOutput("final_conc")
                                 
                        ),
                        
                        
                        # Contact info and rights
                        tabPanel(em("About"),
                                 
                                 br(),
                                 
                                 p("This tool was created by Lucas Fillinger in January 2022 at the Department of Functional and Evolutionary Ecology at the University of Vienna. The tool may be used, redistributed, and modified in any manner by anyone without restrictions."),
                                 
                                 br(),
                                 br(),
                                 p("For questions, remarks, or suggestions, contact:"),
                                 p(a("lucas.fillinger@univie.ac.at",
                                     href = "mailto:lucas.fillinger@univie.ac.at"), "or"),
                                 
                                 p(a("fillinger.microbecol@gmail.com",
                                     href = "mailto:fillinger.microbecol@gmail.com")),
                        )
                        
                        
            ) # end tabs
            , width = 9) # end main panel
    ) # end sidebar layout
    
) # end ui



#################################################
#################################################



# Define server logic
server <- function(input, output) { # begin server function
    
    # Raw data import and sample IDs
    ## extract sample IDs belonging to each data file
    sample_id <- reactive({
        
        sapply(strsplit(basename(input$dat_files$name), split = "-"), '[', 3)
        
    })
    
    
    ## import data files and store as individual dataframes in a list
    data_list <- reactive({
        
        lst <- list()
        for(i in 1:length(input$dat_files[, 1])){
            lst[[i]] <- as.data.frame(read.csv(input$dat_files[[i, "datapath"]], sep = ""))
        }
        
        names(lst) <- sample_id()
        return(lst)
        
    })
    
    
    
    # Interactive user inputs for analysis options
    ## render option for selecting an individual dataframe for the analysis
    output$data_file <- renderUI({
        
        if(is.null(input$dat_files)){return ()}
        
        selectInput("data_file", "Select file for analysis", choices = sample_id())
        
    })
    
    
    ## render options for selecting data columns for the analysis
    ### 12C and 13C CO2 concentrations
    output$co2_12_vals <- renderUI({
        
        if(is.null(input$dat_files)){return ()}
        
        selectInput("co2_12_vals", "Select column with 12C CO2 concentrations (dry)",
                    choices = c("12CO2_dry" = "X12CO2_dry"),
                    selected = "X12CO2_dry")
        
    })
    
    output$co2_13_vals <- renderUI({
        
        if(is.null(input$dat_files)){return ()}
        
        selectInput("co2_13_vals", "Select column with 13C CO2 concentrations (wet)",
                    choices = c("13CO2" = "X13CO2"),
                    selected = "X13CO2")
        
    })
    
    
    ### 12C and 13C CH4 concentrations
    output$ch4_12_vals <- renderUI({
        
        if(is.null(input$dat_files)){return ()}
        
        selectInput("ch4_12_vals", "Select column with 12C CH4 concentrations (dry)",
                    choices = c("HP_12CH4_dry" = "HP_12CH4_dry",
                                "HR_12CH4_dry" = "HR_12CH4_dry"),
                    selected = "HP_12CH4_dry")
        
    })
    
    output$ch4_13_vals <- renderUI({
        
        if(is.null(input$dat_files)){return ()}
        
        selectInput("ch4_13_vals", "Select column with 13C CH4 concentrations (wet)",
                    choices = c("HP_13CH4" = "HP_13CH4",
                                "HR_13CH4" = "HR_13CH4"),
                    selected = "HP_13CH4")
        
    })
    
    
    
    # Data analysis
    ## extract single dataframe from list to analyze, and:
    ### convert 'TIME' to POSIXct timestamp
    ### add additional column for measurement time in seconds starting at 0 for the first data point
    df_raw <- reactive({
        
        # wait until a data file has been selected
        req(input$data_file)
        
        data_list()[[input$data_file]] %>%
            mutate(datetime = as.POSIXct(paste(.$DATE, .$TIME, sep = " ")),
                   meas_time_sec = c(0, cumsum(as.numeric(diff(datetime)))))
        
    })
    
    
    ## define time window for plateau values
    ### plateau time points identified with default method based on variability of CO2 concentrations over time
    default_plateau <- reactive({
        
        req(input$co2_12_vals)
        
        # empty vector
        plat <- rep(0, nrow(df_raw()))
        
        # iterate over 10 consecutive time points;
        # for values at plateau (i.e. above a fixed defined threshold concentration and with a defined tolerance for standard deviation relative to mean values), record 1 in the vector, else 0
        
        for(i in 10:(length(plat))){ # begin loop
            
            if(df_raw()[i, input$co2_12_vals] > input$co2_thresh &&
               sd(df_raw()[i:(i-10), input$co2_12_vals]) / mean(df_raw()[i:(i-10), input$co2_12_vals]) < input$co2_tol) {plat[i] <- 1}
            else {plat[i] <- 0}
            
        } # end loop
        
        return(plat)
    })
    
    
    ### slider widget for adjusting plateau
    output$plateau_window <- renderUI({
        
        if(is.null(input$dat_files)){return ()}
        
        # use the full time range if no plateau was detected automatically
        plat_idx <- which(default_plateau() == 1)
        if(length(plat_idx) == 0){plat_idx <- seq_along(df_raw()$datetime)}
        
        sliderInput("plateau_window", "Adjust plateau time window",
                    min = min(df_raw()$datetime), max = max(df_raw()$datetime),
                    value = c(min(df_raw()[plat_idx, "datetime"]),
                              max(df_raw()[plat_idx, "datetime"])),
                    width = "550px"
                    
        )
        
    })
    
    
    ## create new dataframe for the analysis
    df_analysis <- reactive({
        
        # wait until the plateau time window and column selections have been defined
        req(input$plateau_window, input$co2_12_vals, input$co2_13_vals,
            input$ch4_12_vals, input$ch4_13_vals)
        
        df <- df_raw()
        
        co2_12_col <- input$co2_12_vals
        co2_13_col <- input$co2_13_vals
        ch4_12_col <- input$ch4_12_vals
        ch4_13_col <- input$ch4_13_vals
        
        # wet 12C counterpart of the selected dry 12C CO2 column
        co2_12_wet_col <- sub("_dry$", "", co2_12_col)
        
        # CH4 measurement mode (HP or HR) is defined by the selected 13C column;
        # the 12C dry/wet columns from the same mode are used for the dry correction
        ch4_mode <- sub("_13CH4$", "", ch4_13_col)
        ch4_12_wet_col <- paste0(ch4_mode, "_12CH4")
        ch4_12_dry_col <- paste0(ch4_mode, "_12CH4_dry")
        
        validate(
            need(sub("_12CH4_dry$", "", ch4_12_col) == ch4_mode,
                 "Please select 12C and 13C CH4 columns from the same measurement mode (HP or HR)."),
            need(all(c(co2_12_col, co2_13_col, co2_12_wet_col,
                       ch4_12_col, ch4_13_col, ch4_12_wet_col, ch4_12_dry_col) %in% names(df)),
                 "One or more selected columns are not present in this data file.")
        )
        
        data.frame(
            datetime = df$datetime,
            meas_time_sec = df$meas_time_sec,
            co2_12 = df[[co2_12_col]],
            co2_13 = df[[co2_13_col]] * (df[[co2_12_col]] / df[[co2_12_wet_col]]),
            ch4_12 = df[[ch4_12_col]],
            ch4_13 = df[[ch4_13_col]] * (df[[ch4_12_dry_col]] / df[[ch4_12_wet_col]]),
            plateau = df$datetime >= input$plateau_window[1] &
                      df$datetime <= input$plateau_window[2]
        )
        
    })
    
    
    ## funtion to plot data
    draw_plot <- function(text_large = 14, text_small = 12){
        
        df_analysis() %>%
            pivot_longer(cols = c(co2_12, co2_13, ch4_12, ch4_13),
                         values_to = "vals", names_to = "var") %>%
            mutate(var_lab = case_when(var == "co2_12" ~ "12CO2 (ppm, dry)",
                                       var == "co2_13" ~ "13CO2 (ppm, dry-corrected)",
                                       var == "ch4_12" ~ "12CH4 (ppm, dry)",
                                       var == "ch4_13" ~ "13CH4 (ppm, dry-corrected)"),
                   plateau = case_when(plateau == TRUE ~ "Yes",
                                       plateau == FALSE ~ "No")) %>%
            mutate(var_lab = fct_relevel(.$var_lab, c("12CO2 (ppm, dry)", "13CO2 (ppm, dry-corrected)", "12CH4 (ppm, dry)", "13CH4 (ppm, dry-corrected)")),
                   plateau = fct_relevel(.$plateau, c("Yes", "No"))) %>%
            ggplot(aes(x = datetime, y = vals, colour = plateau)) +
            geom_point(shape = 1, size = 3) +
            geom_vline(xintercept = input$plateau_window, colour = "red", linetype = "dashed", size = 1) +
            scale_colour_manual(values = c("red", "black")) +
            facet_wrap(.~ var_lab, scales = "free_y") +
            labs(x = "Time", y = NULL, colour = "Data points selected\nfor analysis",
                 title = paste("Data file ID:", input$data_file,
                               "; measured on:", unique(as.Date(df_analysis()$datetime)), sep = " ")) +
            theme_bw() +
            theme(panel.border = element_rect(colour = "black", size = 0.8),
                  panel.grid = element_blank(),
                  axis.title = element_text(colour = "black", size = text_large, face = "bold"),
                  axis.text = element_text(colour = "black", size = text_small),
                  axis.ticks = element_line(colour = "black"),
                  legend.title = element_text(colour = "black", size = text_large, face = "bold"),
                  legend.text = element_text(colour = "black", size = text_small),
                  strip.text = element_text(colour = "black", size = text_large, face = "bold"),
                  strip.background = element_rect(colour = "black", size = 0.8),
                  title = element_text(colour = "black", size = text_large, face = "bold"),
                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
        
    }
    
    
    ## render plot output
    output$plots <- renderPlot({
        
        if(is.null(input$dat_files)){return ()}
        draw_plot()
        
    })
    
    
    ## export figures
    ### pdf
    output$plots_pdf <- downloadHandler(
        
        filename = function(){
            paste(input$data_file, "-", gsub(pattern = "-", replacement = "", df_raw()[1, "DATE"]),
                  ".pdf", sep = "")
        },
        
        content = function(file){
            cairo_pdf(file, width = 8, height = 5)
            print(draw_plot(text_large = 12, text_small = 10))
            dev.off()
        }
    )
    
    
    ### png
    output$plots_png <- downloadHandler(
        
        filename = function(){
            paste(input$data_file, "-", gsub(pattern = "-", replacement = "", df_raw()[1, "DATE"]),
                  ".png", sep = "")
        },
        
        content = function(file){
            png(file, width = 8, height = 5, units = "in", res = 600)
            print(draw_plot(text_large = 12, text_small = 10))
            dev.off()
        }
    )
    
    
    ### svg
    output$plots_svg <- downloadHandler(
        
        filename = function(){
            paste(input$data_file, "-", gsub(pattern = "-", replacement = "", df_raw()[1, "DATE"]),
                  ".svg", sep = "")
        },
        
        content = function(file){
            svg(file, width = 8, height = 5)
            print(draw_plot(text_large = 12, text_small = 10))
            dev.off()
        }
    )
    
    ## calculations for data points at plateau for currently selected table
    ## calculated columns:
    ### time in seconds measured at plateau ('time_at_plateau_sec')
    ### number of data points at plateau ('n_data_points')
    ### mean values and standard deviations of concentrations and delta13C ratios for CO2 and CH4
    results_i_table <- reactive({
        
        df_analysis() %>% filter(plateau == TRUE) %>%
            select(-c(datetime, plateau)) %>%
            rename_all(~ c("sec", "co2_12", "co2_13", "ch4_12", "ch4_13")) %>%
            summarise(seconds_plateau = abs(diff(range(sec))),
                      n_data_points = n(),
                      mean_co2_12_ppm = mean(co2_12),
                      sd_co2_12_ppm = sd(co2_12),
                      mean_co2_13_ppm = mean(co2_13),
                      sd_co2_13_ppm = sd(co2_13),
                      mean_ch4_12_ppm = mean(ch4_12),
                      sd_ch4_12_ppm = sd(ch4_12),
                      mean_ch4_13_ppm = mean(ch4_13),
                      sd_ch4_13_ppm = sd(ch4_13),
                      mean_co2_total_ppm = mean(co2_12 + co2_13),
                      sd_co2_total_ppm = sd(co2_12 + co2_13),
                      mean_ch4_total_ppm = mean(ch4_12 + ch4_13),
                      sd_ch4_total_ppm = sd(ch4_12 + ch4_13))
        
    })
    
    
    ## render table with results for currently selected data file
    output$results_i <- renderDT({
        
        if(is.null(input$dat_files)){return ()}
        
        datatable(results_i_table(),
                  selection = "none",
                  rownames = FALSE,
                  options = list(paging = FALSE,
                                 ordering = FALSE,
                                 info = FALSE,
                                 bFilter = 0,
                                 scrollX = TRUE)
        )
    })
    
    
    
    # Save current results to collective results table
    ## reactive counter to record number of action button clicks
    click_count <- reactiveValues(value = 0)
    
    observeEvent(input$add_results, {
        click_count$value <- click_count$value + 1
    })
    
    
    ## create reactive values dataframe to store results
    results_collect <- reactiveValues(df = NULL)
    
    
    ## upon action button click `add_results`, write results from individual sample to the table
    observeEvent(input$add_results, {
        
        # create new table on first click of action button 
        if(click_count$value == 1){
            results_collect$df <- results_i_table()  %>%
                mutate(data_file = input$data_file,
                       date = df_raw()[1, "DATE"],
                       sample_name = input$sample_name,
                       .before = 1)
        }
        
        # if a table already exists after the first click, add new results to existing table
        else{
            results_collect$df <- rbind(results_collect$df,
                                         results_i_table()  %>%
                                             mutate(data_file = input$data_file,
                                                    date = df_raw()[1, "DATE"],
                                                    sample_name = input$sample_name,
                                                    .before = 1))
        }
        
    })
    
    
    ## render choices to exclude individual rows in results table
    output$excl_rows <- renderUI({
        
        if(is.null(input$dat_files) || is.null(results_collect$df)){return ()}
        
        selectInput("excl_rows", "Exclude table rows",
                    choices = c("None" = 0,
                                1:nrow(results_collect$df)),
                    selected = 0,
                    multiple = TRUE)
        
    })
    
    
    ## subset results table based on selected rows to exclude
    results_all_final <- reactive({
        
        if(is.null(results_collect$df)){return (NULL)}
        
        excl <- as.numeric(input$excl_rows)
        excl <- excl[!is.na(excl) & excl != 0]
        
        if(length(excl) == 0){
            results_collect$df
        }
        else{
            results_collect$df[-excl, , drop = FALSE]
        }
        
    })
    
    
    ## render results table
    output$results_all <- renderDT({
        
        if(is.null(results_collect$df)){return ()}
        
        datatable(results_all_final() %>%
                      mutate(across(where(is.numeric), ~ round(.x, digits = 4))),
                  
                  options = list(paging = FALSE,
                                 ordering = FALSE,
                                 info = FALSE,
                                 bFilter = 0,
                                 scrollX = TRUE,
                                 scrollY = "450px"
                  )
        )
        
    })
    
    
    ## upon action button click `clear_results_tab`, set `results_collect$df` back to NULL and click counter back to 0
    observeEvent(input$clear_results_tab, {
        
        results_collect$df <- NULL
        click_count$value <- 0
        
    })
    
    
    ## export results table
    output$table1_exp <- downloadHandler(
        filename = paste("table_A_", format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss"), ".csv", sep = ""),
        
        content = function(file){
            write.table(results_all_final(), file,
                        col.names = TRUE,
                        row.names = FALSE,
                        na = "",
                        sep = input$field_sep,
                        dec = input$dec_sep
            )
        }
    )
    
    
    
    # Input of measurement conditions for subsequent calculation of final concentrations of dissolved gases in the original sample
    
    ## create reactive values dataframe to store user inputs 
    meas_cond_user <- reactiveValues(df = NULL)
    
    
    ## upon click on action button `proceed`, write current selection of samples from `results_all_final` to `meas_cond_user$df`
    ## add new columns for user inputs for:
    ### sample temperature in °C
    ### volumes of liquid, headspace, and injected sample gas in mL
    ### salinity (set to default of 0)
    ### pressure in atm (set to default of 1)
    observeEvent(input$proceed, {
        
        meas_cond_user$df <- isolate(results_all_final()) %>%
            select(c("data_file", "date", "sample_name")) %>%
            mutate(temperature_C = as.numeric(NA),
                   liquid_mL = input$v_liq,
                   headspace_mL = input$v_hs,
                   injected_mL = input$v_inj,
                   salinity_ppt = input$sal,
                   pressure_atm = input$pressure)
        
    })
    
    
    ## render editable table from `meas_cond_user$df`
    output$meas_cond <- renderDT({
        
        if(is.null(meas_cond_user$df)){return ()}
        
        datatable(meas_cond_user$df,
                  selection = "none",
                  editable = list(target = "cell", disable = list(columns = c(1:2))),
                  rownames = TRUE,
                  options = list(paging = FALSE,
                                 ordering = FALSE,
                                 info = FALSE,
                                 bFilter = 0,
                                 scrollX = TRUE,
                                 scrollY = "450px")
                  
        )
        
    })
    
    
    ## create proxy table from user inputs 
    meas_cond_proxy <- dataTableProxy("meas_cond")
    
    
    ## update `meas_cond_user$df` with user inputs
    observeEvent(input$meas_cond_cell_edit, {
        
        info = input$meas_cond_cell_edit
        str(info)
        row_ind = info$row
        col_ind = info$col
        new_val = info$value
        
        meas_cond_user$df[row_ind, col_ind] <- coerceValue(new_val, meas_cond_user$df[row_ind, col_ind])
        
    })
    
    
    ## upon action button click `reset_tab`, reset values back to original
    observeEvent(input$reset_tab, {
        
        meas_cond_user$df <- isolate(results_all_final()) %>%
            select(c("data_file", "date", "sample_name")) %>%
            mutate(temperature_C = as.numeric(NA),
                   liquid_mL = input$v_liq,
                   headspace_mL = input$v_hs,
                   injected_mL = input$v_inj,
                   salinity_ppt = input$sal,
                   pressure_atm = input$pressure)
        
    })
    
    
    
    ## download table
    output$table2_exp <- downloadHandler(
        
        filename = paste("table_B_", format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss"), ".csv", sep = ""),
        
        content = function(file){
            write.table(meas_cond_user$df, file,
                        col.names = TRUE,
                        row.names = FALSE,
                        na = "",
                        sep = input$field_sep,
                        dec = input$dec_sep
            )
        }
        
    ) 
    
    
    
    # Calculation of final concentrations of dissolved gases in the original samples
    
    ## upon action button click `calc_conc`:
    ## print out warning message before calculating final concentrations if temperature data is missing
    observeEvent(input$calc_conc, {
        
        if(sum(is.na(meas_cond_user$df$temperature_C)) != 0){
            
            showNotification("Warning: Make sure to enter temperature values for all samples before proceeding to the calculations!",
                             type = "error",
                             duration = NULL)
            
        }
        
    })
    
    
    ## combine data of measured values from `results_all_final()` with user input data from `meas_con_user$df`
    df_for_calc <- reactiveValues(df = NULL)
    
    observeEvent(input$calc_conc, {
        
        df_for_calc$df <- isolate(results_all_final()) %>%
            select(c("data_file", "date",
                     "mean_co2_12_ppm", "sd_co2_12_ppm",
                     "mean_co2_13_ppm", "sd_co2_13_ppm",
                     "mean_co2_total_ppm", "sd_co2_total_ppm",
                     "mean_ch4_12_ppm", "sd_ch4_12_ppm",
                     "mean_ch4_13_ppm", "sd_ch4_13_ppm",
                     "mean_ch4_total_ppm", "sd_ch4_total_ppm")) %>%
            mutate(sample_name = meas_cond_user$df$sample_name, .before = 3) %>%
            rename(co2_12_meas = mean_co2_12_ppm,
                   sd_co2_12 = sd_co2_12_ppm,
                   co2_13_meas = mean_co2_13_ppm,
                   sd_co2_13 = sd_co2_13_ppm,
                   co2_total_meas = mean_co2_total_ppm,
                   sd_co2_total = sd_co2_total_ppm,
                   ch4_12_meas = mean_ch4_12_ppm,
                   sd_ch4_12 = sd_ch4_12_ppm,
                   ch4_13_meas = mean_ch4_13_ppm,
                   sd_ch4_13 = sd_ch4_13_ppm,
                   ch4_total_meas = mean_ch4_total_ppm,
                   sd_ch4_total = sd_ch4_total_ppm) %>%
            mutate(temp_c = meas_cond_user$df$temperature_C,
                   v_liq = meas_cond_user$df$liquid_mL,
                   v_hs = meas_cond_user$df$headspace_mL,
                   v_inj = meas_cond_user$df$injected_mL,
                   salin = meas_cond_user$df$salinity_ppt,
                   pressure_atm = meas_cond_user$df$pressure_atm)
        
    })
    
    
    ## dataframe with calculations of final concentrations
    calc_conc_raw <- reactive({
        
        df_for_calc$df %>%
            mutate(
                # temperature and water vapour conditions
                temp_k = temp_c + 273.15, # temperature in kelvin
                p_h2o = exp(24.4543 - 67.4509 * (100 / temp_k) - 4.8489 * log(temp_k / 100) - 0.000544 * salin), # water vapour partial pressure in atm
                p_dry = pressure_atm - p_h2o, # 12C and dry-corrected 13C values are dry mole fractions
                
                # solubility at given temperature and salinity
                sol_co2 = solubility_co2(temp_k, salin), # solubility of co2 in mol/L/atm
                beta_ch4 = solubility_ch4(temp_k, salin), # solubility of ch4 as L gas (STP) per L water per atm
                
                # CO2 concentrations (12C, 13C and total)
                co2_12_final_uM = dissolved_co2_uM(co2_12_meas, temp_k, p_dry, v_hs, v_liq, sol_co2),
                co2_13_final_uM = dissolved_co2_uM(co2_13_meas, temp_k, p_dry, v_hs, v_liq, sol_co2),
                co2_total_final_uM = co2_12_final_uM + co2_13_final_uM,
                
                # CH4 concentrations (12C, 13C and total)
                ch4_12_final_uM = dissolved_ch4_uM(ch4_12_meas, temp_k, p_dry, v_hs, v_liq, beta_ch4),
                ch4_13_final_uM = dissolved_ch4_uM(ch4_13_meas, temp_k, p_dry, v_hs, v_liq, beta_ch4),
                ch4_total_final_uM = ch4_12_final_uM + ch4_13_final_uM,
                
                # technical standard deviations for 12C and 13C
                sd_co2_12_final_uM = dissolved_co2_uM(sd_co2_12, temp_k, p_dry, v_hs, v_liq, sol_co2),
                sd_co2_13_final_uM = dissolved_co2_uM(sd_co2_13, temp_k, p_dry, v_hs, v_liq, sol_co2),
                sd_ch4_12_final_uM = dissolved_ch4_uM(sd_ch4_12, temp_k, p_dry, v_hs, v_liq, beta_ch4),
                sd_ch4_13_final_uM = dissolved_ch4_uM(sd_ch4_13, temp_k, p_dry, v_hs, v_liq, beta_ch4),
                sd_co2_total_final_uM = dissolved_co2_uM(sd_co2_total, temp_k, p_dry, v_hs, v_liq, sol_co2),
                sd_ch4_total_final_uM = dissolved_ch4_uM(sd_ch4_total, temp_k, p_dry, v_hs, v_liq, beta_ch4),
                
                # concentrations in nmol/L (1 uM = 1000 nmol/L)
                co2_12_final_nmol_L = co2_12_final_uM * 1000,
                co2_13_final_nmol_L = co2_13_final_uM * 1000,
                co2_total_final_nmol_L = co2_total_final_uM * 1000,
                ch4_12_final_nmol_L = ch4_12_final_uM * 1000,
                ch4_13_final_nmol_L = ch4_13_final_uM * 1000,
                ch4_total_final_nmol_L = ch4_total_final_uM * 1000,
                
                # isotope ratios, delta13C and atom fraction from micromolar concentrations
                R_ch4 = ch4_13_final_uM / ch4_12_final_uM,
                delta13C_ch4_permil = (R_ch4 / 0.0111802 - 1) * 1000,
                AF_13C_ch4 = ch4_13_final_uM / ch4_total_final_uM,
                
                R_co2 = co2_13_final_uM / co2_12_final_uM,
                delta13C_co2_permil = (R_co2 / 0.0111802 - 1) * 1000,
                AF_13C_co2 = co2_13_final_uM / co2_total_final_uM
            )
        
    })
    
    
    ## extract relevant columns from raw calculations table generated above
    final_conc_tab <- reactive({
        
        calc_conc_raw() %>%
            select(c("data_file", "date", "sample_name",
                     "co2_12_final_uM", "sd_co2_12_final_uM", "co2_12_final_nmol_L",
                     "co2_13_final_uM", "sd_co2_13_final_uM", "co2_13_final_nmol_L",
                     "co2_total_final_uM", "sd_co2_total_final_uM", "co2_total_final_nmol_L",
                     "ch4_12_final_uM", "sd_ch4_12_final_uM", "ch4_12_final_nmol_L",
                     "ch4_13_final_uM", "sd_ch4_13_final_uM", "ch4_13_final_nmol_L",
                     "ch4_total_final_uM", "sd_ch4_total_final_uM", "ch4_total_final_nmol_L",
                     "R_ch4", "delta13C_ch4_permil", "AF_13C_ch4",
                     "R_co2", "delta13C_co2_permil", "AF_13C_co2"))
        
    })
    
    
    ## render table with final concentrations
    output$final_conc <- renderDT({
        
        if(is.null(df_for_calc$df)){return ()}
        
        datatable(final_conc_tab(),
                  selection = "none",
                  rownames = TRUE,
                  options = list(paging = FALSE,
                                 ordering = FALSE,
                                 info = FALSE,
                                 bFilter = 0,
                                 scrollX = TRUE,
                                 scrollY = "450px")
                  
        )
        
    })
    
    
    ## download table
    output$table3_exp <- downloadHandler(
        
        filename = paste("table_C_", format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss"), ".csv", sep = ""),
        
        content = function(file){
            write.table(final_conc_tab(), file,
                        col.names = TRUE,
                        row.names = FALSE,
                        na = "",
                        sep = input$field_sep,
                        dec = input$dec_sep
            )
        }
    )
    
    
} # end server function


# Run the application 
shinyApp(ui = ui, server = server)
