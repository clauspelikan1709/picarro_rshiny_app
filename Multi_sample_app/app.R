# R Shiny app to process and analyze Picarro CRDS autosampler output
# Adapted from the single-file app (Single_file_app/app.R) for the
# "..._Results.csv" files produced by the autosampler software.
#
# In these files, plateau detection has already been performed and each row
# contains the measured mean and standard deviation for one sample (vial).
# The app therefore skips plateau detection and starts from those means/SDs,
# but keeps the same water-vapour correction, 12C/13C dry correction,
# solubility calculation and final dissolved-gas output as the single app.

# Packages
library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)


# Allow uploads larger than Shiny's 5 MB default
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


# Helpers for locating columns in the autosampler results file
## build the name of a result block column, e.g. "HP_12CH4_dry measured mean"
result_col <- function(df, prefix, stat = "measured mean") {
    nm <- trimws(names(df))
    target <- paste(prefix, stat)
    hit <- nm[tolower(nm) == tolower(target)]
    if (length(hit) == 0) {
        pattern <- paste0("^", gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", prefix),
                          "\\s+", stat, "$")
        hit <- grep(pattern, nm, ignore.case = TRUE, value = TRUE)
    }
    if (length(hit) == 0) NA_character_ else hit[1]
}

## name of the standard-deviation column matching a measured-mean column
sd_from_mean <- function(col) {
    sub("measured mean", "standard deviation", col, fixed = TRUE)
}

## return a column from a data frame, or NA if it is absent
col_or_na <- function(df, nm) {
    if (nm %in% names(df)) df[[nm]] else rep(NA, nrow(df))
}


# Define UI
ui <- fluidPage( # begin ui
    # Application title
    titlePanel("Processing and analysis tool for Picarro CRDS autosampler data"),
    
    
    # Split layout into side and main panel
    sidebarLayout( # begin sidebar layout
        
        # Sidebar for data import
        sidebarPanel( # begin sidebar panel
            width = 3,
            
            # Upload button for importing data files
            h3("Data import"),
            fileInput("csv_file", "Import autosampler results file (.csv)",
                      multiple = FALSE,
                      accept = c(".csv", "text/plain")),
            
            p(em("Upload the '*_Results.csv' file produced by the autosampler software. Plateau detection is not needed because the file already contains the measured mean and standard deviation for each sample.")),
            
            
            # Reactive dropdown menus to select data columns for analysis
            h3("Select result columns"),
            
            p("12C columns must be dry mole fractions; 13C columns are measured wet and are converted to dry values in the app."),
            
            ## 12C and 13C CO2 columns
            uiOutput("co2_12_vals"),
            uiOutput("co2_12_wet_vals"),
            uiOutput("co2_13_vals"),
            
            ## CH4 column set
            h4("CH4 measurement mode"),
            selectInput("ch4_mode", "Select CH4 columns",
                        choices = c("Automatic (HP below, HR above threshold)" = "auto",
                                    "HP only" = "HP",
                                    "HR only" = "HR"),
                        selected = "auto"),
            numericInput("ch4_switch", "CH4 threshold for automatic HP/HR selection (ppm)",
                         value = 12, min = 0, step = 0.1),
            p(em("Use HR columns above 12 ppm CH4 and HP columns below 12 ppm. In automatic mode the choice is made separately for every sample row.")),
            
            
            # Options for default measurement conditions
            h3("Default values for measurement conditions"),
            
            ## sample temperature
            numericInput("temp", "Temperature (°C)",
                         value = 20, step = 0.1),
            
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
        
        
        # Main panel for plots and tables
        mainPanel( # begin main panel
            
            # Tab layout
            tabsetPanel(type = "tabs", # begin tabs
                        
                        # Figure tab with a summary of the imported samples
                        tabPanel("Figures", # begin figure tab
                                 
                                 br(),
                                 
                                 # download buttons to save figures
                                 downloadButton("plots_pdf", label = "Export figure (.pdf)"),
                                 downloadButton("plots_png", label = "Export figure (.png)"),
                                 downloadButton("plots_svg", label = "Export figure (.svg)"),
                                 
                                 br(),
                                 
                                 # display figures
                                 plotOutput("plots", width = 800, height = 600)
                                 
                        ), # end figure tab
                        
                        
                        # Results tab
                        tabPanel("Results", # begin results tab
                                 
                                 fluidRow(
                                     
                                     br(),
                                     
                                     # table options and buttons
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
                                            ## action button to reset excluded rows
                                            actionButton("clear_excl", "Clear exclusions")
                                     ),
                                     
                                 ),
                                 
                                 p("The table shows the measured mean and standard deviation of each sample row. 13C values are already dry-corrected. Use 'Exclude table rows' to remove samples from the subsequent calculations."),
                                 
                                 # table with results for all sample rows
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
                                        
                                        ## action button to reset values in table to default
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
                                 
                                 p("This autosampler version adapts the same calculations to the '*_Results.csv' output of the Picarro autosampler software. Plateau detection is skipped because the input file already contains per-sample means and standard deviations."),
                                 
                                 p("Because only means and standard deviations are available in the autosampler file, the dry correction uses the ratio of the reported 12C dry and wet means, and the standard deviation of a dry-corrected 13C value is approximated by scaling its wet standard deviation by that same ratio. Total standard deviations are combined as the square root of the sum of the squared 12C and 13C standard deviations."),
                                 
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
server <- function(input, output, session) { # begin server function
    
    # Raw data import
    ## import the autosampler results file and remove a possible UTF-8 BOM from the first column name
    dat_raw <- reactive({
        
        req(input$csv_file)
        
        df <- read.csv(input$csv_file$datapath,
                       check.names = FALSE,
                       stringsAsFactors = FALSE,
                       fileEncoding = "UTF-8-BOM")
        
        names(df) <- trimws(sub("^\ufeff", "", names(df)))
        
        validate(
            need(nrow(df) > 0, "The selected file does not contain any data rows.")
        )
        
        df
        
    })
    
    
    # Interactive user inputs for analysis options
    ## render options for selecting result columns (measured means)
    output$co2_12_vals <- renderUI({
        
        req(dat_raw())
        nm <- names(dat_raw())
        ch <- nm[grepl("12CO2_dry", nm, ignore.case = TRUE) &
                    grepl("measured mean$", nm, ignore.case = TRUE)]
        
        if (length(ch) == 0) {
            return(helpText("No '12CO2_dry ... measured mean' column was found in this file."))
        }
        
        selectInput("co2_12_vals", "Select column with 12C CO2 concentrations (dry)",
                    choices = ch,
                    selected = ch[1])
        
    })
    
    output$co2_12_wet_vals <- renderUI({
        
        req(dat_raw())
        nm <- names(dat_raw())
        ch <- nm[grepl("12CO2", nm, ignore.case = TRUE) &
                    !grepl("_dry", nm, ignore.case = TRUE) &
                    grepl("measured mean$", nm, ignore.case = TRUE)]
        
        if (length(ch) == 0) {
            return(helpText("No '12CO2 ... measured mean' (wet) column was found in this file."))
        }
        
        selectInput("co2_12_wet_vals", "Select column with 12C CO2 concentrations (wet, for dry correction)",
                    choices = ch,
                    selected = ch[1])
        
    })
    
    output$co2_13_vals <- renderUI({
        
        req(dat_raw())
        nm <- names(dat_raw())
        ch <- nm[grepl("13CO2", nm, ignore.case = TRUE) &
                    grepl("measured mean$", nm, ignore.case = TRUE)]
        
        if (length(ch) == 0) {
            return(helpText("No '13CO2 ... measured mean' column was found in this file."))
        }
        
        selectInput("co2_13_vals", "Select column with 13C CO2 concentrations (wet)",
                    choices = ch,
                    selected = ch[1])
        
    })
    
    
    # Data analysis
    ## extract the needed columns and apply the same dry correction as the single-file app:
    ### 12C CO2 is already dry
    ### 13C CO2 (wet) is multiplied by the 12C dry/wet ratio
    ### 12C CH4 is already dry (HP or HR set)
    ### 13C CH4 (wet) is multiplied by the 12C dry/wet ratio of the same HP/HR set
    results_raw <- reactive({
        
        df <- dat_raw()
        
        req(input$co2_12_vals, input$co2_12_wet_vals, input$co2_13_vals)
        
        # required CO2 columns and their SD counterparts
        co2_12_col <- input$co2_12_vals
        co2_12_wet_col <- input$co2_12_wet_vals
        co2_13_col <- input$co2_13_vals
        co2_12_sd_col <- sd_from_mean(co2_12_col)
        co2_13_sd_col <- sd_from_mean(co2_13_col)
        
        # required CH4 columns for both modes
        hp_dry_col <- result_col(df, "HP_12CH4_dry")
        hp_wet_col <- result_col(df, "HP_12CH4")
        hp_13_col <- result_col(df, "HP_13CH4")
        hr_dry_col <- result_col(df, "HR_12CH4_dry")
        hr_wet_col <- result_col(df, "HR_12CH4")
        hr_13_col <- result_col(df, "HR_13CH4")
        
        hp_available <- !anyNA(c(hp_dry_col, hp_wet_col, hp_13_col)) &&
            all(c(hp_dry_col, hp_wet_col, hp_13_col) %in% names(df))
        hr_available <- !anyNA(c(hr_dry_col, hr_wet_col, hr_13_col)) &&
            all(c(hr_dry_col, hr_wet_col, hr_13_col) %in% names(df))
        
        validate(
            need(all(c(co2_12_col, co2_12_wet_col, co2_13_col) %in% names(df)),
                 "One or more selected CO2 columns are not present in this data file."),
            need(hp_available || hr_available,
                 "Could not find a complete HP or HR CH4 column set (12CH4_dry, 12CH4, 13CH4 measured means) in this file.")
        )
        
        num <- function(col) {
            if (is.na(col) || !col %in% names(df)) return(rep(NA_real_, nrow(df)))
            suppressWarnings(as.numeric(df[[col]]))
        }
        
        # choose the CH4 column set for every row
        if (input$ch4_mode == "HP") {
            mode <- rep("HP", nrow(df))
        } else if (input$ch4_mode == "HR") {
            mode <- rep("HR", nrow(df))
        } else {
            # automatic: use HP below the threshold and HR above it
            basis <- if (hp_available) num(hp_dry_col) else num(hr_dry_col)
            mode <- ifelse(!is.na(basis) & basis < input$ch4_switch, "HP", "HR")
            # fall back to the only available set
            if (!hp_available) mode[] <- "HR"
            if (!hr_available) mode[] <- "HP"
        }
        
        if (!hp_available && any(mode == "HP")) {
            mode[mode == "HP"] <- "HR"
        }
        if (!hr_available && any(mode == "HR")) {
            mode[mode == "HR"] <- "HP"
        }
        
        pick <- function(hp_col, hr_col) {
            ifelse(mode == "HP", num(hp_col), num(hr_col))
        }
        
        # CO2 dry correction
        co2_12 <- num(co2_12_col)
        co2_12_wet <- num(co2_12_wet_col)
        co2_13_wet <- num(co2_13_col)
        co2_factor <- co2_12 / co2_12_wet
        
        co2_12_sd <- num(co2_12_sd_col)
        co2_13_sd <- num(co2_13_sd_col) * co2_factor
        
        # CH4 dry correction within the selected HP/HR set
        ch4_12 <- pick(hp_dry_col, hr_dry_col)
        ch4_12_wet <- pick(hp_wet_col, hr_wet_col)
        ch4_13_wet <- pick(hp_13_col, hr_13_col)
        ch4_factor <- ch4_12 / ch4_12_wet
        
        ch4_12_sd <- pick(sd_from_mean(hp_dry_col), sd_from_mean(hr_dry_col))
        ch4_13_sd <- pick(sd_from_mean(hp_13_col), sd_from_mean(hr_13_col)) * ch4_factor
        
        # 13C dry-corrected values
        co2_13 <- co2_13_wet * co2_factor
        ch4_13 <- ch4_13_wet * ch4_factor
        
        data.frame(
            group_name = col_or_na(df, "Group Name"),
            vial = col_or_na(df, "Vial"),
            sample_type = col_or_na(df, "Sample Type"),
            standard_name = col_or_na(df, "Standard Name"),
            sampling_method = col_or_na(df, "Sampling Method"),
            pressure_flagged = col_or_na(df, "pressure_flagged"),
            ch4_mode = mode,
            mean_co2_12_ppm = co2_12,
            sd_co2_12_ppm = co2_12_sd,
            mean_co2_13_ppm = co2_13,
            sd_co2_13_ppm = co2_13_sd,
            mean_co2_total_ppm = co2_12 + co2_13,
            sd_co2_total_ppm = sqrt(co2_12_sd^2 + co2_13_sd^2),
            mean_ch4_12_ppm = ch4_12,
            sd_ch4_12_ppm = ch4_12_sd,
            mean_ch4_13_ppm = ch4_13,
            sd_ch4_13_ppm = ch4_13_sd,
            mean_ch4_total_ppm = ch4_12 + ch4_13,
            sd_ch4_total_ppm = sqrt(ch4_12_sd^2 + ch4_13_sd^2),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        
    })
    
    
    ## render choices to exclude individual rows in results table
    output$excl_rows <- renderUI({
        
        df <- results_raw()
        if (is.null(df)){return ()}
        
        row_choices <- seq_len(nrow(df))
        names(row_choices) <- as.character(row_choices)
        
        selectInput("excl_rows", "Exclude table rows",
                    choices = c("None" = 0, row_choices),
                    selected = 0,
                    multiple = TRUE)
        
    })
    
    
    ## subset results table based on selected rows to exclude
    results_all_final <- reactive({
        
        df <- results_raw()
        if (is.null(df)){return (NULL)}
        
        excl <- as.numeric(input$excl_rows)
        excl <- excl[!is.na(excl) & excl != 0]
        
        if(length(excl) == 0){
            df
        }
        else{
            df[-excl, , drop = FALSE]
        }
        
    })
    
    
    ## reset excluded rows
    observeEvent(input$clear_excl, {
        updateSelectInput(session, "excl_rows", selected = 0)
    })
    
    
    ## render results table
    output$results_all <- renderDT({
        
        if(is.null(results_all_final())){return ()}
        
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
    
    
    ## export results table
    output$table1_exp <- downloadHandler(
        filename = paste("table_A_", format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss"), ".csv", sep = ""),
        
        content = function(file){
            df <- results_all_final()
            if (!is.null(df)) {
                write.table(df, file,
                            col.names = TRUE,
                            row.names = FALSE,
                            na = "",
                            sep = input$field_sep,
                            dec = input$dec_sep
                )
            }
        }
    )
    
    
    ## plot summary of the imported samples
    draw_plot <- function(text_large = 14, text_small = 12){
        
        df <- results_all_final()
        req(df, nrow(df) > 0)
        
        df %>%
            mutate(sample_label = paste0(group_name,
                                         ifelse(is.na(vial), "", paste0(" (vial ", vial, ")")))) %>%
            pivot_longer(cols = c(mean_co2_12_ppm, mean_co2_13_ppm,
                                  mean_ch4_12_ppm, mean_ch4_13_ppm),
                         values_to = "vals", names_to = "var") %>%
            mutate(var_lab = case_when(var == "mean_co2_12_ppm" ~ "12CO2 (ppm, dry)",
                                       var == "mean_co2_13_ppm" ~ "13CO2 (ppm, dry-corrected)",
                                       var == "mean_ch4_12_ppm" ~ "12CH4 (ppm, dry)",
                                       var == "mean_ch4_13_ppm" ~ "13CH4 (ppm, dry-corrected)"),
                   var_lab = factor(var_lab,
                                    levels = c("12CO2 (ppm, dry)", "13CO2 (ppm, dry-corrected)",
                                               "12CH4 (ppm, dry)", "13CH4 (ppm, dry-corrected)"))) %>%
            ggplot(aes(x = sample_label, y = vals)) +
            geom_point(shape = 1, size = 3) +
            facet_wrap(.~ var_lab, scales = "free_y", ncol = 2) +
            labs(x = NULL, y = "Mole fraction (ppm)",
                 title = "Measured concentrations per sample") +
            theme_bw() +
            theme(panel.border = element_rect(colour = "black", size = 0.8),
                  panel.grid = element_blank(),
                  axis.title = element_text(colour = "black", size = text_large, face = "bold"),
                  axis.text = element_text(colour = "black", size = text_small),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
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
        
        if(is.null(results_all_final())){return ()}
        draw_plot()
        
    })
    
    
    ## export figures
    ### pdf
    output$plots_pdf <- downloadHandler(
        
        filename = function(){
            paste("autosampler_overview-", format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss"), ".pdf", sep = "")
        },
        
        content = function(file){
            cairo_pdf(file, width = 8, height = 6)
            print(draw_plot(text_large = 12, text_small = 10))
            dev.off()
        }
    )
    
    
    ### png
    output$plots_png <- downloadHandler(
        
        filename = function(){
            paste("autosampler_overview-", format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss"), ".png", sep = "")
        },
        
        content = function(file){
            png(file, width = 8, height = 6, units = "in", res = 600)
            print(draw_plot(text_large = 12, text_small = 10))
            dev.off()
        }
    )
    
    
    ### svg
    output$plots_svg <- downloadHandler(
        
        filename = function(){
            paste("autosampler_overview-", format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss"), ".svg", sep = "")
        },
        
        content = function(file){
            svg(file, width = 8, height = 6)
            print(draw_plot(text_large = 12, text_small = 10))
            dev.off()
        }
    )
    
    
    # Input of measurement conditions for subsequent calculation of final concentrations of dissolved gases in the original sample
    
    ## create reactive values dataframe to store user inputs 
    meas_cond_user <- reactiveValues(df = NULL)
    
    
    ## build the editable measurement-conditions table from the selected results
    default_meas_cond <- function(){
        isolate(results_all_final()) %>%
            select(c("group_name", "vial", "sample_type", "standard_name")) %>%
            mutate(sample_name = group_name,
                   temperature_C = input$temp,
                   liquid_mL = input$v_liq,
                   headspace_mL = input$v_hs,
                   injected_mL = input$v_inj,
                   salinity_ppt = input$sal,
                   pressure_atm = input$pressure)
    }
    
    
    ## upon click on action button `proceed`, write current selection of samples from `results_all_final` to `meas_cond_user$df`
    observeEvent(input$proceed, {
        
        req(results_all_final())
        
        meas_cond_user$df <- default_meas_cond()
        
    })
    
    
    ## render editable table from `meas_cond_user$df`
    output$meas_cond <- renderDT({
        
        if(is.null(meas_cond_user$df)){return ()}
        
        datatable(meas_cond_user$df,
                  selection = "none",
                  editable = list(target = "cell", disable = list(columns = c(1:4))),
                  rownames = TRUE,
                  options = list(paging = FALSE,
                                 ordering = FALSE,
                                 info = FALSE,
                                 bFilter = 0,
                                 scrollX = TRUE,
                                 scrollY = "450px")
                  
        )
        
    })
    
    
    ## update `meas_cond_user$df` with user inputs
    observeEvent(input$meas_cond_cell_edit, {
        
        info = input$meas_cond_cell_edit
        row_ind = info$row
        col_ind = info$col
        new_val = info$value
        
        # col 0 is the row-name column when rownames = TRUE
        if (col_ind > 0) {
            meas_cond_user$df[row_ind, col_ind] <-
                coerceValue(new_val, meas_cond_user$df[row_ind, col_ind])
        }
        
    })
    
    
    ## upon action button click `reset_tab`, reset values back to original
    observeEvent(input$reset_tab, {
        
        req(results_all_final())
        
        meas_cond_user$df <- default_meas_cond()
        
    })
    
    
    ## download table
    output$table2_exp <- downloadHandler(
        
        filename = paste("table_B_", format(Sys.time(), "%Y-%m-%d-%Hh%Mm%Ss"), ".csv", sep = ""),
        
        content = function(file){
            df <- meas_cond_user$df
            if (!is.null(df)) {
                write.table(df, file,
                            col.names = TRUE,
                            row.names = FALSE,
                            na = "",
                            sep = input$field_sep,
                            dec = input$dec_sep
                )
            }
        }
        
    ) 
    
    
    
    # Calculation of final concentrations of dissolved gases in the original samples
    
    ## upon action button click `calc_conc`:
    ## print out warning message before calculating final concentrations if temperature data is missing
    observeEvent(input$calc_conc, {
        
        if(!is.null(meas_cond_user$df) &&
           sum(is.na(meas_cond_user$df$temperature_C)) != 0){
            
            showNotification("Warning: Make sure to enter temperature values for all samples before proceeding to the calculations!",
                             type = "error",
                             duration = NULL)
            
        }
        
    })
    
    
    ## combine data of measured values from `results_all_final()` with user input data from `meas_cond_user$df`
    df_for_calc <- reactiveValues(df = NULL)
    
    observeEvent(input$calc_conc, {
        
        req(results_all_final(), meas_cond_user$df)
        
        df_for_calc$df <- isolate(results_all_final()) %>%
            select(c("group_name", "vial",
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
        
        req(df_for_calc$df)
        
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
            select(c("group_name", "vial", "sample_name",
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
            if (is.null(df_for_calc$df)){return ()}
            
            df <- final_conc_tab()
            if (!is.null(df)) {
                write.table(df, file,
                            col.names = TRUE,
                            row.names = FALSE,
                            na = "",
                            sep = input$field_sep,
                            dec = input$dec_sep
                )
            }
        }
    )
    
    
} # end server function


# Run the application 
shinyApp(ui = ui, server = server)
