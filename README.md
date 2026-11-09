# Picarro CRDS processing apps

This repository contains two Shiny apps for processing Picarro CRDS
measurements of headspace gas samples:

- `Single_file_app/` – processes raw, high-frequency `.dat` files containing
  one sample per file and performs automatic plateau detection.
- `Multi_sample_app/` – processes autosampler `*_Results.csv` files containing
  one sample per row, with plateau means and standard deviations already
  calculated.

Both apps apply the same dry correction, water-vapour correction, solubility
calculation and final dissolved-gas calculation. The main difference is the
input data and the starting point of the workflow.

## Requirements

- R (tested with R 4.1.2)
- The following R packages:

```r
install.packages(c("shiny", "DT", "ggplot2", "dplyr", "tidyr", "forcats"))
```

Uploads of large Picarro `.dat` files are supported; the upload limit is set to
500 MB in both apps.

## Running the apps

From the repository root:

```r
# Single-sample raw-data app
shiny::runApp("Single_file_app")

# Multi-sample autosampler app
shiny::runApp("Multi_sample_app")
```

Alternatively, open the respective `app.R` file in RStudio and click **Run App**.

## Single_file_app

### Purpose and input

This app processes raw Picarro CRDS `.dat` files in which each file contains a
time series for one sample. Multiple `.dat` files can be imported, and each file
is analysed individually before being added to a combined results table.

The sample ID used by the app is taken from the third dash-separated component
of the file name.

### Workflow

1. Import one or more `.dat` files.
2. Select the file to analyse.
3. Select the raw concentration columns:
   - 12C CO2 (dry)
   - 13C CO2 (wet)
   - 12C CH4 (dry), HP or HR
   - 13C CH4 (wet), HP or HR
4. Adjust the automatic plateau detection if necessary:
   - CO2 threshold concentration
   - Variation tolerance
   - Plateau time window slider
5. Review the time-series figures and the calculated plateau means and standard
   deviations.
6. Submit the current sample to the collected results table.
7. Enter the measurement conditions for each sample (temperature, liquid,
   headspace and injected volumes, salinity, pressure).
8. Calculate and export the final concentrations.

The figures show the measured 12C/13C CO2 and CH4 time series, highlight the
data points selected for analysis, and indicate the plateau window. Figures can
be exported as PDF, PNG or SVG.

## Multi_sample_app

### Purpose and input

This app processes the `*_Results.csv` files exported by the Picarro
autosampler software. Each row represents one sample or vial, and the plateau
mean and standard deviation for each measured species are already present in
the file.

Because plateau detection has already been performed by the autosampler
software, this step is skipped. The app starts from the reported means and
standard deviations and processes all sample rows at once.

### Workflow

1. Upload a `*_Results.csv` file.
2. Check the automatically detected CO2 columns:
   - 12C CO2 (dry) measured mean
   - 12C CO2 (wet) measured mean, used for the dry correction
   - 13C CO2 (wet) measured mean
3. Select the CH4 column set:
   - **Automatic** – HP below and HR above a user-defined threshold (default
     12 ppm), evaluated separately for every sample row.
   - **HP only** – use the HP columns for all samples.
   - **HR only** – use the HR columns for all samples.
4. Review the per-sample results table and exclude rows that should not be
   carried forward.
5. Continue to the editable measurement-conditions table. Temperature, liquid
   volume, headspace volume, injected volume, salinity and pressure are
   pre-filled from the sidebar default values and can be adjusted per sample.
6. Calculate and export the final concentrations.

An overview figure shows the calculated mean concentrations for every sample
row and can be exported as PDF, PNG or SVG.

## Shared calculations

### Dry correction

The 12C CO2 and 12C CH4 values are taken from dry-mole-fraction columns. The
13C values are measured wet and are converted to dry values using the matching
12C dry/wet ratio:

```text
13C_dry = 13C_wet × (12C_dry / 12C_wet)
```

In `Single_file_app`, this correction is applied to every measurement point. In
`Multi_sample_app`, it is applied to the reported means, and the SD of the
corrected 13C value is approximated by scaling the wet SD by the same ratio.

### Final concentrations

For each sample, the apps calculate:

- Final dissolved CO2 and CH4 concentrations in µmol/L and nmol/L, separately
  for 12C, 13C and total carbon.
- Isotope ratios, delta13C (‰) and the atom fraction of 13C (AF_13C).
- Technical standard deviations for the final concentrations.
- Water-vapour partial pressure and dry pressure from temperature, salinity
  and measured pressure.

The measurement conditions required for these calculations are:

- Sample temperature (°C)
- Liquid volume (mL)
- Headspace volume (mL)
- Injected gas volume (mL)
- Salinity (ppt)
- Pressure (atm)

In `Multi_sample_app`, all of these are set as default values in the sidebar
and applied to every sample. The values remain editable per sample in the
measurement-conditions table.

### Exports

Both apps provide CSV export of the results tables with selectable field and
decimal separators:

- Results table (including plateau statistics in the single-file app)
- Measurement-conditions table
- Final-concentrations table

## Notes

- In both apps, HP CH4 columns are recommended below approximately 12 ppm and
  HR columns above approximately 12 ppm.
- `Multi_sample_app` does not interpret the autosampler result flags or the
  "excluded by user" columns; these can be handled by excluding rows in the
  results table before continuing.
- The single-file app applies the dry correction point by point before
  averaging, while the multi-sample app applies it to the reported means
  because the raw time series are not present in the autosampler results file.
