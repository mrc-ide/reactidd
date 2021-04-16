# reactidd

This repository supports epidemiological and disease-dynamic analyses of data from the REal Time Assessment of Community Transmission (REACT) study. It includes both code and data. 

REACT is a program of different studies with multiple rounds, organised into two groups: REACT-1 is a community sample of swab-positivity in which participants are asked to swab themselves, arange for a courier to pickup that swab and then fill out a questionnaire. The swabs are then tested in a lab using PCR. REACT-2 is a community sample of antibody-positivity in which participants are asked to tst themselves using  lateral flow test (LFT) and then report that result at the same time they fill out the questionnaire. Data from both are included here. 

The repository is structured as an R package and is most easily installed using the `install_github` in the package `dev_tools`. The package requires an R enrironment that can build from source. Also, some of the functions rely on the package `rstan` which in trun needs the `stan` library to be installed on your system. However the data are also directly available from the `inst\extdata`as `csv` files.

## Temporal data

The main temporal data for REACT are `inst/extdata`. The file `positive.csv` contains the number of positive swabs collected by day and by region for all currently reported rounds of REACT-1. Similarly, the file `total.csv` contains the total number of swabs collected by day and by region. 

The vignette `TemporalAnalysisREACT.rmd` demonstrates how the REACT data can be loaded, exponential models fit/plotted, estimates of growth rate/R calculated, and p-spline models fit/plotted. The vignette `TemporalAnalysisPHE.rmd` demonstrates how similar analyses can be performed on publically avaialble PHE case data.

# Spatio-temporal geostatistical model

The spatial analyses is contained in the directory `inst/spatial`. Each of the R scripts there can be run in sequence to regenerate spatial output for rounds 1 to 4, using the geospatial modelling framework.

This code was tested with R version 3.6.3 on Ubuntu 18.04 LTS.

To run it download the content of the folder `inst/spatial` and make it your working directory, then run these scripts:

- `01_data-cleaning.R`: cleans and processes the input data in `data/original` and gets it ready for modelling. The output of this script goes into the `data/processed` folder.

- `02_model-fitting.R`: uses the processed data to fit the spatio-temporal model and save it in `output/models/`. 

- `03_predictions.R`: uses the fitted spatio-temporal model to generate predictions and and summarise
them for mapping. All outputs will be saved in `output/predictions/`.

- `04_mappings.R`: create the maps reported in the paper, these will be saved as pdfs in the `figs` folder.

- `05_tables.R`: this creates a summary table with the estimated model parameters.

The `R` folder contains the `functions.R` file that has a set of costum functions needed to run the scripts above.
