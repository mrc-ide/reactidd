# reactidd

This repository supports epidemiological and disease-dynamic analyses of data from the REal Time Assessment of Community Transmission (REACT) study. It includes both code and data. 

REACT is a program of different studies with multiple rounds, organised into two groups: REACT-1 is a community sample of swab-positivity in which participants are asked to swab themselves, arange for a courier to pickup that swab and then fill out a questionnaire. The swabs are then tested in a lab using PCR. REACT-2 is a community sample of antibody-positivity in which participants are asked to tst themselves using  lateral flow test (LFT) and then report that result at the same time they fill out the questionnaire. Data from both are included here. 

The repository is structured as an R package and is most easily installed using the `install_github` in the package `dev_tools`. The package requires an R enrironment that can build from source. Also, some of the functions rely on the package `rstan` which in trun needs the `stan` library to be installed on your system. However the data are also directly available from the `inst\extdata`as `csv` files.

## Temporal analyses

The main temporal data for REACT are `inst/extdata`. The file `positive.csv` contains the number of positive swabs collected by day and by region for all currently reported rounds of REACT-1. Similarly, the file `total.csv` contains the total number of swabs collected by day and by region. 

The vignette `TemporalAnalysisREACT.rmd` demonstrates how the REACT data can be loaded, exponential models fit/plotted, estimates of growth rate/R calculated, and p-spline models fit/plotted. The vignette `TemporalAnalysisPHE.rmd` demonstrates how similar analyses can be performed on publically avaialble PHE case data.

## Notes for individual publications

### [Resurgence of SARS-CoV-2: Detection by Community Viral Surveillance](http://dx.doi.org/10.1126/science.abf0874)

The two vignettes `TemporalAnalysisREACT.rmd` and `TemporalAnalysisPHE.rmd` demonstrate the temporal analysis used in this publication

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

### [REACT-1 round 13 final report: exponential growth, high prevalence of SARS-CoV-2 and vaccine effectiveness associated with Delta variant in England during May to July 2021](https://doi.org/10.1101/2021.09.02.21262979)

The vignette `TemporalAnalysisREACT_rounds12and13.rmd` demonstrates the temporal analysis used in this publication

### [Characterising the persistence of RT-PCR positivity and incidence in a community survey of SARS-CoV-2](http://dx.doi.org/10.12688/wellcomeopenres.17723.1)

The vignette `EstimatingDurationOfSwabPositivity.rmd` in the vignette subfolder `PCR_Positivity_Paper` contains the code used in the analysis for this preprint. The analysis perfomed on the data for repeat tests is demonstrated on simulated data as the indiviudal level data could not be shared due to ethical/security considerations. Also in the vignette subfolder named `PCR_Positivity_Paper` is the extended data to support the submission of the paper to Wellcome Open. The files included are 1) `COG_UK authorship.xlsx`, which contains the author deatils for the Covid-19 Genomics UK (COG-UK) consortium 2) `SupplementaryFigure1.pdf` which contains supplementary figure 1 3) `SupplementaryTable1.xlsx` which contains supplementary table 1 and 4) `Extended data descriptions.docx` which contains the legends for each supplementary materials.

### [Appropriately smoothing prevalence data to inform estimates of growth rate and reproduction number (preprint)](https://doi.org/10.1101/2022.02.04.22270426)

The vignettes `REACT_rounds1-7_analysis.rmd` and `PHE_rounds1-7_analysis.rmd` in the vigentte subfolder `TemporalMethodsPaper` demonstrate the temporal analysis used in this preprint.


### [The new normal? Dynamics and scale of the SARS-CoV-2 variant Omicron epidemic in England (preprint)](https://www.medrxiv.org/content/10.1101/2022.03.29.22273042v1)

The vignette `REACT_rounds14-18_omicron_analysis.rmd` in the vigentte subfolder `TemporalOmicronPaper` demonstrate the temporal analysis used in this preprint. The code allows the analysis of two competing variants when overall prevalence is known and the daily proportion of both competing variants is known.
