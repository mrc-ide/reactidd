# reactidd

This repository supports epidemiological and disease-dynamic analyses of data from the REal Time Assessment of Community Transmission (REACT) study. It includes both code and data. 

REACT is a program of different studies with multiple rounds, organised into two groups: REACT-1 is a community sample of swab-positivity in which participants are asked to swab themselves, arange for a courier to pickup that swab and then fill out a questionnaire. The swabs are then tested in a lab using PCR. REACT-2 is a community sample of antibody-positivity in which participants are asked to tst themselves using  lateral flow test (LFT) and then report that result at the same time they fill out the questionnaire. Data from both are included here. 

The repository is structured as an R package and is most easily installed using the `install_github` in the package `dev_tools`. The package requires an R enrironment that can build from source. Also, some of the functions rely on the package `rstan` which in trun needs the `stan` library to be installed on your system. However the data are also directly available from the `inst\extdata`as `csv` files.

## Temporal data

The main temporal data for REACT are 

## Spatial data

## Risk factor data
