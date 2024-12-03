# Multiple Imputation in Cure Models

This repository aims to make the results from "A multiple imputation approach to distinguish curative from life-prolonging effects in the presence of missing covariates" reproducible. 
There are two directories, with R_Simulations showing the code for the simulation study (exactly reproducible) and R_Analysis showing the code used to obtain results for the BO06 clinical trial data.

The two subdirectories are described in more detail below.


# R_Simulations
In this directory, the code used to perform the simulation studies can be found. There are the following files:
- "imputation_cure_similation" is the main code reproducing the algorithm
- "imputation_cure_similation_4" is the adjusted algorithm code for Scenario D
- "functions" contains all the necessary functions
- "functions_4" contains the customized functions for Scenario D
- "parameters" is the file from which to select the parameters set for each scenario.


# R_Analysis
In this directory is contained the code which was used to perform analyses on the BO06 clinical trial data. 
The results are not exactly reproducible as the data is not publicly available. 
