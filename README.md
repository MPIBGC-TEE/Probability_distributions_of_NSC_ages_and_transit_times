
#title: Probability distributions of non-structural carbon ages and transit times provide insights in carbon allocation dynamics of mature trees

author: "David_Herrera"
date: "1/5/2020"



In this repository you will find the code used to compute the results published in the manuscript "Probability distributions of non-structural carbon ages and transit times provide insights in carbon allocation dynamics of mature trees"


## Packages to install

For our calculations we used open packages in R and Python. For running the scripts included in this repository you should install the following packages. 

### R Packages 

* SoilR

### Python packages 

* CompartmentalSystems. You can find instructions about how to install this package in https://github.com/MPIBGC-TEE/CompartmentalSystems/blob/master/README.md

## Files description 

In the repository there are two main folders: 1) Code, it contains all the files and scripts with the calculations for obtaining the same results we published in the manuscript; 2) Documents_manuscript, it contains all the documents related with the manuscript itself, including all the versions of the manuscript generated during the review process, and also the reviewer comments were included. 

### Code folder

There are two files and three folders. The two files have all the general R functions necessary to run the scripts in each of the species folders. So, these scripts are always sourced from each species script. The script "script_function_senes_uncert_ACGCA.R" contains the functions for running the sensitivity analysis, and the script "script_functions_model_ACGCA.R" contains the functions for running the ACGCA model and calculating the model parameters. 

#### *A. rubrum* and *P. taeda* folders 

Each folder contains the scripts for running the computations for each individual species, the figures presented in the manuscript and the information generated in each simulation. Some of the simulations take hours for running, therefore we load these data in the scripts for the subsequent calculations. But simulations can be run again if desired. The main scripts in each species folder are: 

* (species_name)_ages_and_transit_time.ipynb: This is the Jupyter Notebook where the age and transit time calculations were done. Also the NSC carbohydrates consumption for each year of the simulation is calculated in this document.

* (species_name)_parameters_sensitivty_and_uncertainty_analysis.Rmd: In this R markdown you can find the code for running the sensitivity and uncertainty analyses related to the species name.

Please make sure to adjust the working directory in each script, for loading the data and sourcing the functions. 

#### *P. helepensis* folder 

Here you will find four scripts:

* Functions_Phalepensis.R: Here some specific functions for running the model for _P. halepensis_ are written 

* Klain_2015_annual_parameter_estimation.R: Here are the calculations for estimating the annual transfer coefficients (parameters) for the _P. halepens_ model based on the values reported by Klein and Hoch (2015) paper. 

* Phalepensis_parameters_sensitivity_and_uncertainty_analysis.R: Here the sensitivity and uncertainty analyses for the _P. halepensis_ model are written. 

* P.halepenesis_age_and_transit_times.ipynb: This is the Jupyter Notebook where the age and transit time calculations were computed. 

The .Rda files contain the results from the different simulations coded in the R scripts. These results are all loaded in the scripts to avoid the long computation time that some simulations take. 
