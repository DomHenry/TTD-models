# TTD-models

The material provided here is the supplementary material for *Henry, D.A.W., Lee,. A.T.K. & Altwegg, R. (2019). Can time-to-detection models with fewer survey replicates provide a robust alternative to traditional site-occupancy models? Methods in Ecology and Evolution* 

The R code and JAGS model definitions are provided to reproduce the simulations and graphical outputs. We also define the models used to analyses the field data.

For the bird point count data used in the occupancy models, see the Dryad data repository: XXX

If you will be using this code, please cite this GitHub repository: XXX

Please cite this work as: XXX


## Run the simulation analysis

`sim_runs.R` is the primary file to recreate simulations  

`sim_ttd_data.R` is a function to generate time-to-detection data used in `sim_runs.R`.  


## Plot the simulation outputs 

`array_to_df.R` is a function to convert results of simulations (stored in arrays) to tidy data frames used for plotting.  

`sim_result_plots.R` is the primary file for generating Figs. 2, 3, S1, S2, S3 and S4 from the manuscript.  

                           
    
## JAGS models used in the simulations (single species)

`JAGS_sim_DND.R` is the model definition of the detection/non-detection occupancy model.

`JAGS_sim_TTD.R` is the model definition of the multi-survey time-to-detection occupancy model.  

`JAGS_sim_TTD_single.R` is the model definition of the single-survey time-to-detection occupancy model.  


## JAGS models used in analysis of bird field data (multi-species)

`JAGS_field_DND.R` is the model definition of the detection/non-detection occupancy model.  

`JAGS_field_TTD.R` is the model definition of the multi-survey time-to-detection occupancy model.  

`JAGS_field_TTD_single.R` is the model definition of the single-survey time-to-detection occupancy model.      

