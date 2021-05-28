# measles-regional-immunity
This repository contains the code and files to reproduce the analysis of our paper
> Alexis Robert, Adam J Kucharski, Helen C Johnson, Nick Bundle, 
> and Sebastian Funk, *The impact of local vaccine coverage and 
> recent incidence on measles transmission in France between 2009 
> and 2018*

This repository also contains scripts to generate a simulation study, aiming to assess the impact of aggregation on the Epidemic-Endemic framework. 

# Reproduce the analysis from the paper

The analysis presented in the paper uses case counts data reported to the European Center for Disease Prevention and Control, which is not publicly available. Therefore, in this repository, we first simulate the number of daily cases in each French department over a similar time span. Then, we fit the models developed using the Epidemic-Endemic framework with the same covariates as the models presented in the paper, and generate the figures. The data used to generate each covariate is located in the `Data` repository, and come from publicly available datasets. The details of each file and their sources is presented at the end of the README file. 

Two scripts are needed to generate all the figures from the analysis (Main and Supplement): 

*	`R/generate_analysis_paper.R`: Simulate all the data and fits needed to generate the figures in the main and the supplement:
    * Simulate the daily case counts.
    * Apply the report rate (i.e. keep 70% of the cases).
    *	Fit both models (i.e. where all regions can infect each other, or only neighbours).
    *	Simulate 1-year predictions using the model fits.
    *	Generate the sensitivity analyses (impact of the inferred vaccine coverage + impact of the proportion of direct transmission in the model fits).
    *	Compute the calibration scores 3,7,10, and 14 days ahead for both models, and the aggregated models.
    *	Save all the elements listed above.
*	`R/generate_all_figures_paper.R`: Use the outputs from the script “generate_analysis_paper.R” to generate all the figures from the Main Paper, and from the Supplement.

Given the various analyses needed to generate all the figures, the script `R/generate_analysis_paper.R` can take several hours to run. Therefore, all the outputs are saved in the list `Output/all_analysis.RDS`, which is imported at the beginning of `R/generate_all_figures_paper.R`.

Two other scripts are needed to run the models and the analysis:

*	`R/import_library.R`: list and import all the libraries required in the analysis
*	`R/import_covariates.R`: Import all functions, generate the time series of each covariate (coverage, population, surface), compute the serial interval, compute the distance matrices. The covariates and the distance matrices were computed using publicly available datasets, imported from the folder “Data”.

The functions required to run each script are contained in the following files:

*	`R/function_location_centroids.R`:  Defines the function `function_centroids`, which computes the weighted population centre of each French department using the R function zonal from the package raster, and the 1km<sup>2</sup> European Grid dataset.
*	`R/function_distance_population.R`: Defines the function `importation_distance` and `importation_pop_area`. The function `importation_distance` generates the distance between population centroids using the files computed in `function_centroids`, or computes which regions are neighbours (if the input parameter `neighbour` is TRUE).  The function `importation_pop_area` generates the number of inhabitant and the surface of each department using the Data file `Data/age_structure.csv`.
*	`R/function_coverage.R`: Contains the function `importation_coverage`, used to generate the three year average vaccine coverage in each department. The local vaccine uptake is imported from the Data file `Data/Vacc_coverage_departement.csv`. The missing values are then inferred using a Beta Mixed Model, and the function returns the time series of the local vaccine coverage.
*	`R/function_serial_interval.R`: Contains the function `conv` and `serial_interval`. `conv` is used to computed to convolution of a given function with itself. `serial_interval` is used to generate the serial interval of the disease, using mean and standard deviation as inputs, along with the proportion of transmission without missing generation. This function returns the distribution of the serial interval.
*	`R/function_generate_outbreak.R`: Contains the function `generate_sim` which generates one simulated outbreak using the Epidemic-Endemic framework, with the covariates used in the paper (coverage, population, recent incidence, surface, seasonality). Returns a list containing 8 elements: `n_cases` the number of daily cases per department, `incid` the level of recent incidence per department at each date, `r0_ar`, `r0_ne`, `r0_en` the value of the local predictor in each component, for each day, `n_cases_ar`, `n_cases_ne`, `n_cases_en` the average local number of cases stemming from each component at each date.
*	`R/function_generate_all_outbreaks.R`: Contains the function `generate_all_sim`, which generates a number of simulated outbreaks (defined in the input parameter `n_sim`). This function draws `n_sim` parameter sets, and uses `generate_sim` to simulate each outbreak. It returns a list containing 5 elements: `sim` a list containing the number of daily cases per department in each simulation, `r_ar`, `r_ne`, `r_end` a data frame describing the last value of each local predictor in each simulation, `params_sim` the parameter set used to generate each simulation, `prop_comp`: The proportion of cases stemming from each component for each simulation.
*	`R/function_report_rate.R`: Contains the function `report_rate`, used to remove a fraction of the generated cases, in order to mimic the partial detection of cases during an outbreak. Returns the number of daily cases reported per department.
*	`R/function_interpret_control_daily.R`: Contains a new version of the function `surveillance:::interpretControl`, re-defined to adapt the Epidemic-Endemic framework to the use of non-aggregated data. 
*	`R/function_exp_gravity.R`: Contains the function `W_exp_gravity_tot`, used to implement an exponential gravity model to construct the weight matrix connecting the regions included in the Epidemic-Endemic model. Follows the structure of the function `surveillance::W_powerlaw`.
*	`R/function_prepare_data_hhh4.R`: Contains the function `prep_data`, used to prepare the covariate before running the Epidemic-Endemic models. This function generates the category of incidence at each date and the transmission potential using the simulated case counts data, and ensures the columns and rows order are the same in each covariate. It returns a list containing all the covariates used in the Epidemic-Endemic models
*	`R/function_analysis_hhh4.R`: Contains the functions `function_hhh4_daily` and `function_hhh4_aggreg`, used to fit a daily or aggregated Epidemic-Endemic model to a set of simulated data. The functions use `prep_data` to prepare the covariates. The models are then run. These functions return a list containing `n_sim` elements, with each element corresponding to one model fit. 
*	`R/function_predict_1y.R`: Contains the functions `term_predict`, `generate_1y_outbreak`, `simulation_main_figures`, and `coordinates`, which are used to simulate future outbreaks after fitting the model, which will eventually be used to generate Figures 5 and 6. `term_predict` is used to integrate different conditions on the covariates (i.e. variations in vaccine coverage). `generate_1y_outbreak` is used to simulate one year of outbreak and returns the daily case counts per region. `simulation_mains_figures` uses `generate_1y_outbreak` to simulate the four simulation sets needed to generate Figures 5 and 6 (i.e. using the previous incidence, decreasing and increasing the vaccine coverage, and setting the level of previous incidence to the minimum). The simulations sets are run using the mean parameters and the covariance matrix from the model fits: the parameter `n_param` designates the number of parameter sets drawn, and `n_sim_param` describes the number of simulations per parameter set.  Finally, `coordinate` is used to simulate the impact of group importations on transmission. It returns the coordinates of importation, using the regions of the imported cases as inputs. The coordinates can then be used as an input to `figure_5_6` (parameter `import_loc`). 
*	`R/function_calibration.R`: Contains the functions `simulate_calib` and `calibration_model`: In the daily model, the calibration study relies on simulating the number of overall cases over the calibration period at each date. `simulate_calib` is used to simulate the number of daily cases per department over the calibration period, and `calibration_model` calls this function for each date of calibration, and compute the scores describing the match between the simulations and the data. `calibration_model` is also used to generate the calibration scores for aggregated model, in this case it uses the function `OneStepAhead` from the surveillance package.
*	`R/function_sensitivity.R`: Contains the functions `sens_SI`, `sens_vax`, and `sens_weekday`: Generate the model fits obtained under different conditions. In `sens_SI` the proportion of direct transmission in the composite serial interval varies (parameter `prop_gen1`); in `sens_vax` various values of the missing coverage data are generated before fitting the model; in `sens_weekday` a covariate is added to each component, quantifying the impact of weekends on the number of daily cases.
*	`R/function_generate_map.R`: Contains the function `generate_map`, used to import the shapefile of the metropolitan French department.

The objects generated by all these functions are used to generate the figures from the Main Paper and the Supplement. The functions used to generate the figures are contained in the following scripts:

*	`R/function_Main_figures.R`: Contains the functions `figure_1`, `figure_2`, figure_3`, figure_4`, and `figure_5_6`.
*	`R/function_Supplement_figures.R`: Contains the functions `figures_coverage`, `figures_sensitivity`, `figure_s9`, `figure_s11`, `figure_s12`, `figure_s15`, figure_s18`, and `figure_s19`

The analyses generated by the functions and scripts described above rely on various publicly available datasets, which are stored in the folder `Data`:

*	`Data/Vacc_coverage_departement.csv`: Vaccine coverage reported per department between 2004 and 2017 (excluding 2009). Aggregated from https://www.santepubliquefrance.fr/determinants-de-sante/vaccination/articles/donnees-departementales-2013-2017-de-couverture-vaccinale-rougeole-rubeole-oreillons-a-24-mois (accessed 28th May 2021), https://www.santepubliquefrance.fr/determinants-de-sante/vaccination/articles/donnees-departementales-2007-2012-de-couverture-vaccinale-rougeole-rubeole-oreillons-a-24-mois (accessed 28th May 2021), and https://www.santepubliquefrance.fr/determinants-de-sante/vaccination/articles/donnees-departementales-2013-2017-de-couverture-vaccinale-rougeole-rubeole-oreillons-a-24-mois (accessed 28th May 2021).
*	`Data/age_structure.csv`: Number of inhabitants (stratified by age group), per department between 2008 and 2019. Collected from: https://www.insee.fr/fr/statistiques/1893198, section « Estimation de population par département, sexe et grande classe d`âge » (accessed 28th May 2021).
*	`Data/location_centroids.RDS`: Generated using the function `function_centroids`` describes above.
*	`Data/nuts_to_dep.csv`: Table showing the correspondence between the NUTS ID and the name of each department. 
*	`Data/params_model_day_tot.RDS`: Mean estimates and covariance matrix of the parameters generated in the analysis using the French data (i.e. the parameter fits obtained in the analysis presented in the Main Paper). This object is used to generate the simulated case counts.

One file generated in the scripts is saved in the folder `Output`. 

*	`Output /all_analysis.RDS`: List generated using the script `R/generate_analysis_paper.R`, contains all objects needed to generate the figures from the Main paper and the Supplement. Since running the entirety of this script can take several hours, this object can be used to generate the figures directly.

# Analyse the impact of aggregation on the Epidemic-Endemic framework using a simulation study

The functions described in this folder were also used to generate the paper “Impact of aggregation on the Epidemic-Endemic framework: A simulation study”. The objective of this paper is to compare the fits obtained using a daily and an aggregated Epidemic-Endemic model, in order to assess the impact of aggregation on the parameter estimates, and the calibration. We generated 100 simulated case counts, fitted a daily and an aggregated model to each simulation, and compared the fits obtained with the input parameters. Most of the functions used in this project were described in the previous section. The scripts and functions developed specifically for this paper are:

*	`R/generate_simulations.R`: R script, can be run to generate and save a list containing 100 simulated case counts (saved as the object `Output/simulation_set.RDS`). 
*	`R/generate_analysis_simulations.R`: R script used to fit aggregated and daily Epidemic-Endemic models to the 100 simulations. This can take several hours, the models are grouped in lists and saved in the `Output` folder. 
*	`R/generate_calibration_scores.R`: R script used to compute the calibration scores of each model. Since we used a 10-day aggregation, we only compute the 10 days ahead calibration scores. The scores are computed over the last year of data (35 dates of calibration). This script can take several hours to run. The scores are therefore grouped in two lists, which are saved as RDS files in the `Output` folder.
*	`R/function_figures.R`: Contains the functions `plot_simulations`, `transp` `distance_params`, `in_CI`, `pit_param`, `plot_analysis`, `plot_prop_comp`, and `plot_calib_scores` used to generate the different figures from the paper. `plot_simulations` is used to generate figures 1 and 2; `plot_analysis` to generate figures 3, 4, 5, 6, and 8; `plot_prop_comp` to generate figure 7, and `plot_calib_scores` to generate figure 9.
*	`R/generate_plots_simulations.R`: R script used to generate all the figures included in the paper, using the simulations, fits, and calibrations computed during the previous scripts. 

Given some scripts take several hours to run, some of the files generated are saved in the `Output` folder:

*	`Output/simulation_set.RDS`: List containing 6 elements describing the simulated data: `sim` the list containing all the daily case counts; `r_ar`, `r_ne`, and `r_end` the last local value of each predictor for each simulation, `params_sim` a data frame containing the parameter set used to generate each simulation, and `prop_comp` the proportion of cases stemming from each component in each simulation.
*	`Output/models_aggregated.RDS`: List containing the 100 aggregated models fitted to the simulations.
*	`Output/models_daily1.RDS`, `models_daily2.RDS`, `models_daily3.RDS`, `models_daily4.RDS`: Lists containing the 100 daily models fitted to the simulations. The 100 models were split in four lists, otherwise the object was too heavy for Github commits.
*	`Output/scores_aggreg`: List containing 6 data frames describing the calibration scores of the aggregated models: `bias` the value of bias for each calibration point, `sharpness` the value of sharpness for each calibration point, `rps` the ranked probability score for each calibration point,  `log` the log score for each calibration point, `px` the probability that the number of cases in the simulation is equal to the data, `pxm1` the probability that the number of cases in the simulation is equal to the data - 1.
*	`Output/scores_day`: List containing 7 data frames describing the calibration scores of the aggregated models: `bias` the value of bias for each calibration point, `sharpness` the value of sharpness for each calibration point, `rps` the ranked probability score for each calibration point,  `log` the log score for each calibration point, `px` the probability that the number of cases in the simulation is equal to the data, `pxm1` the probability that the number of cases in the simulation is equal to the data - 1, `tot` the median, 95% and 50% credible intervals for the number of cases generated in the calibration.
