$Data:$

GBRTENW_d_interm_idr_agg.csv: death counts by cause from HMD for England and Wales from 2001 to 2016, intermediate classifications of cause.
mltper_5x1_UK.txt: male life tables for UK (England and Wales). 5 year age bands, annual data
fltper_5x1_UK.txt: females life tables for UK (England and Wales). 5 year age bands, annual data

USA_d_long_idr_agg.csv: death counts by cause from HMD for US from 1979 to 2021, grouped causes disaggregated for cardiovascular related deaths.
mltper_5x1_USA.txt: male life tables for US. 5 year age bands, annual data
fltper_5x1_USA.txt: female life tables for US. 5 year age bands, annual data
USA_population_5.txt: US population data

$Scripts:$

01F_CoDA_Mortality_UK_IntermAgg_LC_clean.R: this file is used to produce results using the CLR, ILR, and alpha transformations on UK and Wales mortality by cause data for females
02F_CoDA_LCplots_clean.R: run this code after 01F to produce plots for females

01M_CoDA_Mortality_UK_IntermAgg_LC_clean.R: this file is used to produce results using the CLR, ILR, and alpha transformations on UK and Wales mortality by cause data for males
02M_CoDA_LCplots_clean.R: run this code after 01F to produce plots for males

CoDA_Mortality_US_LongAgg_LC_v2_male.R: this file is used to produce results for US data, males
CoDA_Mortality_US_LongAgg_LC_v2_female.R: this file is used to produce results for US data, females
CoDA_MLR_Mortality_US_v1_male.R: this file is used to produce comparative results using the MLR approach

$Functions:$

alfainv_mod.R: modified alfainv function for use with this data set which has a number of border cases
save_function.R: funciton used to save charts and plots produced using scripts

Interval_Forecasts_alpha_v0.1.R: function to produce interval forecasts using the alpha-transformation
Interval_Forecasts_CLR_v0.1.R: function to produce interval forecasts using CLR
Interval_Forecasts_ILR_v0.1.R: function to produce interval forecasts using ILR

$Misc:$

session_info.rtf: information about the original set up used to run the code and produce results
