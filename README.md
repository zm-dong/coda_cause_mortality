$Data:$

GBRTENW_d_interm_idr_agg.csv: mortality by case from HMD for England and Wales from 2001 to 2016, intermediate classifications of cause.
mltper_5x1_UK.txt: male life tables for UK (England and Wales). 5 year age bands, annual data
fltper_5x1_UK.txt: females life tables for UK (England and Wales). 5 year age bands, annual data

$Scripts:$

01F_CoDA_Mortality_UK_IntermAgg_LC_clean.R: this file is used to produce results using the CLR, ILR, and alpha transformations on UK and Wales mortality by cause data for females
02F_CoDA_LCplots_clean.R: run this code after 01F to produce plots for females

01M_CoDA_Mortality_UK_IntermAgg_LC_clean.R: this file is used to produce results using the CLR, ILR, and alpha transformations on UK and Wales mortality by cause data for males
02M_CoDA_LCplots_clean.R: run this code after 01F to produce plots for males

$Functions:$

alfainv_mod.R: modified alfainv function for use with this data set which has a number of border cases
save_function.R: funciton used to save charts and plots produced using scripts

$Misc:$

session_info.rtf: information about the original set up used to run the code and produce results
