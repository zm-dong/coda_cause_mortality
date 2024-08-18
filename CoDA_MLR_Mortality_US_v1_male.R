# This analysis looks at the US mortality by cause data
# Long list of causes (aggregated)
# We apply the alpha-transformation to this data set, where alpha = 0, we get the isometric log ratios
rm(list=ls())

#load packages
setwd("~/Documents/PhD/03_CODA/Model_US")
source("function.r")
source("save_function.r")
source("generate_mortality_tensor.r")

library(demography)
library(stats)
library(strucchange)
library(tables)
library(calibrate)
library(compositions)
library(abind)
library(MLmetrics)
library(tseries)
library(vars)
library(rTensor)
library(r.jive)
library(FactoMineR)
library(Metrics)
library(rainbow)
library(shape)
library(ggplot2)
library(reshape2)
library(viridis)
library(patchwork)
library(RColorBrewer)
library(rcartocolor)
library(paletteer)

library(nnet)
library(tidyverse)

library(reshape2)
library(data.table)
library(rTensor)

options(scipen=999)

# ================================================================
## read original data: deaths by cause
# ================================================================
dati0 <- read.csv("USA_d_long_idr_agg.csv", header=TRUE)

## take out total all causes (cause==0)
dati_step0 <- subset(dati0, cause!=0)
dati_step1 <- dati_step0[,-1]
dati_step2 <- dati_step1[,-3]
dati_step3 <- dati_step2[,-3]
dati_step4 <- dati_step3[,-4]
dati_step5 <- dati_step4[,-14]
dati_step6 <- dati_step5[,-13]
dati <- dati_step6

## domains and dimensions
ages <- c(0,seq(30, 90, 10), 100)

# ================================================================
# for use with the old version / England and Wales data
#aggregate causes: sum everything that is not within cause 48 to 57, keep 48 to 57 disaggregated
#aggregate years: sum everything below age 25, and then have 10 year age bands from 25 to 95, then keep age 95+
# ================================================================

## take the mean of each interval,
## and 95 for the last open interval
m <- length(ages)
years <- unique(dati$year)
n <- length(years)
CoD <- unique(dati$cause)
k <- length(CoD)

## names of the CoD (sorted as CoD)
# Refer to data formats file from cause of death database for ICD10 mappings
## data separate by sex, in arrays
#males
datiM0 <- as.matrix(subset(dati, sex=="1")[,1:m+3])
#females
datiF0 <- as.matrix(subset(dati, sex=="2")[,1:m+3])
#total
datiT0 <- as.matrix(subset(dati, sex=="3")[,1:m+3])

# ================================================================
# required for MLR: reshaping data
# ================================================================
usa_population_0 <-(read.table("USA_population_5.txt",skip=1,header=TRUE,stringsAsFactors=F))

usa_population_m <- as.matrix(subset(usa_population_0)[,c(1:2, 4)])
usa_population_f <- as.matrix(subset(usa_population_0)[,c(1:2, 3)])

usa_pop_m <- matrix(as.numeric(usa_population_m[,"Male"]),length(unique(usa_population_m[,"Age"])),length(unique(usa_population_m[,"Year"])))
colnames(usa_pop_m) <- unique(usa_population_m[,"Year"])
rownames(usa_pop_m) <- unique(usa_population_m[,"Age"])

# Take a subset of the years based on what is available by cause
usa_pop_m1 <- usa_pop_m[,colnames(usa_pop_m) %in% years]
rownames(usa_pop_m1) <- unique(usa_population_m[,"Age"])

age.temp<- unique(usa_population_m[,"Age"])

usa_pop_m2 <- rbind(usa_pop_m1[1,]+usa_pop_m1[2,]+usa_pop_m1[3,]+usa_pop_m1[4,]+usa_pop_m1[5,]+usa_pop_m1[6,]+usa_pop_m1[7,],
                   usa_pop_m1[8,]+usa_pop_m1[9,],
                   usa_pop_m1[10,]+usa_pop_m1[11,],
                   usa_pop_m1[12,]+usa_pop_m1[13,],
                   usa_pop_m1[14,]+usa_pop_m1[15,],
                   usa_pop_m1[16,]+usa_pop_m1[17,],
                   usa_pop_m1[18,]+usa_pop_m1[19,],
                   usa_pop_m1[20,]+usa_pop_m1[21,],
                   usa_pop_m1[22,]+usa_pop_m1[23,]+usa_pop_m1[24,])

rownames(usa_pop_m2) <- ages

usa_pop_m3 <- reshape2::melt(data = usa_pop_m2, id.vars = c(ages), measure.vars  = c(years), value.name    = "population")
colnames(usa_pop_m3) <- c("ages", "years", "popm")

datiM_0 <- as.matrix(subset(dati, sex=="1"))
datiM_1 <- datiM_0[,c(1, 3:12)]
colnames(datiM_1) <- c("years", "cause", ages)

datiM_2 <- as.data.table(datiM_1)
datiM_3 <- reshape2::melt(data = datiM_2, id = 1:2, measure = 3:11, variable.name = "ages")

datiM_c1 <- subset(datiM_3, cause == 1, select = c(years, cause, ages, value))
datiM_c102 <- subset(datiM_3, cause == 102, select = c(years, cause, ages, value))
datiM_c103 <- subset(datiM_3, cause == 103, select = c(years, cause, ages, value))
datiM_c104 <- subset(datiM_3, cause == 104, select = c(years, cause, ages, value))
datiM_c105 <- subset(datiM_3, cause == 105, select = c(years, cause, ages, value))
datiM_c106 <- subset(datiM_3, cause == 106, select = c(years, cause, ages, value))
datiM_c107 <- subset(datiM_3, cause == 107, select = c(years, cause, ages, value))
datiM_c108 <- subset(datiM_3, cause == 108, select = c(years, cause, ages, value))
datiM_c109 <- subset(datiM_3, cause == 109, select = c(years, cause, ages, value))
datiM_c110 <- subset(datiM_3, cause == 110, select = c(years, cause, ages, value))
datiM_c111 <- subset(datiM_3, cause == 111, select = c(years, cause, ages, value))
datiM_c112 <- subset(datiM_3, cause == 112, select = c(years, cause, ages, value))

usa_pop_m4 <- transform(usa_pop_m3, ages = as.factor(ages))

deaths_1 <- left_join(usa_pop_m4, datiM_c1, by = c("ages", "years"))
colnames(deaths_1) <- c("ages", "years", "popm", "cause", "d_oth")

deaths_2 <- left_join(deaths_1, datiM_c102, by = c("ages", "years"))
deaths_2 <- deaths_2[,c(1:3, 5, 7)]
colnames(deaths_2) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute")

deaths_3 <- left_join(deaths_2, datiM_c103, by = c("ages", "years"))
deaths_3 <- deaths_3[,c(1:5,7)]
colnames(deaths_3) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic")

deaths_4 <- left_join(deaths_3, datiM_c104, by = c("ages", "years"))
deaths_4 <- deaths_4[,c(1:6,8)]
colnames(deaths_4) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension")

deaths_5 <- left_join(deaths_4, datiM_c105, by = c("ages", "years"))
deaths_5 <- deaths_5[,c(1:7,9)]
colnames(deaths_5) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension", "d_hypertensive_heart")

deaths_6 <- left_join(deaths_5, datiM_c106, by = c("ages", "years"))
deaths_6 <- deaths_6[,c(1:8,10)]
colnames(deaths_6) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension", "d_hypertensive_heart", 
                        "d_hypertensive_renal")

deaths_7 <- left_join(deaths_6, datiM_c107, by = c("ages", "years"))
deaths_7 <- deaths_7[,c(1:9,11)]
colnames(deaths_7) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension", "d_hypertensive_heart", 
                        "d_hypertensive_renal", "d_hypertensive_hr")

deaths_8 <- left_join(deaths_7, datiM_c108, by = c("ages", "years"))
deaths_8 <- deaths_8[,c(1:10,12)]
colnames(deaths_8) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension", "d_hypertensive_heart", 
                        "d_hypertensive_renal", "d_hypertensive_hr", "d_myoc_inf")

deaths_9 <- left_join(deaths_8, datiM_c109, by = c("ages", "years"))
deaths_9 <- deaths_9[,c(1:11,13)]
colnames(deaths_9) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension", "d_hypertensive_heart", 
                        "d_hypertensive_renal", "d_hypertensive_hr", "d_myoc_inf", "d_IHD_acute")

deaths_10 <- left_join(deaths_9, datiM_c110, by = c("ages", "years"))
deaths_10 <- deaths_10[,c(1:12,14)]
colnames(deaths_10) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension", "d_hypertensive_heart", 
                        "d_hypertensive_renal", "d_hypertensive_hr", "d_myoc_inf", "d_IHD_acute", "d_IHD_chronic")


deaths_11 <- left_join(deaths_10, datiM_c111, by = c("ages", "years"))
deaths_11 <- deaths_11[,c(1:13,15)]
colnames(deaths_11) <- c("ages", "years", "popm", "d_oth", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension", "d_hypertensive_heart", 
                         "d_hypertensive_renal", "d_hypertensive_hr", "d_myoc_inf", "d_IHD_acute", "d_IHD_chronic", "d_pulmonary")

deaths_12 <- left_join(deaths_11, datiM_c112, by = c("ages", "years"))
deaths_12 <- deaths_12[,c(1:14,16)]
colnames(deaths_12) <- c("ages", "years", "popm", "d_oth_all", "d_rheumatic_acute", "d_rheumatic_chronic", "d_hypertension", "d_hypertensive_heart"
                         , "d_hypertensive_renal", "d_hypertensive_hr", "d_myoc_inf", "d_IHD_acute", "d_IHD_chronic", "d_pulmonary"
                         , "d_oth_cardio")


deaths_12$all <- rowSums(deaths_12[4:15])
deaths_12$surv <- deaths_12$popm - deaths_12$all

deaths_m <- as.matrix(deaths_12[, c(4:17)])

# fit mlr
years <- (1:43) %>%
  scale(center = FALSE, scale = FALSE)
ages <- c(ages)
dat <- data.frame(years = rep(years, each = length(unique(deaths_12$ages))), agegroups = rep(ages, length(unique(deaths_12$years))))
head(dat)

make_offset <- matrix(deaths_12$popm, nrow = nrow(deaths_12), ncol = ncol(deaths_12)-3)

# simplest
fit_mlr1s <- multinom(deaths_m[1:387, ] ~ years + agegroups + offset(make_offset[1:387, ]), data = dat[1:387, ], maxit = 1000)
# single
fit_mlr1 <- multinom(deaths_m[1:387, ] ~ years + agegroups + years:agegroups + offset(make_offset[1:387, ]), data = dat[1:387, ], maxit = 1000)
# quadratic
fit_mlr2 <- multinom(deaths_m[1:387, ] ~ years + agegroups + I(years):agegroups + I(years^2) + I(years^2):agegroups + offset(make_offset[1:387, ]), data = dat[1:387, ], maxit = 1000)
# cubic
fit_mlr3 <- multinom(deaths_m[1:387, ] ~ years + agegroups + I(years^2) + I(years^2):agegroups + I(years):agegroups + I(years^3) + I(years^3):agegroups + offset(make_offset[1:387, ]), data = dat[1:387, ], maxit = 5000)

summary(fit_mlr1s)
summary(fit_mlr1)
summary(fit_mlr2)
summary(fit_mlr3)

coefs1s <- coef(fit_mlr1s)
coefs1 <- coef(fit_mlr1)
coefs2 <- coef(fit_mlr2)
coefs3 <- coef(fit_mlr3)


## --------------------
## --------------------
## fit
## --------------------
## --------------------
fit_dat <- data.frame(years = rep(c(1:43), each = length(ages)), agegroups = rep(ages, 43))
X_fit1s <- model.matrix(~ years + agegroups, data = fit_dat)
X_fit1 <- model.matrix(~ years + agegroups + years:agegroups, data = fit_dat)
X_fit2 <- model.matrix(~ years + agegroups + I(years):agegroups + I(years^2) + I(years^2):agegroups, data = fit_dat)
X_fit3 <- model.matrix(~ years + agegroups + I(years^2) + I(years^2):agegroups + I(years):agegroups + I(years^3) + I(years^3):agegroups, data = fit_dat)

## --------------------
## simplest
## --------------------
fit_probs1s <- tcrossprod(X_fit1s, coefs1s) %>%
  exp()

fit_probs1s <- fit_probs1s / (1 + rowSums(fit_probs1s))
  
  #cbind(1 / (1 + rowSums(fit_probs1s)), fit_probs1s / (1 + rowSums(fit_probs1s)))
#colnames(fit_probs1s)[1] <- "Survival"

fit_probs1s <- cbind(fit_dat, fit_probs1s[, 1:12])
colnames(fit_probs1s)[1] <- "Year"

df_wide1s <- as.data.table(fit_probs1s) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(fit_probs1s))[3:14])

fit1s <- df_wide1s %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

## --------------------
## single
## --------------------
fit_probs1 <- tcrossprod(X_fit1, coefs1) %>%
  exp()
fit_probs1 <- fit_probs1 / (1 + rowSums(fit_probs1))
fit_probs1 <- cbind(fit_dat, fit_probs1[, 1:12])
colnames(fit_probs1)[1] <- "Year"

df_wide1 <- as.data.table(fit_probs1) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(fit_probs1))[3:14])

fit1 <- df_wide1 %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

## --------------------
## quadratic
## --------------------
fit_probs2 <- tcrossprod(X_fit2, coefs2) %>%
  exp()
fit_probs2 <- fit_probs2 / (1 + rowSums(fit_probs2))
fit_probs2 <- cbind(fit_dat, fit_probs2[, 1:12])
colnames(fit_probs2)[1] <- "Year"

df_wide2 <- as.data.table(fit_probs2) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(fit_probs2))[3:14])

fit2 <- df_wide2 %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

#fit2 <- fit2[, colnames(usam[, c(2:22)])]
##### difference for quadratic
#fnorm(generate_mortality_tensor(as.data.table(fit2)) - generate_mortality_tensor(as.data.table(actuals)))

## --------------------
## cubic
## --------------------
fit_probs3 <- tcrossprod(X_fit3, coefs3) %>%
  exp()
fit_probs3 <- fit_probs3 / (1 + rowSums(fit_probs3))
fit_probs3 <- cbind(fit_dat, fit_probs3[, 1:12])
colnames(fit_probs3)[1] <- "Year"

df_wide3 <- as.data.table(fit_probs3) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(fit_probs3))[3:14])

fit3 <- df_wide3 %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

#fit3 <- fit3[, colnames(usam[, c(2:22)])]


## --------------------
## --------------------
## predict
## --------------------
## --------------------
pre_dat <- data.frame(years = rep(c(44:48), each = length(unique(ages))), agegroups = rep(ages, 5))
X_pred1s <- model.matrix(~ years + agegroups, data = pre_dat)
X_pred1 <- model.matrix(~ years + agegroups + years:agegroups, data = pre_dat)
X_pred2 <- model.matrix(~ years + agegroups + I(years):agegroups + I(years^2) + I(years^2):agegroups, data = pre_dat)
X_pred3 <- model.matrix(~ years + agegroups + I(years^2) + I(years^2):agegroups + I(years):agegroups + I(years^3) + I(years^3):agegroups, data = pre_dat)

## --------------------
## simplest
## --------------------
pred_probs1s <- tcrossprod(X_pred1s, coefs1s) %>%
  exp()
pred_probs1s <- pred_probs1s / (1 + rowSums(pred_probs1s))
  #cbind(1 / (1 + rowSums(pred_probs1s)), pred_probs1s / (1 + rowSums(pred_probs1s)))
#colnames(pred_probs1s)[1] <- "Survival"
pred_probs1s <- cbind(pre_dat, pred_probs1s[, 1:12])
colnames(pred_probs1s)[1] <- "Year"

df_wide1ps <- as.data.table(pred_probs1s) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(pred_probs1s))[3:14])

pred1s <- df_wide1ps %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

#pred1s <- pred1s[, colnames(usam[, c(2:22)])]
##### difference
#fnorm(generate_mortality_tensor(pred1s) - generate_mortality_tensor(usam[(53 * 6 + 1):(58 * 6), ]))

## --------------------
## single
## --------------------
pred_probs1 <- tcrossprod(X_pred1, coefs1) %>%
  exp()
pred_probs1 <- pred_probs1 / (1 + rowSums(pred_probs1))
pred_probs1 <- cbind(pre_dat, pred_probs1[, 1:12])
colnames(pred_probs1)[1] <- "Year"

df_wide1p <- as.data.table(pred_probs1) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(pred_probs1))[3:14])

pred1 <- df_wide1p %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

#pred1 <- pred1[, colnames(usam[, c(2:22)])]
##### difference
#fnorm(generate_mortality_tensor(pred1) - generate_mortality_tensor(usam[(53 * 6 + 1):(58 * 6), ]))

## --------------------
## quadratic
## --------------------
pred_probs2 <- tcrossprod(X_pred2, coefs2) %>%
  exp()
pred_probs2 <- pred_probs2 / (1 + rowSums(pred_probs2))
pred_probs2 <- cbind(pre_dat, pred_probs2[, 1:12])
colnames(pred_probs2)[1] <- "Year"

df_wide2p <- as.data.table(pred_probs2) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(pred_probs2))[3:14])

pred2 <- df_wide2p %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

#pred2 <- pred2[, colnames(usam[, c(2:22)])]
##### difference
#fnorm(generate_mortality_tensor(pred2) - generate_mortality_tensor(usam[(53 * 6 + 1):(58 * 6), ]))

## --------------------
## cubic
## --------------------
pred_probs3 <- tcrossprod(X_pred3, coefs3) %>%
  exp()
pred_probs3 <- pred_probs3 / (1 + rowSums(pred_probs3))
pred_probs3 <- cbind(pre_dat, pred_probs3[, 1:12])
colnames(pred_probs3)[1] <- "Year"

df_wide3p <- as.data.table(pred_probs3) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(pred_probs3))[3:14])

pred3 <- df_wide3p %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

#pred3 <- df_wide3p %>%
#  gather("variable", "value", -Year) %>%
#  mutate(
#    ages = sub(".*_", "", variable)
#    , variable = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable)
#    , variable = sub("_[[:digit:]][[:digit:]]", "", variable)
#    , variable = sub("_[[:digit:]]", "", variable)
#  ) %>%
#  spread(variable, value)

#pred3 <- pred3[, colnames(usam[, c(2:22)])]
##### difference
#fnorm(generate_mortality_tensor(pred3) - generate_mortality_tensor(deaths_12[(43 * 12 + 1):(48 * 12), ]))

#causes <- colnames(pred3[,3:14])

#actual death rates by cause
deaths_nozero <- deaths_12
deaths_nozero[deaths_nozero == 0] <- 0.01

actual_mort <- cbind(Year = deaths_12[, 2], agegroups = deaths_12[, 1], deaths_nozero[, 5:15] / deaths_nozero$popm, all = deaths_nozero[, 4] / deaths_nozero$popm)
actual_wide <- as.data.table(actual_mort) %>% data.table::dcast(Year ~ agegroups, value.var = colnames(as.data.table(actual_mort))[3:14])

actuals <- actual_wide %>%
  gather("variable", "value", -Year) %>%
  mutate(
    Cause = sub("_[[:digit:]][[:digit:]][[:digit:]]", "", variable),
    Cause = sub("_[[:digit:]][[:digit:]]", "", Cause),
    Cause = sub("_[[:digit:]]", "", Cause),
    variable = sub(".*_", "", variable)
  ) %>%
  spread(variable, value)

actuals$Year <- actuals$Year - 1978


##### difference for single
fnorm(generate_mortality_tensor(as.data.table(fit1)) - generate_mortality_tensor(as.data.table(actuals)))

actuals$total <- rowSums(actuals[,3:9])
fit1s$total <- rowSums(fit1s[,3:9])
fit1$total <- rowSums(fit1[,3:9])
fit2$total <- rowSums(fit2[,3:9])
fit3$total <- rowSums(fit3[,3:9])

simple_rmse <- sqrt(mean((actuals[,12] - fit1s[,12])^2 , na.rm = TRUE))
simple_mae <- sum(abs(actuals[,12] - fit1s[,12]), na.rm = TRUE)/nrow(actuals)

single_rmse <- sqrt(mean((actuals[,12] - fit1[,12])^2 , na.rm = TRUE))
single_mae <- sum(abs(actuals[,12] - fit1[,12]), na.rm = TRUE)/nrow(actuals)

quadratic_rmse <- sqrt(mean((actuals[,12] - fit2[,12])^2 , na.rm = TRUE))
quadratic_mae <- sum(abs(actuals[,12] - fit2[,12]), na.rm = TRUE)/nrow(actuals)

cubic_rmse <- sqrt(mean((actuals[,12] - fit3[,12])^2 , na.rm = TRUE))
cubic_mae <- sum(abs(actuals[,12] - fit3[,12]), na.rm = TRUE)/nrow(actuals)

print(simple_rmse*100)
print(single_rmse*100)
print(quadratic_rmse*100)
print(cubic_rmse*100)

print(simple_mae*100)
print(single_mae*100)
print(quadratic_mae*100)
print(cubic_mae*100)

#fnorm(generate_mortality_tensor(as.data.table(pred2)) - generate_mortality_tensor(as.data.table(actuals)))
##### difference for quadratic
fnorm(generate_mortality_tensor(as.data.table(fit3)) - generate_mortality_tensor(as.data.table(actuals)))



generate_mortality_tensor(as.data.table(pred2))

# Create plots
# prepare actual observations for charts
CoD.label <- c("Other non-cardio", 
               "Rheumatic acute","Rheumatic chronic", "Hypertension", 
               "Hypertensive heart","Hypertensive renal","Hypertensive heart renal", "Myocardial infarction","IHD acute","IHD chronic",
               "Pulmonary", "Other cardio")
# Select causes to include in plot (drop all other non-cardio)
include.causes <- c("Rheumatic acute","Rheumatic chronic", "Hypertension", "Hypertensive heart","Hypertensive renal","Hypertensive heart renal", "Myocardial infarction","IHD acute","IHD chronic","Pulmonary", "Other cardio")

col.cause <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#28E2E5", "#000000")

# actuals
actual_plot <- actuals
actual_plot$years.fit <- rowSums(actual_plot[, 3:11])
actual_plot <- actual_plot[, c(1:2, 12)]
actual_plot$Year <- actual_plot$Year + 1978
colnames(actual_plot) = c("year", "CoD", "proportion")

  actual_plot$CoD[actual_plot$CoD == 'all'] <- 'Other non-cardio'
  actual_plot$CoD[actual_plot$CoD == 'd_hypertension'] <- 'Hypertension'
  actual_plot$CoD[actual_plot$CoD == 'd_hypertensive_heart'] <- 'Hypertensive heart'
  actual_plot$CoD[actual_plot$CoD == 'd_hypertensive_hr'] <- 'Hypertensive heart renal'
  actual_plot$CoD[actual_plot$CoD == 'd_hypertensive_renal'] <- 'Hypertensive renal'
  actual_plot$CoD[actual_plot$CoD == 'd_IHD_acute'] <- 'IHD acute'
  actual_plot$CoD[actual_plot$CoD == 'd_IHD_chronic'] <- 'IHD chronic'
  actual_plot$CoD[actual_plot$CoD == 'd_myoc_inf'] <- 'Myocardial infarction'
  actual_plot$CoD[actual_plot$CoD == 'd_oth_cardio'] <- 'Other cardio'
  actual_plot$CoD[actual_plot$CoD == 'd_pulmonary'] <- 'Pulmonary'
  actual_plot$CoD[actual_plot$CoD == 'd_rheumatic_acute'] <- 'Rheumatic acute'
  actual_plot$CoD[actual_plot$CoD == 'd_rheumatic_chronic'] <- 'Rheumatic chronic'

# simple fit
fit1s_plot <- fit1s
fit1s_plot$years.fit <- rowSums(fit1s_plot[, 3:11])
fit1s_plot <- fit1s_plot[, c(1:2, 12)]
fit1s_plot$Year <- fit1s_plot$Year + 1978
colnames(fit1s_plot) = c("year", "CoD", "proportion")

  fit1s_plot$CoD[fit1s_plot$CoD == 'all'] <- 'Other non-cardio'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_hypertension'] <- 'Hypertension'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_hypertensive_heart'] <- 'Hypertensive heart'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_hypertensive_hr'] <- 'Hypertensive heart renal'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_hypertensive_renal'] <- 'Hypertensive renal'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_IHD_acute'] <- 'IHD acute'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_IHD_chronic'] <- 'IHD chronic'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_myoc_inf'] <- 'Myocardial infarction'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_oth_cardio'] <- 'Other cardio'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_pulmonary'] <- 'Pulmonary'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_rheumatic_acute'] <- 'Rheumatic acute'
  fit1s_plot$CoD[fit1s_plot$CoD == 'd_rheumatic_chronic'] <- 'Rheumatic chronic'

# single fit
fit1_plot <- fit1
fit1_plot$years.fit <- rowSums(fit1_plot[, 3:11])
fit1_plot <- fit1_plot[, c(1:2, 12)]
fit1_plot$Year <- fit1_plot$Year + 1978
colnames(fit1_plot) = c("year", "CoD", "proportion")

  fit1_plot$CoD[fit1_plot$CoD == 'all'] <- 'Other non-cardio'
  fit1_plot$CoD[fit1_plot$CoD == 'd_hypertension'] <- 'Hypertension'
  fit1_plot$CoD[fit1_plot$CoD == 'd_hypertensive_heart'] <- 'Hypertensive heart'
  fit1_plot$CoD[fit1_plot$CoD == 'd_hypertensive_hr'] <- 'Hypertensive heart renal'
  fit1_plot$CoD[fit1_plot$CoD == 'd_hypertensive_renal'] <- 'Hypertensive renal'
  fit1_plot$CoD[fit1_plot$CoD == 'd_IHD_acute'] <- 'IHD acute'
  fit1_plot$CoD[fit1_plot$CoD == 'd_IHD_chronic'] <- 'IHD chronic'
  fit1_plot$CoD[fit1_plot$CoD == 'd_myoc_inf'] <- 'Myocardial infarction'
  fit1_plot$CoD[fit1_plot$CoD == 'd_oth_cardio'] <- 'Other cardio'
  fit1_plot$CoD[fit1_plot$CoD == 'd_pulmonary'] <- 'Pulmonary'
  fit1_plot$CoD[fit1_plot$CoD == 'd_rheumatic_acute'] <- 'Rheumatic acute'
  fit1_plot$CoD[fit1_plot$CoD == 'd_rheumatic_chronic'] <- 'Rheumatic chronic'

# quadratic fit
fit2_plot <- fit2
fit2_plot$years.fit <- rowSums(fit2_plot[, 3:11])
fit2_plot <- fit2_plot[, c(1:2, 12)]
fit2_plot$Year <- fit2_plot$Year + 1978
colnames(fit2_plot) = c("year", "CoD", "proportion")

  fit2_plot$CoD[fit2_plot$CoD == 'all'] <- 'Other non-cardio'
  fit2_plot$CoD[fit2_plot$CoD == 'd_hypertension'] <- 'Hypertension'
  fit2_plot$CoD[fit2_plot$CoD == 'd_hypertensive_heart'] <- 'Hypertensive heart'
  fit2_plot$CoD[fit2_plot$CoD == 'd_hypertensive_hr'] <- 'Hypertensive heart renal'
  fit2_plot$CoD[fit2_plot$CoD == 'd_hypertensive_renal'] <- 'Hypertensive renal'
  fit2_plot$CoD[fit2_plot$CoD == 'd_IHD_acute'] <- 'IHD acute'
  fit2_plot$CoD[fit2_plot$CoD == 'd_IHD_chronic'] <- 'IHD chronic'
  fit2_plot$CoD[fit2_plot$CoD == 'd_myoc_inf'] <- 'Myocardial infarction'
  fit2_plot$CoD[fit2_plot$CoD == 'd_oth_cardio'] <- 'Other cardio'
  fit2_plot$CoD[fit2_plot$CoD == 'd_pulmonary'] <- 'Pulmonary'
  fit2_plot$CoD[fit2_plot$CoD == 'd_rheumatic_acute'] <- 'Rheumatic acute'
  fit2_plot$CoD[fit2_plot$CoD == 'd_rheumatic_chronic'] <- 'Rheumatic chronic'


# cubic fit
fit3_plot <- fit3
fit3_plot$years.fit <- rowSums(fit3_plot[, 3:11])
fit3_plot <- fit3_plot[, c(1:2, 12)]
fit3_plot$Year <- fit3_plot$Year + 1978
colnames(fit3_plot) = c("year", "CoD", "proportion")

  fit3_plot$CoD[fit3_plot$CoD == 'all'] <- 'Other non-cardio'
  fit3_plot$CoD[fit3_plot$CoD == 'd_hypertension'] <- 'Hypertension'
  fit3_plot$CoD[fit3_plot$CoD == 'd_hypertensive_heart'] <- 'Hypertensive heart'
  fit3_plot$CoD[fit3_plot$CoD == 'd_hypertensive_hr'] <- 'Hypertensive heart renal'
  fit3_plot$CoD[fit3_plot$CoD == 'd_hypertensive_renal'] <- 'Hypertensive renal'
  fit3_plot$CoD[fit3_plot$CoD == 'd_IHD_acute'] <- 'IHD acute'
  fit3_plot$CoD[fit3_plot$CoD == 'd_IHD_chronic'] <- 'IHD chronic'
  fit3_plot$CoD[fit3_plot$CoD == 'd_myoc_inf'] <- 'Myocardial infarction'
  fit3_plot$CoD[fit3_plot$CoD == 'd_oth_cardio'] <- 'Other cardio'
  fit3_plot$CoD[fit3_plot$CoD == 'd_pulmonary'] <- 'Pulmonary'
  fit3_plot$CoD[fit3_plot$CoD == 'd_rheumatic_acute'] <- 'Rheumatic acute'
  fit3_plot$CoD[fit3_plot$CoD == 'd_rheumatic_chronic'] <- 'Rheumatic chronic'

#actual plot
g <- ggplot(data = actual_plot, aes(x = CoD, y = proportion, group = year, colour = year))
theme_set(theme_bw())

actual_plot_graph <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  #scale_x_continuous(label = CoD.label)
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Actual") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none") 

#simple fit plot
g <- ggplot(data = fit1s_plot, aes(x = CoD, y = proportion, group = year, colour = year))
theme_set(theme_bw())

fit1s_plot_graph <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Simple Fit (MLR)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none") +
  theme(axis.title.y=element_blank())


#single fit plot
g <- ggplot(data = fit1_plot, aes(x = CoD, y = proportion, group = year, colour = year))
theme_set(theme_bw())

fit1_plot_graph <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Single Fit (MLR)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none") +
  theme(axis.title.y=element_blank())

#quadratic fit plot
g <- ggplot(data = fit2_plot, aes(x = CoD, y = proportion, group = year, colour = year))
theme_set(theme_bw())

fit2_plot_graph <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Quadratic Fit (MLR)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none") +
  theme(axis.title.y=element_blank())

#cubic fit plot
g <- ggplot(data = fit3_plot, aes(x = CoD, y = proportion, group = year, colour = year))
theme_set(theme_bw())

fit3_plot_graph <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Cubic Fit (MLR)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none") +
  theme(axis.title.y=element_blank())

#savepdf("us_MLR_comparison_m", width = 32, height = 12, toplines = 0.8)
par(mar=c(5, 4, 4, 8), xpd=TRUE)
actual_plot_graph+ fit1s_plot_graph + fit1_plot_graph + fit2_plot_graph + fit3_plot_graph + 
  plot_layout(ncol = 5)
#dev.off()


# plot 2: time series by cause comparing alpha and observed
# Create plot

theme_set(theme_bw())

simple_proj_m <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  #geom_line(data = subset(fit1s_plot, CoD %in% include.causes), linetype = "dashed") + 
  #geom_line(data = subset(actual_plot, CoD %in% include.causes)) +  
  geom_line(data = fit1s_plot, linetype = "dashed") + 
  geom_line(data = actual_plot) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2031, 4)) +
  ggtitle("Simple Fit (MLR)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none")

single_proj_m <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  #geom_line(data = subset(fit1s_plot, CoD %in% include.causes), linetype = "dashed") + 
  #geom_line(data = subset(actual_plot, CoD %in% include.causes)) +  
  geom_line(data = fit1_plot, linetype = "dashed") + 
  geom_line(data = actual_plot) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2031, 4)) +
  ggtitle("Single Fit (MLR)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none") +
  theme(axis.title.y=element_blank())

quadratic_proj_m <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  #geom_line(data = subset(fit1s_plot, CoD %in% include.causes), linetype = "dashed") + 
  #geom_line(data = subset(actual_plot, CoD %in% include.causes)) +  
  geom_line(data = fit2_plot, linetype = "dashed") + 
  geom_line(data = actual_plot) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2031, 4)) +
  ggtitle("Quadratic Fit (MLR)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none") +
  theme(axis.title.y=element_blank())

cubic_proj_m <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  #geom_line(data = subset(fit1s_plot, CoD %in% include.causes), linetype = "dashed") + 
  #geom_line(data = subset(actual_plot, CoD %in% include.causes)) +  
  geom_line(data = fit3_plot, linetype = "dashed") + 
  geom_line(data = actual_plot) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2031, 4)) +
  ggtitle("Cubic Fit (MLR)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(axis.title.y=element_blank())


savepdf("us_MLR_all_proj_m", width = 32, height = 12, toplines = 0.8)
par(mar=c(2, 2, 2, 2), xpd=TRUE)
simple_proj_m + single_proj_m + quadratic_proj_m + cubic_proj_m + 
  plot_layout(ncol = 4)
dev.off()