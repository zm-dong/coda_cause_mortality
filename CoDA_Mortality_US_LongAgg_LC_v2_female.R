# This analysis looks at the US mortality by cause data
# Long list of causes (aggregated)
# We apply the alpha-transformation to this data set, where alpha = 0, we get the isometric log ratios

rm(list=ls())
#load packages
setwd("~/Documents/PhD/03_CODA/Model_US")
source("function.r")
source("save_function.r")
source("alfainv_mod.r")
source("Interval_Forecasts_alpha_v0.1.R")
source("Interval_Forecasts_CLR_v0.1.R")
source("Interval_Forecasts_ILR_v0.1.R")

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
library(rcartocolor)
library(paletteer)

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

# setting up arrays
datiM1 <- array(c(datiM0), dim=c(k,n,m))
  dimnames(datiM1) <- list(CoD, years, ages)  

# Different from the original, where the array is arranged as (ages, causes, years) rather than (ages, years, causes)
datiF1 <- array(c(datiF0), dim=c(k,n,m))
  dimnames(datiF1) <- list(CoD, years, ages)
  
datiT1 <- array(c(datiT0), dim=c(k,n,m))
  dimnames(datiT1) <- list(CoD, years, ages)

#Kjaergaard 2019 only selected ages 25+ and for 1955 onwards
#For UK mortality by cause analysis, select all ages and years to test alpha-transformation
datiM<- datiM1
datiF<- datiF1

# ================================================================
# read national life table and prepare data 
# ================================================================
LT.F1 <-(read.table("fltper_5x1_USA.txt",skip=1,header=TRUE,stringsAsFactors=F))
LT.M1 <-(read.table("mltper_5x1_USA.txt",skip=1,header=TRUE,stringsAsFactors=F))

#Create a table with only deaths (dx) to get life table deaths for each year
F.dx.nat <- matrix(LT.F1$dx,length(unique(LT.F1$Age)),length(unique(LT.F1$Year)))
colnames(F.dx.nat) <- unique(LT.F1$Year)
rownames(F.dx.nat) <- unique(LT.F1$Age)

# Repeat life table deaths for males
M.dx.nat <- matrix(LT.M1$dx,length(unique(LT.M1$Age)),length(unique(LT.M1$Year)))
colnames(M.dx.nat) <- unique(LT.M1$Year)
rownames(M.dx.nat) <- unique(LT.M1$Age)

#Create a table with only exposure (Lx) to get life table exposure for each year 
#useful for Lee-Carter Comparisons later
F.Lx.nat <- matrix(LT.F1$Lx,length(unique(LT.F1$Age)),length(unique(LT.F1$Year)))
colnames(F.Lx.nat) <- unique(LT.F1$Year)
M.Lx.nat <- matrix(LT.M1$Lx,length(unique(LT.M1$Age)),length(unique(LT.M1$Year)))
colnames(M.Lx.nat) <- unique(LT.M1$Year)

# Take a subset of the years based on what is available by cause
M.dx.nat1 <- M.dx.nat[,colnames(M.dx.nat) %in% years]
M.Lx.nat1 <- M.Lx.nat[,colnames(M.Lx.nat) %in% years]
F.dx.nat1 <- F.dx.nat[,colnames(F.dx.nat) %in% years]
F.Lx.nat1 <- F.Lx.nat[,colnames(F.Lx.nat) %in% years]

rownames(F.dx.nat1) <- unique(LT.F1$Age)
rownames(F.Lx.nat1) <- unique(LT.F1$Age)
rownames(M.dx.nat1) <- unique(LT.M1$Age)
rownames(M.Lx.nat1) <- unique(LT.M1$Age)

age.temp<- unique(LT.M1$Age)

M.dx.nat2 <- M.dx.nat1
M.Lx.nat2 <- M.Lx.nat1
F.dx.nat2 <- F.dx.nat1
F.Lx.nat2 <- F.Lx.nat1

# set the used ages from the life tables
# Kjaergaard 2019 selected ages 25+
# USA data has aggregated age bands from age 100 onwards (i.e. last 3 age bands)

M.dx.nat3 <- rbind(M.dx.nat2[1,]+M.dx.nat2[2,]+M.dx.nat2[3,]+M.dx.nat2[4,]+M.dx.nat2[5,]+M.dx.nat2[6,]+M.dx.nat2[7,],
                   M.dx.nat2[8,]+M.dx.nat2[9,],
                   M.dx.nat2[10,]+M.dx.nat2[11,],
                   M.dx.nat2[12,]+M.dx.nat2[13,],
                   M.dx.nat2[14,]+M.dx.nat2[15,],
                   M.dx.nat2[16,]+M.dx.nat2[17,],
                   M.dx.nat2[18,]+M.dx.nat2[19,],
                   M.dx.nat2[20,]+M.dx.nat2[21,],
                   M.dx.nat2[22,]+M.dx.nat2[23,]+M.dx.nat2[24,])


F.dx.nat3 <- rbind(F.dx.nat2[1,]+F.dx.nat2[2,]+F.dx.nat2[3,]+F.dx.nat2[4,]+F.dx.nat2[5,]+F.dx.nat2[6,]+F.dx.nat2[7,],
                   F.dx.nat2[8,]+F.dx.nat2[9,],
                   F.dx.nat2[10,]+F.dx.nat2[11,],
                   F.dx.nat2[12,]+F.dx.nat2[13,],
                   F.dx.nat2[14,]+F.dx.nat2[15,],
                   F.dx.nat2[16,]+F.dx.nat2[17,],
                   F.dx.nat2[18,]+F.dx.nat2[19,],
                   F.dx.nat2[20,]+F.dx.nat2[21,],
                   F.dx.nat2[22,]+F.dx.nat2[23,]+F.dx.nat2[24,])

M.Lx.nat3 <- rbind(M.Lx.nat2[1,]+M.Lx.nat2[2,]+M.Lx.nat2[3,]+M.Lx.nat2[4,]+M.Lx.nat2[5,]+M.Lx.nat2[6,]+M.Lx.nat2[7,],
                   M.Lx.nat2[8,]+M.Lx.nat2[9,],
                   M.Lx.nat2[10,]+M.Lx.nat2[11,],
                   M.Lx.nat2[12,]+M.Lx.nat2[13,],
                   M.Lx.nat2[14,]+M.Lx.nat2[15,],
                   M.Lx.nat2[16,]+M.Lx.nat2[17,],
                   M.Lx.nat2[18,]+M.Lx.nat2[19,],
                   M.Lx.nat2[20,]+M.Lx.nat2[21,],
                   M.Lx.nat2[22,]+M.Lx.nat2[23,]+M.Lx.nat2[24,])

F.Lx.nat3 <- rbind(F.Lx.nat2[1,]+F.Lx.nat2[2,]+F.Lx.nat2[3,]+F.Lx.nat2[4,]+F.Lx.nat2[5,]+F.Lx.nat2[6,]+F.Lx.nat2[7,],
                   F.Lx.nat2[8,]+F.Lx.nat2[9,],
                   F.Lx.nat2[10,]+F.Lx.nat2[11,],
                   F.Lx.nat2[12,]+F.Lx.nat2[13,],
                   F.Lx.nat2[14,]+F.Lx.nat2[15,],
                   F.Lx.nat2[16,]+F.Lx.nat2[17,],
                   F.Lx.nat2[18,]+F.Lx.nat2[19,],
                   F.Lx.nat2[20,]+F.Lx.nat2[21,],
                   F.Lx.nat2[22,]+F.Lx.nat2[23,]+F.Lx.nat2[24,])

# define total deaths for each age band and year
# create a table which has total deaths - only for by cause analysis
# sum across all causes (cause, years, age bands)
datiM.total <- colSums(datiM, dims = 1)
datiF.total <- colSums(datiF, dims = 1)

# ================================================================
# prepare arrays 
# ================================================================

# the pre-transposed state is (causes, years, age bands)
# transpose the array such that it aligns with the stacked matrix configuration (years, ages, causes)
# Each row should be the deaths over a single year for a selected age band and cause combination  
datiF.rela1_0 <- datiM.rela1_0 <- datiF.rela1.5 <- datiM.rela1.5 <- 
  datiF.rela3 <- datiM.rela3 <- 
  datiF.rela2 <- datiM.rela2<- 
  datiF.rela1 <- datiM.rela1 <- 
  array(NA,dim=c(nrow(datiM[1,,]),ncol(datiM[1,,]),12))
# since it is fixed by cause:
# where there are 10 causes, so nrow(datiM[1,,]) is the number of years for cause 1
# ncol(datiM[1,,]) is the number of age bands for cause 1

# dimensions are (number of years = 43, number of age bands = 9, number of causes of death = 12)
m.M.cause.i_0 <- m.F.cause.i_0 <- 
  m.M.cause.i <- m.F.cause.i <- 
  datiM.temp <- datiF.temp <- 
  array(NA,dim=c(nrow(datiM[1,,]), ncol(datiM[1,,]),12))

# calculate life table deaths without zeroes (i.e. replaced by a small dummy number)
dummy <- 0.25 # this was 0.25 originally, to remove the zeroes

# for loop is run for the number of causes (12)
   for(i in 1:12){
     datiM.temp[,,i] <- ifelse(datiM[i,,]==0, dummy, datiM[i,,])
     datiF.temp[,,i] <- ifelse(datiF[i,,]==0, dummy, datiF[i,,]) #dealing with zeroes in the data
     
     datiM.rela1[,,i] <- t(M.dx.nat3)*(datiM[i,,]/datiM.total)
     datiF.rela1[,,i] <- t(F.dx.nat3)*(datiF[i,,]/datiF.total)
     
     datiM.rela1_0[,,i] <- t(M.dx.nat3)*(datiM.temp[,,i]/datiM.total)
     datiF.rela1_0[,,i] <- t(F.dx.nat3)*(datiF.temp[,,i]/datiF.total)
     
     m.M.cause.i[,,i] <- (datiM.rela1_0[,,i]) / t(M.Lx.nat3)  
     m.F.cause.i[,,i] <- (datiF.rela1_0[,,i]) / t(F.Lx.nat3)  
   }

# ================================================================
# life table deaths for females
# ================================================================
dx.com.F <- cbind(datiF.rela1[,,1],datiF.rela1[,,2],datiF.rela1[,,3],datiF.rela1[,,4]
                  ,datiF.rela1[,,5],datiF.rela1[,,6],datiF.rela1[,,7],datiF.rela1[,,8]
                  ,datiF.rela1[,,9],datiF.rela1[,,10]
                  ,datiF.rela1[,,11],datiF.rela1[,,12])#,datiF.rela1[,,13],datiF.rela1[,,14],datiF.rela1[,,15],datiF.rela1[,,16])

#zeroes are replaced in this data:
dx.com.F_0 <- cbind(datiF.rela1_0[,,1],datiF.rela1_0[,,2],datiF.rela1_0[,,3],datiF.rela1_0[,,4]
                  ,datiF.rela1_0[,,5],datiF.rela1_0[,,6],datiF.rela1_0[,,7],datiF.rela1_0[,,8]
                  ,datiF.rela1_0[,,9],datiF.rela1_0[,,10]
                  ,datiF.rela1_0[,,11],datiF.rela1_0[,,12]) 

# Calculate average deaths to be used as weights 
# averaging across causes, the last variable should be equal to the number of causes (12)
dataF.colm <- dataF.colm_0 <- array(NA,dim=c(1,length(ages),12))
dataM.colm <- dataM.colm_0 <- array(NA,dim=c(1,length(ages),12))

dea.F.total <- datiF.rela1[,,1]+datiF.rela1[,,2]+datiF.rela1[,,3]+datiF.rela1[,,4]+datiF.rela1[,,5]+datiF.rela1[,,6]+datiF.rela1[,,7]+datiF.rela1[,,8]+datiF.rela1[,,9]+datiF.rela1[,,10]+datiF.rela1[,,11]+datiF.rela1[,,12]#+datiF.rela1[,,13]+datiF.rela1[,,14]+datiF.rela1[,,15]+datiF.rela1[,,16]
dea.M.total <- datiM.rela1[,,1]+datiM.rela1[,,2]+datiM.rela1[,,3]+datiM.rela1[,,4]+datiM.rela1[,,5]+datiM.rela1[,,6]+datiM.rela1[,,7]+datiM.rela1[,,8]+datiM.rela1[,,9]+datiM.rela1[,,10]+datiM.rela1[,,11]+datiM.rela1[,,12]#+datiM.rela1[,,13]+datiM.rela1[,,14]+datiM.rela1[,,15]+datiM.rela1[,,16]

dati.F.mean <- colMeans(dea.F.total)
dati.M.mean <- colMeans(dea.M.total)

# for loop run for the number of causes (12)
for(i in 1:12){
  dataF.colm[,,i] <- colMeans(datiF.rela1[,,i])/sum(dati.F.mean)  
  dataM.colm[,,i] <- colMeans(datiM.rela1[,,i])/sum(dati.M.mean)  
}

dataF.colm1 <- c(dataF.colm[,,1],dataF.colm[,,2],dataF.colm[,,3],dataF.colm[,,4],dataF.colm[,,5],dataF.colm[,,6],dataF.colm[,,7],dataF.colm[,,8]
                 ,dataF.colm[,,9],dataF.colm[,,10]
                 ,dataF.colm[,,11],dataF.colm[,,12])#,dataF.colm[,,13],dataF.colm[,,14],dataF.colm[,,15],dataF.colm[,,16])

# weight plot array dimensions should be 1, number of age bands (9), number of causes (12)
weight.plot <- array(dataF.colm1,dim=c(1,9,12))

dea.rela.dist <-  dati.F.mean/sum(dati.F.mean)
dea.rela.dist1 <- rep(dea.rela.dist,times=length(CoD))

# ================================================================
# cross validation [rerun this for US data]
# ================================================================
# cross-validation
# define years used for the analysis 
# for US and UK mortality by cause, use all years

# training data: years 2001 to 2008 (US: years 1979 to 2001)
   input_k <- 33 # years of training data (UK: 8-11 for optimising alpha, and 12 for projecting; US: 23-32 for optimising alpha, 33 for projecting)
   start_k <- 1 
   years.fit.1 <- years[start_k:input_k] 
   
# validation data: years 2009 to 2012 (US: 2002 to 2011)
# testing data: years 2013 to 2016 (US: 2012 to 2021)
ih = 100 #forecast window is 10 years for validation, and 20 years for testing. Forecast window is 30 years for projection.
# rerun for revision on 100 years to include in response

# to test the methodology, a compositional approach to forecasting lifetimes is used
# based on the singular value decomposition (SVD)
years.for <- c(max(years.fit.1)+1:ih)
# set the data to the correct set: dx.com.F_0 for no zeroes, dx.com.F for with zeroes
dx <- as.matrix(dx.com.F)
#dx <- as.matrix(dx.com.F_0)

  dx <- dx / rowSums(dx)
  close.dx <- dx[start_k:input_k,]
  
  rownames(close.dx) <- as.character(years.fit.1)
  gx <- geometricmeanCol(close.dx)

# ================================================================
# Fit the models
# ================================================================

# forecasting total life table deaths by age and gender (F.dx.nat3)
F.dx.nat3_fore = matrix(NA, 9, ih)
for (j in 1:9)
{
  for(ik in 1:ih)
  {
    F.dx.nat3_fore[j,ik] = exp(forecast(auto.arima(log(F.dx.nat3[j,1:(ik)])), h = 1)$mean)
  }
}
# forecast life table deaths (datiF.rela1) and then divide the result by exposure to get mortality rate
colnames(F.dx.nat3_fore) <- as.character(years.for) #all causes, by age and year
rownames(F.dx.nat3_fore) <- as.character(ages)

dx.cent <- Rfast::eachrow(close.dx, gx, oper = "/")

# prepare actual observations for charts

CoD.label <- c("Other non-cardio", 
               "Rheumatic acute","Rheumatic chronic", "Hypertension", 
               "Hypertensive heart","Hypertensive renal","Hypertensive heart renal", "Myocardial infarction","IHD acute","IHD chronic",
               "Pulmonary", "Other cardio")

# Select causes to include in plot (drop all other non-cardio)
include.causes <- c("Rheumatic acute","Rheumatic chronic", "Hypertension", "Hypertensive heart","Hypertensive renal","Hypertensive heart renal", "Myocardial infarction","IHD acute","IHD chronic","Pulmonary", "Other cardio")

model.fit.for.lc <- COD.LC(mx=m.F.cause.i, t=ih, k=length(CoD),ages=ages,years=years)
pro.dea.lc <- t(apply(model.fit.for.lc$dx.forcast, MARGIN=c(1,3),FUN=sum)/100000)
years.lc <- c(min(years):2121)
#years.lc <- c(1979:2051)
colnames(pro.dea.lc) = years.lc
rownames(pro.dea.lc) = CoD

plot_lc <- pro.dea.lc
colnames(plot_lc) = years.lc
rownames(plot_lc) = CoD.label

datLong_lc <- melt(data          = plot_lc,
                   id.vars       = c(CoD.label),
                   measure.vars  = c(years.lc),
                   variable.name = "year",
                   value.name    = "proportion")
colnames(datLong_lc) = c("CoD", "year", "proportion")

#---

obs.ar <- array(dx,dim=c(43,9,12))
obs.ar.mx <- array(m.F.cause.i,dim=c(43,9,12))

pro.dea.obs <- t(apply(obs.ar, MARGIN=c(1,3),FUN=sum))
pro.dea.obs.mx <- t(apply(obs.ar.mx, MARGIN=c(1,3),FUN=sum))

colnames(pro.dea.obs) = years
rownames(pro.dea.obs) = CoD

colnames(pro.dea.obs.mx) = years
rownames(pro.dea.obs.mx) = CoD

total_years <- c(years.fit.1, years.for)

plot_data_obs <- pro.dea.obs

colnames(plot_data_obs) = years
rownames(plot_data_obs) = CoD.label

datLong_obs <- melt(data          = plot_data_obs,
                id.vars       = c(CoD.label),
                measure.vars  = c(years),
                variable.name = "year",
                value.name    = "proportion")
colnames(datLong_obs) = c("CoD", "year", "proportion")

# ================================================================
# Fit the models: CLR
# This section replicates the CT-CODA part from above
# It does not leverage the functions introduced by Kjaergaard et al. (2019)
# ================================================================

# This was the original code used by Kjaergaard et al. (2019). 
# The CLR was applied to centred data, rather than the original composition
# This is because cause-specific death matrices are stacked horizontally to form a T * NK giant matrix
# Every row still sums to the life table radix because it represents the number of deaths in year t
#dx.cent <- Compositional::perturbation(close.dx, gx, oper = "-")

clr.cent <- clr(dx.cent)
rownames(clr.cent) <- as.character(years.fit.1)
  
#SVD on CLR
par <- svd(clr.cent, nu=3, nv=3)
U <- par$u 
V <- par$v
S <- diag(par$d)
  
bx<- V[,1]
kt<- S[1,1]*U[,1] 
bx2 <- V[,2]
kt2 <- S[2,2]*U[,2]
bx3 <- V[,3]
kt3 <- S[3,3]*U[,3]
  
R2.one <- cbind(par$d[1]^2/sum(par$d^2), par$d[2]^2/sum(par$d^2),par$d[3]^2/sum(par$d^2),par$d[4]^2/sum(par$d^2)) 
  
# fit Arima model to kt
order.one=c(1,1,1)
fcst.one <- forecast(auto.arima(kt,max.d=1),h=ih)
kt.fit.one <- fcst.one$fitted
kt.for.one <- fcst.one$mean

order.two=c(1,0,1)
fcst.two <- forecast(auto.arima(kt2,max.d=1),h=ih)
kt.fit.two <- fcst.two$fitted
kt.for.two <- fcst.two$mean

order.3=c(1,0,1)
fcst.3 <- forecast(auto.arima(kt3,max.d=1),h=ih)
kt.fit.3 <- fcst.3$fitted
kt.for.3 <- fcst.3$mean

gx.all <- matrix(gx,length(ages),k)
bx.all <- matrix(bx,length(ages),k)
bx2.all <- matrix(bx2,length(ages),k)

#projections (CLR)
clr.proj.fit.one <- matrix(kt,length(years.fit.1),1) %*% t(bx) + matrix(kt2,length(years.fit.1),1) %*% t(bx2) + matrix(kt3,length(years.fit.1),1) %*% t(bx3)
clr.proj.for.one <- matrix(c(kt, kt.for.one),length(years.fit.1)+ih,1) %*% t(bx) + matrix(c(kt2, kt.for.two),length(years.fit.1)+ih,1) %*% t(bx2) + matrix(c(kt3, kt.for.3),length(years.fit.1)+ih,1) %*% t(bx3)

#Inv CLR
BK.clr.proj.fit.one <- clrInv(clr.proj.fit.one)
BK.clr.proj.for.one <- clrInv(clr.proj.for.one)

#Add back geometric mean
clr.proj.fit.one.gx.per <- Rfast::eachrow(BK.clr.proj.fit.one, gx, oper = "*")
clr.proj.for.one.gx.per <- Rfast::eachrow(BK.clr.proj.for.one, gx, oper = "*")

clr.proj.fit.one.gx <- clr.proj.fit.one.gx.per / rowSums(clr.proj.fit.one.gx.per)
clr.proj.for.one.gx <- clr.proj.for.one.gx.per / rowSums(clr.proj.for.one.gx.per)

clr.proj.fit.one.ar <- array((unclass(clr.proj.fit.one.gx)), dim=c(length(years.fit.1),m,length(CoD)),dimnames=list(as.character(years.fit.1[1]:(years.fit.1[length(years.fit.1)])),as.character(ages),CoD))
clr.proj.for.one.ar <- array((unclass(clr.proj.for.one.gx)), dim=c(length(years.fit.1)+ih,m,length(CoD)), dimnames=list(as.character(years.fit.1[1]:(years.for[length(years.for)])), as.character(ages),CoD))

#The array of deaths (compositional)  
dx.ar <- array(dx.com.F, dim=c(length(years.fit.1),m,length(CoD)),dimnames=list(as.character(years.fit.1[1]:(years.fit.1[length(years.fit.1)])),as.character(ages),CoD))
dx.ar1 <- array(acomp(dx.com.F), dim=c(length(years.for),m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages),CoD))
  
clr.for.mx.one <- array(NA, dim=c(length(years.for),m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages)))
  for.dx.one <- array(NA, dim=c(length(years.for),m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages)))
  fit.dx.one <- array(NA, dim=c(length(years.fit.1),m,length(CoD)),dimnames=list(as.character(years.fit.1[1]:(years.fit.1[length(years.fit.1)])),as.character(ages)))

clr.proj.for.one.tot <- clr.proj.for.one.ar[,,1]
for(i in 2:k){clr.proj.for.one.tot <- clr.proj.for.one.tot + clr.proj.for.one.ar[,,i]}

  #clr.proj.tot <- array(NA,dim=c(n,m,k),dimnames=list(as.character(years),as.character(ages)))
  #for(y in 1:k){
  #  for(x in 1:j){
  #    clr.proj.tot <- cbind(clr.proj.fit.one.ar[,x,y], clr.proj.for.one.ar[,x,y])
  #  }
  #}
  
  # prepare forecast for charts
  
  pro.dea.clr.LC <- t(apply(clr.proj.for.one.ar, MARGIN=c(1,3),FUN=sum)) 
  colnames(pro.dea.clr.LC) = total_years
  rownames(pro.dea.clr.LC) = CoD
  
  plot_data_clr <- pro.dea.clr.LC[,1:43]
  
  colnames(plot_data_clr) = years
  rownames(plot_data_clr) = CoD.label
  
  datLong_clr <- melt(data          = plot_data_clr,
                      id.vars       = c(CoD.label),
                      measure.vars  = c(years),
                      variable.name = "year",
                      value.name    = "proportion")
  colnames(datLong_clr) = c("CoD", "year", "proportion")
  
  # CLR Data for plot 2
  plot2_data_clr <- pro.dea.clr.LC
  colnames(plot2_data_clr) = total_years
  rownames(plot2_data_clr) = CoD.label
  datLong2_clr <- melt(data          = plot2_data_clr,
                       id.vars       = c(CoD.label),
                       measure.vars  = c(total_years),
                       variable.name = "year",
                       value.name    = "proportion")
  colnames(datLong2_clr) = c("CoD", "year", "proportion")
  
  # ========================================================
  # interval forecast for CLR
  # ========================================================
  
  #constructing interval forecast
  interval_forecast_clr_full <- CLR_fun_int(object = dx.cent, ncomp_method = "fixed", ncomp_tuning = 0.001, fh = ih,
                                            fore_method = "rwf", B = 399, K = 1, sig = 20)
  
  interval_forecast_lower_clr_gx_per <- perturbation(interval_forecast_clr_full[1,,], gx, oper = "*")
  interval_forecast_upper_clr_gx_per <- perturbation(interval_forecast_clr_full[2,,], gx, oper = "*")
  
  interval_forecast_lower_clr_ar <- array((unclass(interval_forecast_lower_clr_gx_per)), dim=c(ih,m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages),CoD))
  interval_forecast_upper_clr_ar <- array((unclass(interval_forecast_upper_clr_gx_per)), dim=c(ih,m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages),CoD))
  
  # and interval forecasts
  interval.lower.dea.CLR.LC <- t(apply(interval_forecast_lower_clr_ar, MARGIN=c(1,3),FUN=sum, na.rm = TRUE))
  interval.upper.dea.CLR.LC <- t(apply(interval_forecast_upper_clr_ar, MARGIN=c(1,3),FUN=sum, na.rm = TRUE))
  
  colnames(interval.lower.dea.CLR.LC) = years.for
  rownames(interval.lower.dea.CLR.LC) = CoD
  colnames(interval.upper.dea.CLR.LC) = years.for
  rownames(interval.upper.dea.CLR.LC) = CoD
  
  # plot interval forecasts
  # lower bound
  interval_lower_CLR <- interval.lower.dea.CLR.LC
  colnames(interval_lower_CLR) = years.for
  rownames(interval_lower_CLR) = CoD.label
  datLong_interval_lower_CLR <- melt(data          = interval_lower_CLR,
                                     id.vars       = c(CoD.label),
                                     measure.vars  = c(years.for),
                                     variable.name = "year",
                                     value.name    = "proportion")
  colnames(datLong_interval_lower_CLR) = c("CoD", "year", "proportion")
  
  # upper bound
  interval_upper_CLR <- interval.upper.dea.CLR.LC
  colnames(interval_upper_CLR) = years.for
  rownames(interval_upper_CLR) = CoD.label
  datLong_interval_upper_CLR <- melt(data          = interval_upper_CLR,
                                     id.vars       = c(CoD.label),
                                     measure.vars  = c(years.for),
                                     variable.name = "year",
                                     value.name    = "proportion")
  colnames(datLong_interval_upper_CLR) = c("CoD", "year", "proportion")
  
  cpd(datLong_obs[397:516,], datLong_interval_lower_CLR[1:120,], datLong_interval_upper_CLR[1:120,], 90)
  interval_score(datLong_obs[397:516,], datLong_interval_lower_CLR[1:120,], datLong_interval_upper_CLR[1:120,], 90)
  
  # create plot clr
  interval_CLR_f <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
    geom_line(data = subset(datLong_interval_lower_CLR, CoD %in% include.causes), linetype = "dashed") + 
    geom_line(data = subset(datLong_interval_upper_CLR, CoD %in% include.causes), linetype = "dashed") +
    geom_line(data = subset(datLong2_clr, CoD %in% include.causes), linetype = "solid") +
    geom_line(data = subset(datLong_obs, CoD %in% include.causes), linetype = "dotted") +  
    labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
    scale_y_continuous(trans='log10') +
    #scale_color_brewer(palette = "Paired") +
    scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
    theme(axis.text = element_text(size = 9),
          axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
    scale_x_continuous("Year", breaks = seq(1979, 2051, 4)) +
    #ggtitle("CLR") + 
    theme(plot.title = element_text(face = "bold",
                                    margin = margin(10, 0, 10, 0),
                                    size = 14)) +
    theme(legend.position = "none") #remove the legend to combine CLR, ILR, and alpha in one chart
  #interval_CLR_f
  
# ================================================================
# Fit the models: ILR
# ================================================================
ilr.cent <- ilr(dx.cent)
rownames(ilr.cent) <- as.character(years.fit.1)

#SVD on ILR
par <- svd(ilr.cent, nu=3, nv=3)
U <- par$u 
V <- par$v
S <- diag(par$d)
  
bx<- V[,1]
kt<- S[1,1]*U[,1] 
bx2 <- V[,2]
kt2 <- S[2,2]*U[,2]
bx3 <- V[,3]
kt3 <- S[3,3]*U[,3]
  
R2.one <- cbind(par$d[1]^2/sum(par$d^2), par$d[2]^2/sum(par$d^2),par$d[3]^2/sum(par$d^2),par$d[4]^2/sum(par$d^2)) 
  
# fit Arima model to kt
order.one=c(1,1,1)
fcst.one <- forecast(auto.arima(kt,max.d=1),h=ih)
kt.fit.one <- fcst.one$fitted
kt.for.one <- fcst.one$mean

order.two=c(1,0,1)
fcst.two <- forecast(auto.arima(kt2,max.d=1),h=ih)
kt.fit.two <- fcst.two$fitted
kt.for.two <- fcst.two$mean
  
order.3=c(1,0,1)
fcst.3 <- forecast(auto.arima(kt3,max.d=1),h=ih)
kt.fit.3 <- fcst.3$fitted
kt.for.3 <- fcst.3$mean
  
gx.all <- matrix(gx,length(ages),k)
bx.all <- matrix(bx,length(ages),k)
bx2.all <- matrix(bx2,length(ages),k)

#projections (ILR)
ilr.proj.fit.one <- matrix(kt,length(years.fit.1),1) %*% t(bx) + matrix(kt2,length(years.fit.1),1) %*% t(bx2) + matrix(kt3,length(years.fit.1),1) %*% t(bx3)
ilr.proj.for.one <- matrix(c(kt, kt.for.one),length(years.fit.1)+ih,1) %*% t(bx) + matrix(c(kt2, kt.for.two),length(years.fit.1)+ih,1) %*% t(bx2) + matrix(c(kt3, kt.for.3),length(years.fit.1)+ih,1) %*% t(bx3)
  
#Inv ILR
BK.ilr.proj.fit.one <- ilrInv(ilr.proj.fit.one)
BK.ilr.proj.for.one <- ilrInv(ilr.proj.for.one)
  
#Add back geometric mean
ilr.proj.fit.one.gx.per <- Rfast::eachrow(BK.ilr.proj.fit.one, gx, oper = "*")
ilr.proj.for.one.gx.per <- Rfast::eachrow(BK.ilr.proj.for.one, gx, oper = "*")

ilr.proj.fit.one.gx <- ilr.proj.fit.one.gx.per / rowSums(ilr.proj.fit.one.gx.per)
ilr.proj.for.one.gx <- ilr.proj.for.one.gx.per / rowSums(ilr.proj.for.one.gx.per)

ilr.proj.fit.one.ar <- array((unclass(ilr.proj.fit.one.gx)), dim=c(length(years.fit.1),m,length(CoD)),dimnames=list(as.character(years.fit.1[1]:(years.fit.1[length(years.fit.1)])),as.character(ages),CoD))
ilr.proj.for.one.ar <- array((unclass(ilr.proj.for.one.gx)), dim=c(length(years.fit.1)+ih,m,length(CoD)), dimnames=list(as.character(years.fit.1[1]:(years.for[length(years.for)])), as.character(ages),CoD))
  
#The array of deaths (compositional)  
dx.ar <- array(dx.com.F, dim=c(length(years.fit.1),m,length(CoD)),dimnames=list(as.character(years.fit.1[1]:(years.fit.1[length(years.fit.1)])),as.character(ages),CoD))
dx.ar1 <- array(acomp(dx.com.F), dim=c(length(years.for),m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages),CoD))
  
ilr.for.mx.one <- array(NA, dim=c(length(years.for),m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages)))
for.dx.one <- array(NA, dim=c(length(years.for),m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages)))
fit.dx.one <- array(NA, dim=c(length(years.fit.1),m,length(CoD)),dimnames=list(as.character(years.fit.1[1]:(years.fit.1[length(years.fit.1)])),as.character(ages)))
  
ilr.proj.for.one.tot <- ilr.proj.for.one.ar[,,1]
for(i in 2:k){ilr.proj.for.one.tot <- ilr.proj.for.one.tot + ilr.proj.for.one.ar[,,i]}
  
#ilr.proj.tot <- array(NA,dim=c(n,m,k),dimnames=list(as.character(years),as.character(ages)))
#for(y in 1:k){
#  for(x in 1:j){
#      ilr.proj.tot <- cbind(ilr.proj.fit.one.ar[,x,y], ilr.proj.for.one.ar[,x,y])
#    }
#  }  
  
pro.dea.ilr.LC <- t(apply(ilr.proj.for.one.ar, MARGIN=c(1,3),FUN=sum)) 
colnames(pro.dea.ilr.LC) = total_years
rownames(pro.dea.ilr.LC) = CoD


# prepare forecasts for charts
plot_data_ilr <- pro.dea.ilr.LC[,1:43]

colnames(plot_data_ilr) = years
rownames(plot_data_ilr) = CoD.label

datLong_ilr <- melt(data          = plot_data_ilr,
                    id.vars       = c(CoD.label),
                    measure.vars  = c(years),
                    variable.name = "year",
                    value.name    = "proportion")
#datLong_ilr
colnames(datLong_ilr) = c("CoD", "year", "proportion")

# ILR Data for plot 2
plot2_data_ilr <- pro.dea.ilr.LC
colnames(plot2_data_ilr) = total_years
rownames(plot2_data_ilr) = CoD.label
datLong2_ilr <- melt(data          = plot2_data_ilr,
                     id.vars       = c(CoD.label),
                     measure.vars  = c(total_years),
                     variable.name = "year",
                     value.name    = "proportion")
colnames(datLong2_ilr) = c("CoD", "year", "proportion")

# ========================================================
# interval forecast for ILR
# ========================================================

#constructing interval forecast
interval_forecast_ilr_full <- ILR_fun_int(object = dx.cent, ncomp_method = "fixed", ncomp_tuning = 0.001, fh = ih,
                                          fore_method = "rwf", B = 399, K = 1, sig = 20)

interval_forecast_lower_ilr_gx_per <- perturbation(interval_forecast_ilr_full[1,,], gx, oper = "*")
interval_forecast_upper_ilr_gx_per <- perturbation(interval_forecast_ilr_full[2,,], gx, oper = "*")

interval_forecast_lower_ilr_ar <- array((unclass(interval_forecast_lower_ilr_gx_per)), dim=c(ih,m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages),CoD))
interval_forecast_upper_ilr_ar <- array((unclass(interval_forecast_upper_ilr_gx_per)), dim=c(ih,m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages),CoD))

# and interval forecasts
interval.lower.dea.ILR.LC <- t(apply(interval_forecast_lower_ilr_ar, MARGIN=c(1,3),FUN=sum, na.rm = TRUE))
interval.upper.dea.ILR.LC <- t(apply(interval_forecast_upper_ilr_ar, MARGIN=c(1,3),FUN=sum, na.rm = TRUE))

colnames(interval.lower.dea.ILR.LC) = years.for
rownames(interval.lower.dea.ILR.LC) = CoD
colnames(interval.upper.dea.ILR.LC) = years.for
rownames(interval.upper.dea.ILR.LC) = CoD

# plot interval forecasts
# lower bound
interval_lower_ILR <- interval.lower.dea.ILR.LC
colnames(interval_lower_ILR) = years.for
rownames(interval_lower_ILR) = CoD.label
datLong_interval_lower_ILR <- melt(data          = interval_lower_ILR,
                                   id.vars       = c(CoD.label),
                                   measure.vars  = c(years.for),
                                   variable.name = "year",
                                   value.name    = "proportion")
colnames(datLong_interval_lower_ILR) = c("CoD", "year", "proportion")

# upper bound
interval_upper_ILR <- interval.upper.dea.ILR.LC
colnames(interval_upper_ILR) = years.for
rownames(interval_upper_ILR) = CoD.label
datLong_interval_upper_ILR <- melt(data          = interval_upper_ILR,
                                   id.vars       = c(CoD.label),
                                   measure.vars  = c(years.for),
                                   variable.name = "year",
                                   value.name    = "proportion")
colnames(datLong_interval_upper_ILR) = c("CoD", "year", "proportion")

#cpd(datLong_obs[397:516,], datLong_interval_lower_ILR[1:120,], datLong_interval_upper_ILR[1:120,], 90)
#interval_score(datLong_obs[397:516,], datLong_interval_lower_ILR[1:120,], datLong_interval_upper_ILR[1:120,], 90)

# create plot ILR
interval_ILR_f <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong_interval_lower_ILR, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_interval_upper_ILR, CoD %in% include.causes), linetype = "dashed") +
  geom_line(data = subset(datLong2_ilr, CoD %in% include.causes), linetype = "solid") +
  geom_line(data = subset(datLong_obs, CoD %in% include.causes), linetype = "dotted") +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2051, 4)) +
  #ggtitle("ILR") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none") + #remove the legend to combine CLR, ILR, and alpha in one chart
  theme(axis.title.y=element_blank()) 
#interval_ILR_f

# ================================================================
# Fit the models: alpha
# applying the same model to alpha-transformed data
# ================================================================
# CT-CODA
# nr.rank = 3 = nu = nv

# rows = years by agebands
# columns = causes
options(rgl.useNULL=TRUE)
library(rgl)
library(Compositional)

# should rec to CT CODA using dx = dx.com.F
# the alpha-transformed data
# females, 0 values kept
# this wont work with the dataset containing zeroes, so we select the original data set here
# when  alpha = 0 and there is no multiplication by Helmert sub-matrix to reconcile back to the centred log ratio transformation
data <- dx.cent
data[!is.finite(data)] <- 0
alpha <- 0.9 #optimal alpha for females cardio is 0.9
helmert <- TRUE 
#while (alpha < 1){
  alpha.F <- alfa(data, alpha, helmert)$aff
  # this was checked and reconciles to clr.F as defined above
  rownames(alpha.F) <- as.character(years.fit.1)
  # Apply SVD to the alpha transformed data
  par <- svd(alpha.F, nu=3, nv=3)
  U <- par$u 
  V <- par$v
  S <- diag(par$d)
  
  bx<- V[,1]
  kt<- S[1,1]*U[,1] 
  bx2 <- V[,2]
  kt2 <- S[2,2]*U[,2]
  bx3 <- V[,3]
  kt3 <- S[3,3]*U[,3]
  
  R2.one <- cbind(par$d[1]^2/sum(par$d^2), par$d[2]^2/sum(par$d^2),par$d[3]^2/sum(par$d^2),par$d[4]^2/sum(par$d^2)) 
  
  # fit Arima model to kt
  order.one=c(1,1,1)
  fcst.one <- forecast(auto.arima(kt,max.d=1),h=ih)
  kt.fit.one <- fcst.one$fitted
  kt.for.one <- fcst.one$mean
  
  order.two=c(1,0,1)
  fcst.two <- forecast(auto.arima(kt2,max.d=1),h=ih)
  kt.fit.two <- fcst.two$fitted
  kt.for.two <- fcst.two$mean
  
  order.3=c(1,0,1)
  fcst.3 <- forecast(auto.arima(kt3,max.d=1),h=ih)
  kt.fit.3 <- fcst.3$fitted
  kt.for.3 <- fcst.3$mean
  
  gx.all <- matrix(gx,length(ages),k)
    
  #projections
  alpha.proj.fit.one <- matrix(kt,length(years.fit.1),1) %*% t(bx) + matrix(kt2,length(years.fit.1),1) %*% t(bx2) + matrix(kt3,length(years.fit.1),1) %*% t(bx3)
  alpha.proj.for.one <- matrix(c(kt, kt.for.one), length(years.fit.1)+ih,1) %*% t(bx) + matrix(c(kt2, kt.for.two), length(years.fit.1)+ih,1) %*% t(bx2) + matrix(c(kt3, kt.for.3), length(years.fit.1)+ih,1) %*% t(bx3)
  
  #The alpha transformation is a 1-1 transformation which maps the data inside the simplex to a subset of R(D-1) and vice versa, for values of alpha <> 0
  #The inverse of the alpha transformation will return NA for some values if the data lie outside the alpha-space.
  #The inverse transformation goes from A(D-1) to S(D-1) again. But once we fit the LC model, there is a probability that the data will be outside the alpha space

  #Invert the transformed data using modified alphainv (alfainv.mod)
  BK.alpha.proj.fit.one <- alfainv_mod(alpha.proj.fit.one, alpha, helmert)
  BK.alpha.proj.for.one <- alfainv_mod(alpha.proj.for.one, alpha, helmert)
  
  # ---
  alpha.proj.fit.one.gx.per <- perturbation(BK.alpha.proj.fit.one, gx, oper = "*")
  alpha.proj.for.one.gx.per <- perturbation(BK.alpha.proj.for.one, gx, oper = "*")
  
  alpha.proj.fit.one.gx <- alpha.proj.fit.one.gx.per / rowSums(alpha.proj.fit.one.gx.per)
  alpha.proj.for.one.gx <- alpha.proj.for.one.gx.per / rowSums(alpha.proj.for.one.gx.per)
  
  alpha.proj.fit.one.ar <- array((unclass(alpha.proj.fit.one.gx)), dim=c(length(years.fit.1),m,length(CoD)),dimnames=list(as.character(years.fit.1[1]:(years.fit.1[length(years.fit.1)])),as.character(ages),CoD))
  alpha.proj.for.one.ar <- array((unclass(alpha.proj.for.one.gx)), dim=c(length(years.fit.1)+ih,m,length(CoD)), dimnames=list(as.character(years.fit.1[1]:(years.for[length(years.for)])), as.character(ages),CoD))
  
  alpha.proj.for.one.tot <- alpha.proj.for.one.ar[,,1]
  for(i in 2:k){alpha.proj.for.one.tot <- alpha.proj.for.one.tot + alpha.proj.for.one.ar[,,i]}
  
  # prepare forecasts for plots
  pro.dea.alpha.LC <- t(apply(alpha.proj.for.one.ar, MARGIN=c(1,3),FUN=sum, na.rm = TRUE))
  
  colnames(pro.dea.alpha.LC) = total_years
  rownames(pro.dea.alpha.LC) = CoD
  plot_data_alpha <- pro.dea.alpha.LC[,1:43]
  
  colnames(plot_data_alpha) = years
  rownames(plot_data_alpha) = CoD.label
  
  datLong_alpha <- melt(data          = plot_data_alpha,
                        id.vars       = c(CoD.label),
                        measure.vars  = c(years),
                        variable.name = "year",
                        value.name    = "proportion")
  #datLong_alpha
  colnames(datLong_alpha) = c("CoD", "year", "proportion")
  
  # Alpha Data for plot 2
  plot2_data_alpha <- pro.dea.alpha.LC
  colnames(plot2_data_alpha) = total_years
  rownames(plot2_data_alpha) = CoD.label
  datLong2_alpha <- melt(data          = plot2_data_alpha,
                         id.vars       = c(CoD.label),
                         measure.vars  = c(total_years),
                         variable.name = "year",
                         value.name    = "proportion")
  colnames(datLong2_alpha) = c("CoD", "year", "proportion")
  
  # ========================================================
  # Results: Proportions of Deaths by cause
  # dimensions should be (number of years (16), number of age bands (9), number of causes (11))
  # Forecast using CLR, ILR, alpha-transformations
  # ========================================================
  rmse_clr_t <- Metrics::rmse(pro.dea.obs[,1:33],pro.dea.clr.LC[,1:33])
  mae_clr_t <- Metrics::mae(pro.dea.obs[,1:33],pro.dea.clr.LC[,1:33])
  
  rmse_ilr_t <- Metrics::rmse(pro.dea.obs[,1:33],pro.dea.ilr.LC[,1:33])
  mae_ilr_t <- Metrics::mae(pro.dea.obs[,1:33],pro.dea.ilr.LC[,1:33])
  
  rmse_alpha_t <- sqrt(mean((pro.dea.obs[,1:33] - pro.dea.alpha.LC[,1:33])^2 , na.rm = TRUE))
  mae_alpha_t <- sum(abs(pro.dea.obs[,1:33] - pro.dea.alpha.LC[,1:33]), na.rm = TRUE)/length(pro.dea.obs)
  
  #in sample RMSE
  print(rmse_clr_t)
  print(mae_clr_t)
  print(rmse_ilr_t)
  print(mae_ilr_t)
  print(rmse_alpha_t)
  print(mae_alpha_t)
  print(alpha)
  #alpha <- alpha + 0.1
  #}
  
  #error using all testing data
  rmse_alpha <- sqrt(mean((pro.dea.obs - pro.dea.alpha.LC[,1:43])^2 , na.rm = TRUE))
  rmse_clr <- Metrics::rmse(pro.dea.obs,pro.dea.clr.LC[,1:43])
  rmse_ilr <- Metrics::rmse(pro.dea.obs,pro.dea.ilr.LC[,1:43])
  mae_alpha <- sum(abs(pro.dea.obs - pro.dea.alpha.LC[,1:43]), na.rm = TRUE)/length(pro.dea.obs)
  mae_clr <- Metrics::mae(pro.dea.obs,pro.dea.clr.LC[,1:43])
  mae_ilr <- Metrics::mae(pro.dea.obs,pro.dea.ilr.LC[,1:43])
  
  #error from applying lee-carter (ignoring compositional data)
  rmse_lc <- Metrics::rmse(pro.dea.obs,pro.dea.lc[,1:43])
  mae_lc <- Metrics::mae(pro.dea.obs,pro.dea.lc[,1:43])
  print(rmse_lc*100)
  print(mae_lc*100)
  
  print(rmse_alpha*100)
  print(rmse_clr*100)
  print(rmse_ilr*100)
  print(mae_alpha*100)
  print(mae_clr*100)
  print(mae_ilr*100)

  # ========================================================
  # Plotting the forecast
  # ========================================================  

#reset the plot (if needed)
dev.off()
graphics.off()
#plot the forecast
col.cause <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC", "#28E2E5", "#000000")

g <- ggplot(data = datLong_obs, aes(x = CoD, y = proportion, group = year, colour = year))

theme_set(theme_bw())

observed_f <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  #scale_x_continuous(label = CoD.label)
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Observed") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none")


g <- ggplot(data = datLong_clr, aes(x = CoD, y = proportion, group = year, colour = year))

theme_set(theme_bw())

clr_f <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("CLR") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) + 
  theme(legend.position = "none")


g <- ggplot(data = datLong_ilr, aes(x = CoD, y = proportion, group = year, colour = year))

theme_set(theme_bw())

ilr_f <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("ILR") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none")

g <- ggplot(data = datLong_alpha, aes(x = CoD, y = proportion, group = year, colour = year))

theme_set(theme_bw())

alpha_f <- 
  g + geom_line() +
  labs(x = "Cause of Death", y = "Proportion", colour = "Year") +
  scale_y_continuous(trans='log10') +
  scale_color_viridis_c() +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  ggtitle("Alpha") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14))

#savepdf("us_comparison_f_0.9", width = 32, height = 12, toplines = 0.8)
par(mar=c(5, 4, 4, 8), xpd=TRUE)
observed_f + clr_f + ilr_f + alpha_f + 
  plot_layout(ncol = 4)
#dev.off()

# plot 2: time series by cause comparing alpha and observed

# Create plot
theme_set(theme_bw())

clr_proj_f <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong2_clr, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_obs, CoD %in% include.causes)) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2111, 15)) +
  #scale_x_continuous("Year", breaks = seq(1979, 2051, 4)) +
  #ggtitle("CLR (zeros replaced by 0.25 count)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none") #+
  #theme(axis.title.y=element_blank()) 

ilr_proj_f <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong2_ilr, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_obs, CoD %in% include.causes)) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2111, 15)) +
  #scale_x_continuous("Year", breaks = seq(1979, 2051, 4)) +
  #ggtitle("ILR (zeros replaced by 0.25 count)") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none") +
  theme(axis.title.y=element_blank()) 

alpha_proj_f<- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong2_alpha, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_obs, CoD %in% include.causes)) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2111, 15)) +
  #scale_x_continuous("Year", breaks = seq(1979, 2051, 4)) +
  #ggtitle("Alpha") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(axis.title.y=element_blank()) 
  #theme(legend.position = "none")

#lee-carter plot
lc_proj_f <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong_lc, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_obs, CoD %in% include.causes)) +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  scale_x_continuous("Year", breaks = seq(1979, 2111, 15)) +
  #scale_x_continuous("Year", breaks = seq(1979, 2051, 4)) +
  #ggtitle("Lee-Carter") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(legend.position = "none")

#savepdf("us_lc_alpha_proj_f", width = 32, height = 12, toplines = 0.8)
lc_proj_f + alpha_proj_f
plot_layout(ncol = 2)
#dev.off()


# Group plots together
#savepdf("us_all_proj_f_0.9_100yrs", width = 30, height = 12, toplines = 0.8)
par(mar=c(2, 2, 2, 2), xpd=TRUE)
clr_proj_f + ilr_proj_f + alpha_proj_f + 
  plot_layout(ncol = 3)
#dev.off()

# zero replacements: group plots together for CLR
#savepdf("us_clr_zeros_f", width = 33, height = 13, toplines = 0.8)
par(mar=c(2, 2, 2, 2), xpd=TRUE)
clr_proj_f + clr_proj_f_0.25 + clr_proj_f_0.5 + 
  plot_layout(ncol = 3)
#dev.off()

# zero replacements: group plots together for ILR
#savepdf("us_ilr_zeros_f", width = 33, height = 13, toplines = 0.8)
par(mar=c(2, 2, 2, 2), xpd=TRUE)
ilr_proj_f + ilr_proj_f_0.25 + ilr_proj_f_0.5 + 
  plot_layout(ncol = 3)
#dev.off()


# ========================================================
# interval forecast for the alpha-transformation
# using simple arima instead of LC-CODA
# ========================================================

#constructing interval forecast
interval_forecast_full <- alpha_fun_int(data, a = 0.9, ncomp_method = "fixed", ncomp_tuning = 0.001, fh = ih,
                                        fore_method = "rwf", B = 399, K = 1, sig = 20)
#optimised using cpd = 49.66667, alpha_val = 0.3

#interval_forecast_full <- apply(d_x_t_star_fore, c(1, 3), quantile, c((100 - sig)/200, 1 - (100 - sig)/200), na.rm = TRUE)

interval_forecast_lower_gx_per <- perturbation(interval_forecast_full[1,,], gx, oper = "*")
interval_forecast_upper_gx_per <- perturbation(interval_forecast_full[2,,], gx, oper = "*")

interval_forecast_lower_ar <- array((unclass(interval_forecast_lower_gx_per)), dim=c(ih,m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages),CoD))
interval_forecast_upper_ar <- array((unclass(interval_forecast_upper_gx_per)), dim=c(ih,m,length(CoD)),dimnames=list(as.character(years.for[1]:(years.for[length(years.for)])),as.character(ages),CoD))

# and interval forecasts
interval.lower.dea.alpha.LC <- t(apply(interval_forecast_lower_ar, MARGIN=c(1,3),FUN=sum, na.rm = TRUE))
interval.upper.dea.alpha.LC <- t(apply(interval_forecast_upper_ar, MARGIN=c(1,3),FUN=sum, na.rm = TRUE))

colnames(interval.lower.dea.alpha.LC) = years.for
rownames(interval.lower.dea.alpha.LC) = CoD
colnames(interval.upper.dea.alpha.LC) = years.for
rownames(interval.upper.dea.alpha.LC) = CoD

# plot interval forecasts
# lower bound
interval_lower_alpha <- interval.lower.dea.alpha.LC
colnames(interval_lower_alpha) = years.for
rownames(interval_lower_alpha) = CoD.label
datLong_interval_lower_alpha <- melt(data          = interval_lower_alpha,
                                     id.vars       = c(CoD.label),
                                     measure.vars  = c(years.for),
                                     variable.name = "year",
                                     value.name    = "proportion")
colnames(datLong_interval_lower_alpha) = c("CoD", "year", "proportion")

# upper bound
interval_upper_alpha <- interval.upper.dea.alpha.LC
colnames(interval_upper_alpha) = years.for
rownames(interval_upper_alpha) = CoD.label
datLong_interval_upper_alpha <- melt(data          = interval_upper_alpha,
                                     id.vars       = c(CoD.label),
                                     measure.vars  = c(years.for),
                                     variable.name = "year",
                                     value.name    = "proportion")
colnames(datLong_interval_upper_alpha) = c("CoD", "year", "proportion")

cpd(datLong_obs[397:516,], datLong_interval_lower_alpha[1:120,], datLong_interval_upper_alpha[1:120,], 90)
interval_score(datLong_obs[397:516,], datLong_interval_lower_alpha[1:120,], datLong_interval_upper_alpha[1:120,], 90)


# create plot
interval_alpha_f <- ggplot(data = NULL, aes(x = year, y = proportion, group = CoD, colour = CoD)) + 
  geom_line(data = subset(datLong_interval_lower_alpha, CoD %in% include.causes), linetype = "dashed") + 
  geom_line(data = subset(datLong_interval_upper_alpha, CoD %in% include.causes), linetype = "dashed") +
  geom_line(data = subset(datLong2_alpha, CoD %in% include.causes), linetype = "solid") +
  geom_line(data = subset(datLong_obs, CoD %in% include.causes), linetype = "dotted") +  
  labs(x = "Year", y = "Proportion", colour = "Cause of Death") +
  scale_y_continuous(trans='log10') +
  #scale_color_brewer(palette = "Paired") +
  scale_color_paletteer_d(palette = "rcartocolor::Vivid") +
  theme(axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 9)) + 
  #scale_x_continuous("Year", breaks = seq(1979, 2111, 4)) +
  scale_x_continuous("Year", breaks = seq(1979, 2051, 4)) +
  #ggtitle("Alpha") + 
  theme(plot.title = element_text(face = "bold",
                                  margin = margin(10, 0, 10, 0),
                                  size = 14)) +
  theme(axis.title.y=element_blank()) 


#savepdf("us_interval_f_0.9", width = 32, height = 10, toplines = 0.8)
par(mar=c(2, 2, 2, 2), xpd=TRUE)
interval_CLR_f + interval_ILR_f + interval_alpha_f +
plot_layout(ncol = 3)
#dev.off()
