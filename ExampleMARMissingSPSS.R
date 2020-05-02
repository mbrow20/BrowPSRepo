##This function simulates with 5 imputed datasets from SPSS imputation. MAR missing were created on the variables
##SC048Q01 and ST034Q01##
rm(list=ls())
gc()
setwd("C:/Users/mbrow20/Desktop/BrowDiss_Fall2019")
source("Functions8SR.R")
library(MCMCpack)
library(MatchIt)
library(Zelig)
library(caret)
library(AppliedPredictiveModeling)
library(MASS)
library(lattice)
library(dplyr)
library(plyr)

CovMatrixWt<-read.csv("CovMat15Vars2LevZWt.csv",header=TRUE,sep=",")
data1<-read.csv("CAN2015_15Vars2LevWithNormWt.csv",header=TRUE,sep=",")
mu_beta<<-c(1.1,.02,-1,2,1.5,0.04,1.2,-1.3,0.8,-0.5,-1.4,0.3,-0.3,0.75,1.6,-0.98)#The 4th variable has a true beta coefficient of '2'###
SimulationWithCovMat(CovMatrixWt,data1,5,mu_beta)#This uses the PS model estimated without weights
