##This sample code uses normalized final weights as recommeded by Method 2 of Pfefferman et al.(1998). 
##This technique allows for more accuracy variance estimates from data collected from multi-stage
##sampling. This is to be used specifically for CEM, which does not allow variable
##'treat' so it has been changed to 'treated'
rm(list=ls())
gc()
#setwd("C:/Users/mbrow20/Desktop/BrowDiss_Fall2019")
source("Functions8SR.R")
library(MCMCpack)
library(MatchIt)
library(Zelig)
library(caret)
library(dplyr)
library(plyr)
library(parallel)
library(MASS)
library(cem)
memory.size(5000000)
CovMatrixWt<-read.csv("CovMatrixWt1.csv",header=TRUE,sep=",")
data1<-read.csv("CAN2015_15Vars2LevWithNormWt.csv",header=TRUE,sep=",")
mu_beta<<-c(1.1,.02,-1,2,1.5,0.04,1.2,-1.3,0.8,-0.5,-1.4,0.3,-0.3,0.75,1.6,-0.98)#The 4th variable has a true beta coefficient of '2'###
SimulationWithCovMat(CovMatrixWt,data1,500,mu_beta)#This uses the PS 
