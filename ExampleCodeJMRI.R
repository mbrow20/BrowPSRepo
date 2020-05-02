########################################################################################
##This is the meta-function to impute values for missing data using the Joint Modeling
##Random Intercepts (JMRI) model (van Buuren, 2012). Be sure to change the outcome
##variable to 'Y' and include as the second to last variable of the dataset. The last
##variable will be the student's id. We will use only one of the plausible values for this procedure.
####################################################################################
source(ImbensRubinFunction.R)
source("JointModelRandomInterceptsFULL.R")
source("FunctionsJMRI.R")
source("FunctionsCovSelection.R")
library(MCMCpack)
library(MatchIt)
library(Zelig)
library(caret)
library(mitml)
library(AppliedPredictiveModeling)
library(MASS)
library(lattice)
library(mice)
library(dplyr)
library(plyr)
library(lme4)
data1<-read.csv("USA2015_PSMODELREL.csv",header=TRUE,sep=",")
#Outcome is labeled as 'Y' is third to last, 'id' as second to last, and 'SCHID' as last.
#Outcome must be labeled 'Y', school id must be labeled 'SCHID', and dichotomous treatment
## variable as 'treat'.
dataFULL<-read.csv("USA2015_FULLMODELREL.csv",header=TRUE, sep=",")
JointModelRandomInterceptsPS(data1, varNames)
JointModelRandomInterceptsFULL(dataFULL)