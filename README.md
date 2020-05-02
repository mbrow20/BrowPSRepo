# BrowPSRepo
This is a repository for the propensity score analysis with 2015 Programme for International Student Assessment (PISA) data
# Functions and Data Files for Simulation
<p>These are the functions and the .csv data and datasets used for the simulation of the dissertation. <strong> Be sure to download all these files into a single directory (folder) on your computer.</strong> Below is a description of the different functions. Many of the functions are run through "ExampleCodeSim15Vars2LevNormWt.R". To utilize the various functions within this R script, specify using the source(), as in <strong>source("Functions13Sub.R")</strong>. In the last line of the 'ExampleCodeSim15Vars2Lev'-class of launch applications you can specify the number of iterations to run for the simulation. Here is an example,</p>
<p> <strong>SimulationWithCovMat(CovMatrixWt,data1,10,mu_beta)</strong>,</p> 

<p>where the number of simulations to run is ten. Here is an example where the launch' application "ExampleCodeSim15Vars2LevlNormWt.R' calls the function <strong><i>Functions6SR.R</i></strong> (a functions that performs data trimming) with 500 iterations:</p>

>rm(list=ls())<br>
>gc()<br>
>setwd("C:/Users/Mark V. Brow/Desktop/PISA2015/PropensityScoreAnalysis")#Change to your file directory<br>
>source(<strong>"Functions6SR.R"</strong>)<br>
>library(MCMCpack)<br>
>library(MatchIt)<br>
>library(Zelig)<br>
>library(caret)<br>
>library(dplyr)<br>
>library(plyr)<br>
>CovMatrixWt<-read.csv("CovMat15Vars2LevZWt.csv",header=TRUE,sep=",")<br>
>data1<-read.csv("CAN2015_15Vars2LevWithNormWt.csv",header=TRUE,sep=",")<br>
>mu_beta<<-c(1.1, .02, -1, 2, 1.5, 0.04, 1.2, -1.3, 0.8, -0.5, -1.4, 0.3, -0.3, 0.75, 1.6, -0.98)<br>
>#The 4th variable has a true beta coefficient of '2'###<br>
>SimulationWithCovMat(CovMatrixWt, data1, <strong>500</strong>, mu_beta)<br> 
<p>(Also, within the functions themselves, it is possible to run a truncated version of the dataset (for speed of computation; not recommended for analysis) that creates a random sample of the data proportional to treatment status. If <strong>.2</strong> is chosen, for example, it will randomly select 20% of the data, maintaining the ratio of treated and control observations. Within the functions themselves, you can comment out (# before a command line) code lines. Here is code for the full data set:</p>
<blockquote>
S2<<-as.matrix(CovMatrixWt,nrow=nrow(CovMatrixWt),ncol=ncol(CovMatrixWt))<br>
set.seed(6)#ensures same results for random components (e.g., partitioning data, random variable selection, etc.)<br>
colnames(data1)[1]<<-"SC048Q01"<br>
<strong>#</strong>dpart<<-createDataPartition(data1$treat, p=0.2, list=F)<br>
<strong>#</strong>data2<<-data1 [dpart,]<br>
data2<<-data1<br>
   </blockquote>
or for the partial data set (in this case, 20 percent):<br>
<blockquote>
dpart<<-createDataPartition(data1$treat, <strong>p=0.2</strong>, list=F)<br>
data2<<-data1 [dpart,]<br>
#data2<<-data1<br>
</blockquote>
<p>Search the function themselves to find these lines. With some functions, e.g.,<strong><i> Functions13Sub.R</strong></i>, <strong><i>Functions18.R</strong></i>, and<strong><i> Functions5VAR500.R</strong></i>, it is <strong>not possible to use the partition option as these functions require the full number of observations</strong>. Here is a brief explanation of the functions. A full explanation can be found in the <i>Appendix</i> of the dissertation. Most functions estimate RMSE, provide a table of balance statistics, display the iteration count, and output a t-stat(<i>df</i>) for sensitivity analysis (null being no confounding due to unobserved covariates).</p>
<ol>
<li><strong>Functions5SR.R</strong>--Performs matching on the data and is used with "ExampleCodeSim15Vars2LevNormWt.R".</li>
<li><strong>Functions5SR_2020CEM.R</strong>--Performs Coarsened Exact Matching (CEM) with various estimators (e.g., random forest, etc.). Used in conjunction with data set <strong>CAN2015_15Vars2LevNormWtCEM.csv</strong> and covariance matrix <strong> CovMatrixWtCEM.csv</strong>.matching on the data and is used with "ExampleCodeSim15Vars2LevNormWt.R".</li>
<li><strong>ExampleCodeSim15Vars2LevNormWtCEM.R</strong>--Estimates various Coarsened Exact Matching (CEM) models. Be sure to change the source to <strong>Functions5SR_2020CEM.R</strong>, <strong>data1<-</strong> to <strong>CAN2015_15Vars2LevWithNormWtCEM.csv</strong>, and the <strong>CovMatrix</strong> to <strong> CovMatrixWtCEM.csv</strong>.</li>
<li><strong>Functions6SR.R</strong>--Performs trimming on data based on algorithm by Imbens and Rubin (2015). Also used with  
   "ExampleCodeSim15Vars2LevNormWt.R".</li>
<li><strong>Functions7SR.R</strong>--Performs joint modeling imputation based on van Buuren (2012). The lauching application for this function is "ExampleMARMissingJOMO.R".</li>
<li><strong>Functions8SR.R</strong>--Performs matching on multiply imputed data sets (n = 5) through the SPSS multiple imputation function. Missing-at-random (MAR) data was created through a van Buuren (2012) algorithm.</li>
<li><strong>ExampleCodeSim15Vars2LevNormWtEMImpute</strong>--Used in conjunction with <strong>Functions5SR.R</strong> for performing mathcing on expectation maximization (EM) imputed data (SPSS).</li>
<li><strong>Functions12.R</strong>--Performs a bootstrapping procedure outlined by Zou et al. (2016) for determine the S.E. for the average treatment effect (ATE). This saves the treatment estimates and s.e. of the treatment effect in an Excel file in the current directory names "boot.est2.csv".</li>
<li><strong>Functions13Sub.R</strong>--Analyzes data subclassified by propensity score strata or ranges. Used in junction with its launch application, "ExampleCodeSime15Vars2LevSub.R". It is not allowed to used the partition feature on this data.</li>
<li><strong>Functions14.R</strong>--Performs mathcing with propensity score estimated through generalized boosted modeling (GBM; Ridgeway, 1999).Used with "ExampleCodeSim15Vars2LevNormWt.R".</li>
<li><strong>Functions16.R</strong>--Performs coarsened exact matching (CEM; Iacus et al., 2004) and is used with "ExampleCodeSim15Vars2LevNormWt.R".</li>
<li><strong>Functions18.R</strong>--Performs Bayesian Additive Regression Trees (BART; Hill, 2011) for estimating ATE. Used with "ExampleCodeSim15Vars2LevBART.R". The partition feature is not allowed with this function.</li>
<li><strong>FunctionsRobustRegression.R</strong>--Performs robust regression (i.e., with interactions and polynomials) with matched and unmatched data. Used with "ExampleCodeSim15Vars2LevNormWt.R".</li>
<li><strong>FunctionsDR.R</strong>--Performs "doubly"-robust estimation (Robins et al., 1995). Used with "ExampleCodeSim15Vars2LevNormWt.R".</li>
<li><strong>Functions5VAR500.R</strong>--Generates 500 matched data sets that are saved to the global environment and must be run prior to running <strong>FunctionsImbensRubinVarEst500.R</strong>. This function must be run with 500 iterations, i.e. <strong>SimulationWithCovMat(CovMatrixWt,data1,500,mu_beta)</strong>. </li>
<li><strong>FunctionsImbensRubinVarEst500.R</strong>--Estimates variance of ATE based on algorithm by Imbens and Rubin (2015). Must be run after running <strong>Functions5VAR500.R</strong>.</li>
<li><strong>SensitivityAnalysis</strong>--Performs sensitivity analysis on simulated data. Based on design-based model proposed by Imbens and Rubin (2015).</li>
<li><strong>SensitivityForApplicationData.R</strong>--Performs sensitivity analysis on the applied data sets, i.e., US and Canadian public and private school students, and US and Canadian religious private school and public school students.</li>
   <li><strong>ImbensRubin</strong>--Function to be used with <i>ExampleCodeJMRI.R</i>. After running <i>ExampleCodeJMRI.R</i>, "MatchIt" output files (named 'm.out1' through 'm.out5') are saved to the global environment. To run this function, include the matched output file along with the corresponding imputed (i.e., complete) dataset. For example, <strong>ImbensRubin</strong>(m.out3, USA2015_FULLIMPUTE3).</li>
   </ol>
Code for generating final imputed files along with Imbens and Rubin variance estimation is in <i>ImputationCovSelectionVarEst</i>. To generate the imputed finals and Imbens and Rubin variance estimation, use <i>ExampleCodeJMRI.R</i> and make sure include <i>ImbensRubinFunction.R</i>, <i>FunctionsCovSelection.R</i>, <i>JointModelRandomInterceptsFULL.R</i>, and <i>FunctionsJMRI.R</i> in the same directory. After running the code and generating the files in the directory, run <b>ImbensRubin(m.out2, CAN2015_FULLIMPUTE1)</b> (this is an example; modify code accordingly) to generate the Imbens and Rubin (2015) variance estimation for matched samples.</p>
