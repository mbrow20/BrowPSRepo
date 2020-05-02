##This R code will use BART for estimating average treatment effect (ATE). Code is
##based on example1.R code from Hill, though Hill estimates ATT.
SimulationWithCovMat=function(CovMatrixWt,data1,Sims,mu_beta,numTreat){
  options(warn=-1)
  library(MCMCpack)
  library(MatchIt)
  library(Zelig)
  library(caret)
  library(dplyr)
  library(plyr)
  library(BayesTree)
  source("Functions5SR.R")
  S2<<-as.matrix(CovMatrixWt,nrow=nrow(CovMatrixWt),ncol=ncol(CovMatrixWt))
  #S2<<-as.matrix(CovMatrixWt,nrow=5,ncol=5)
  set.seed(6)#ensures same results for random components (e.g., partitioning data, random variable selection, etc.)
  colnames(data1)[1]<<-"SC048Q01"
  #dpart<<-createDataPartition(data1$treat,p=0.2,list=F)
  #data2<<-data1[dpart,]
  data2<<-data1
  varNames<<-colnames(data2)
  n=nrow(data2)
  Q<<-length(unique(data2$STRATUM)) #Number of Explicit strata; Should be less than 'n'
  g2<<-as.data.frame(rep(0,n),ncol=1)
  g2[,1]<<-data2$STRATUM
  
  #These are actual Strata indicators used in the PISA data. Each stratum is composed is composed
  #of three explicit stratification variables: Province, Language, and School size (see OECD 2017, p. 72).
  g2$SCHID<<-data2$SCHID#These school ids correspond to actual schools. Make sure schools in dataset correspond to a variable labeled 'SCHID'.
  g2[,3]<<-data2$W_FSTUWT
  g2[,4]<<-data2$CNTSTUID
  
  g2<<-g2[order(g2[,1]),]
  colnames(g2)<<-c("g","SCHID","W_FSTUWT", "CNTSTUID")
  lis<<-as.data.frame(unique(g2$g))#These are unique STRATA in this sample
  lis2<<-as.vector(lis[,1])
  lis3<<-as.character(lis2)
  lis4<<-as.character(seq(1,length(lis3),1))#STRATA are re-labeled in sequential order. The 'mapvalues' function requires class 'character'
  g3<<-mapvalues(g2$g,c(lis3),c(lis4))#a case in the first vector that matches a case in the second vector will be replace by item in third vector
  g2<<-cbind(as.data.frame(as.numeric(g3)),g2$SCHID, g2$W_FSTUWT, g2$CNTSTUID)
  colnames(g2)<<-c("g","SCHID", "W_FSTUWT","CNTSTUID")
  listSCH1<<-as.data.frame(unique(g2$SCHID))
  listSCH2<<-as.vector(listSCH1[,1])
  listSCH3<<-as.character(listSCH2)
  listSCH4<<-as.character(seq(1,length(listSCH3),1))
  g4<<-mapvalues(g2$SCHID,c(listSCH3),c(listSCH4))
  g2<<-cbind(as.data.frame(as.numeric(g2[,1])),as.data.frame(as.numeric(g4)), g2[,3],g2[,4])
  colnames(g2)<<-c("g","SCHID","W_FSTUWT","CNTSTUID")
  g2 <<- g2[order(g2$CNTSTUID),]
  data2<<-data2[order(data2$CNTSTUID),]
  
  ##Create weights for implicit stratification (corresponding to the three levels of urbanicity) within the explicit strata
  for (i in 1:nrow(g2)){
    for (j in 1:Q){
      g2[i,j+4]<<-ifelse(g2$g[i]==j,sample((1:3)/4,1,replace=TRUE),0)}}
  ExplicitWT<<-rgamma(length(unique(g2$g)),1,1)
  g2$ExplWt<<-ExplicitWT[g2$g]
  SCHWT<<-rgamma(unique(g2$SCHID),1,1)
  g2$SCHWT<<-SCHWT[g2$SCHID]
  g2$FSCHWT<<-rowSums(g2[,c(3:ncol(g2))])
  g2$scaleweights<<-g2$FSCHWT*10#increase weights by a factor of ten to correspond to actual PISA weights
  g2$Female<<-data2$Female##As stated earlier, 'female' is one of the cells in the within school non-response; therefore, it is included as weight.
  g2$FSTUWT<<-g2$Female+g2$scaleweights#the weight of the 'female' variable is 1.
  #g2$CNTSTUID<<-data2$CNTSTUID
  #g2<<-g2[order(g2$SCHID),]
  colnames(data2)[1]="SC048Q01"
  #mu_beta<<-c(1.1,2,-1,1,1.5,0.04,1.2,-1.3,0.8,-0.5,-1.4)##The actual population beta coefficient is '2.0'. Recall, we are using wieghts, so inferences are made to the population.
  K<<-length(mu_beta)-2
  #make sigma^2 from inverse gamma distribution with parameters 1 and 1 BayesM package
  sigma_2<<-rgamma(nrow(data2),1,1)
  g2 <<- g2[order(g2$CNTSTUID),]#order the g2 dataframe by student ID. 
  #sigma_2<<-1
  sigma_2<<-sigma_2/g2$W_FSTUWT
  v<<-length(mu_beta)
  varNames<<-varNames[c(1:v-1)]
  Sigma_beta<<-riwish(v,S2)
  Beta<<-mvrnorm(1,mu_beta,Sigma_beta)
  U<<-mvrnorm(length(unique(g2$SCHID)),rep(0,K+2),diag(K+2))
  U<<-U[g2$SCHID,]
  length(unique(U[,1]))
  Ones<<-matrix(1,nrow(data2),1)
  data2<<-data2[order(data2$CNTSTUID),]#order data2 to match the order in the g2 dataframe
  Z<<-as.matrix(data2[,1:v-1])
  
  #Z<<-as.matrix(data2[,1:4])
  Z<<-cbind(Ones,Z)
  Y=Z %*% Beta + diag(Z %*% t(U)) + rnorm(n,0,sqrt(sigma_2))#One Simulated data of Y
  test_data<<-cbind(Z,Y,g2$FSTUWT,g2$SCHID,g2$W_FSTUWT,g2$CNTSTUID)
  varNames[1]="SC048Q01"
  colnames(test_data)=c("V1",varNames, "Y", "FSTUWT", "SCHID","W_FSTUWT","CNTSTUID")
  test_data<<-as.data.frame(test_data)
  write.table(test_data,file="test_data_NEW.csv",row.names=FALSE,col.names=TRUE, sep=",",append=FALSE)
  ##Examine the 'beta' of the treatment effect through a linear model
  data4<<-test_data[,c(1:length(varNames)+1)]
  data4<<-as.data.frame(data4)
  test_data<<-as.data.frame(test_data)
  indxer=which( colnames(data4)=="treat" )#--in case yourdataset has the variable 'treat' in a different column than one
  
  ###The following functions reorder 'data4' so the first variable is 'treat'##########################
  data4NotIndxer<<-data4[,-(indxer)]
  data4Indxer<<-data4[,indxer]
  data4<<-cbind(data4Indxer,data4NotIndxer)
  data4<<-as.data.frame(data4)
  colnames(data4)[1]="treat"
  View(data4)
  
  Sims=Sims
  R<-as.data.frame(rep(0,Sims),ncol=1)
  Beta_1<-matrix(rep(0,Sims*(K+2)),Sims,v)
  Beta_1<-as.data.frame(Beta_1)
  Beta_1=mvrnorm(Sims,mu_beta,Sigma_beta)
  View(Beta_1)
  int<-NULL
  mr<-NULL
  TFrame<-NULL
  indxer2=which(colnames(Z)=="treat" )
  Zt<-Z[,-c(1,indxer2)]
  
  for(i in 1:Sims){
    assign(paste0("Y",i), Z %*% Beta_1[i,] + diag(Z %*% t(U))+rnorm(nrow(Z),0,sqrt(sigma_2)))
    as.data.frame(assign(paste0("test_data",i),cbind(Z,get(paste0("Y",i)), test_data[,c(v+2,v+3, v+4)])))
    yp<-paste0("test_data",i)
    tmp2<-get(yp)
    colnames(tmp2)[v+1]<-"Y"
    assign(yp,tmp2)
    assign(paste0("test_data",i),as.data.frame(get(paste0("test_data",i)))[,-1])
    assign(colnames(get(paste0("test_data",i)))[ncol(Z)],"Y")
    assign(paste0("test_data",i), get(paste0("test_data",i))[,c(3,1:2,4:16)])
  }
  View(test_data1)
  
  for (i in 1:Sims){
    assign(paste0("data.bart",i), get(paste0("test_data",i))[get(paste0("test_data",i))$treat==0,])
    randSampCont<-sample(row.names(get(paste0("data.bart",i))),numTreat,replace=FALSE)
    assign(paste0("data.bart.Cont",i), get(paste0("data.bart",i))[randSampCont,])
    y=as.numeric(get(paste0("test_data",i))[,length(test_data1)])
    assign(paste0("bart.tot",i),bart(x.train=get(paste0("test_data",i)),   y.train=y,  x.test=get(paste0("data.bart.Cont",i))))
  }
  for (i in 1:Sims){
    assign(paste0("diffs",i),get(paste0("bart.tot",i))$yhat.train[,get(paste0("test_data",i))$treat==1]-get(paste0("bart.tot",i))$yhat.test)
    assign(paste0("mndiffs",i),apply(get(paste0("diffs",i)),1,mean))
    assign(paste0("bart.mean",i),mean(get(paste0("mndiffs",i))))
    assign(paste0("bart.sd",i), sd(get(paste0("mndiffs",i))))}
  bart.totONE<<-bart.tot1
  bart.totTWO<<-bart.tot2
  bart.meanOne<<-bart.mean1
  mndiffsOne<<-mndiffs1
  bart.meanTwo<<-bart.mean2
  mndiffsTwo<<-mndiffs2
  
  TFrame<-matrix(rep(0,Sims*2),ncol=2)
  TFrame<-as.data.frame(TFrame)
  for (i in 1:Sims){
    TFrame[i,1]<-get(paste0("bart.mean",i))
    TFrame[i,2]<-get(paste0("bart.sd",i))
    }
  names(TFrame)<-c("ATE","SD")
  View(TFrame)
  for(i in 1:Sims){
    TFrame$bias[i]=abs(2-TFrame$ATE[i])
  }
  TFrameBART<<-TFrame
  #View(TFrame)
  range.bias<-range(TFrame$bias)
  W<-sprintf("The mean bias from the BART model is %f, with a mean SD of  %f", mean(TFrame$bias),mean(TFrame$SD))
  W1<-sprintf("The range of the bias estimates is from %f to %f", range.bias[[1]], range.bias[[2]])
  print(W)
  print(W1)
  #bart.tot.One<<-bart.tot1
  actual.m<-as.vector(TFrame$ATE)
  expected.m<-as.vector(rep(2,Sims))
  RootMean.m<-RMSE(actual.m,expected.m)
  
  W2<-sprintf("The RMSE of the matched units is %f",RootMean.m)
  print(W2)
  TFrameOne<<-TFrame
  View(data.bart1)
  View(data.bart.Cont1)
  View(y)
}