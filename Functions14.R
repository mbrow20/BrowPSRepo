##Functions14 uses the generalized boosted model (GBM) to estimate propensity scores. It then
##estimates the 'ATE' using the 'survey' package, with GBM-estimated propensity scores
##as weights in the model. You only need one 'ps.binary' model to performbalance diagnostics,
##since the model for balance does not require the outcome variable. See 'twang' User Guide
##for all the graphs and tables offered.
##This function is to be used with the Normalized final student weights recommended by
##Pfefferman et al. (1998) for multi-stage sampling designs
##Functions5.R is for the 15 Variable 2 Level Simulation##

PS.CovSelection.Function<-function(data4){
  dev=NULL
  max_loc=NULL
  min_loc=NULL
  dev_sub=NULL
  factors=NULL
  formula=NULL
  dev_matrix=matrix(0,ncol(data4),ncol(data4))
  dev_matrix1=NULL
  storeMaxLoc=NULL
  indexMax=NULL
  iter=1:ncol(data4)
  devPrior=NULL
  dev=NULL
  max_loc=NULL
  min_loc=NULL
  dev_sub=NULL
  factors=NULL
  formula=NULL
  dev_matrix=matrix(0,ncol(data4),ncol(data4))
  dev_matrix1=NULL
  storeMaxLoc=NULL
  indexMax=NULL
  iter=1:ncol(data4)
  devPrior=NULL
  listnames=NULL
  data4<-na.omit(data4)
  
  ##This will give you the deviance statistic--subtract the deviance from the nested (smaller) model with the larger model for the likelihood ratio statistic--the largest value has greatest effect
  for (j in 2:ncol(data4)) {
    assign(paste0("dev", j), dev)
    assign(paste0("max_loc", j), max_loc)
    assign(paste0("dev_sub", j), dev_sub)
    assign(paste0("min_loc", j), min_loc)
    assign(paste0("factors", j), factors)
    assign(paste0("formula", j), formula)
    assign(paste0("storeMaxLoc", j), storeMaxLoc)
    assign(paste0("indexMax",j), indexMax)
  }
  for( i in 2:ncol(data4)) {
    factors[i]<-colnames(data4)[i]
    assign(paste0("formula",1), as.formula(paste("treat~", paste(factors[i], collapse="+"))))
    assign(paste0("Variable", i), glm(formula1,data=data4,family=binomial(link="logit"))$deviance)
    dev[i]=get(paste0("Variable",i))
  }
  
  for (i in 2:ncol(data4)) {
    assign(paste0("dev",1),dev)
    assign(paste0("min_loc",1),which.min(get(paste0("dev",1))))
    dev_matrix[,1]=get(paste0("dev",1))
    storeMaxLoc[1]=factors[get(paste0("min_loc",1))]
    Variable1=dev[min_loc1]
  }
  
  for (j in 2:2){
    for (i in 2:ncol(data4)){
      assign(paste0("Variable", i), glm(data4[,1]~data4[,get(paste0("min_loc",1))]+data4[,i],data=data4,family=binomial(link="logit"))$deviance)
      dev[i]=get(paste0("Variable",i))
      assign(paste0("dev",2),dev)
      dev_sub2=Variable1-dev2
      dev_matrix[,2]=dev_sub2
      max_loc2=which.max(dev_sub2)
      assign(paste0("indexMax", j),get(paste0("dev_sub",j))[get(paste0("max_loc",j))])
      storeMaxLoc[2]=factors[get(paste0("max_loc",2))]
      assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
    }
  }
  #####################################################################################################
  for(t in 3:ncol(data4)) {
    for (j in t:t) {
      for (i in 2:ncol(data4)){
        assign(paste0("storeMaxLoc",j),storeMaxLoc)
        y=paste(get(paste0("storeMaxLoc",j)),collapse="+")
        assign(paste0("formula",j), as.formula(paste("treat~", paste(y,paste("+"),paste(factors[i], collapse="+")))))
        assign(paste0("Variable", i), glm(get(paste0("formula", j)), data=data4,family=binomial(link="logit"))$deviance)
        dev[i]=get(paste0("Variable",i))
        assign(paste0("dev",j),dev)
        num=j-1
        devPrior=get(paste0("indexMax",num))
        assign(paste0("dev_sub", j),devPrior-get(paste0("dev",j)))
        dev_matrix[,j]=get(paste0("dev_sub",j))
        assign(paste0("max_loc", j), which.max(get(paste0("dev_sub",j))))
        assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
      }
      storeMaxLoc[j]=factors[get(paste0("max_loc",j))]
    }
  }
  
  ##The difference in the deviance statistic from the previous best nested model should be greater than '1'###
  MatrixLen<-length(unique(storeMaxLoc))
  a<-matrix(0,MatrixLen,MatrixLen)
  a<-ifelse(dev_matrix<=1,1,0)
  a<-as.data.frame(a)
  a<-a[-1,-ncol(a)]
  c<-colSums(a)
  d=length(which(c<ncol(a)))
  StoreMaxLocGreaterThanOne=unique(storeMaxLoc)[c(1:d)]
  
  ################################################Check the interaction terms. 'CL' level is 2.71 for interaction terms (see p. 288 in Imbens and Rubin (2015))#######
  ###################################################################################################################################################################
  #b<-unique(storeMaxLoc)
  b<-unique(StoreMaxLocGreaterThanOne)
  ###Need to create a data set with no categorical variables. Quadratics of categorical variables yield the same results the original variables############
  ###The best option is to create a separate data set named 'DataNonCat' in which you remove the categorical variables. This code will remove##############
  ###'integer' class variables, though if the 'integer' class is not simply binary but ordinal, it will unnecessarily remove these variables also##########
  
  #data4$age<-as.numeric(data4$age)
  #data4$educ<-as.numeric(data4$educ)
  
  DataNonCat<-data4[,b]
  DataNonCat1<-NULL
  for (i in 1:ncol(DataNonCat)){
    DataNonCat1[i]=is.integer(DataNonCat[,i])}
  DataNonCat<-DataNonCat[,c(!DataNonCat1)]
  c<-paste(b,collapse="+")
  z=paste(c,")")
  z=paste("(",z)
  z=paste(z,"^2")
  FormulaNew=as.formula(paste("treat", z, sep=" ~ "))
  VarNamesInteraction=variable.names(glm(FormulaNew,data=data4,family="binomial"(link="logit")))
  listnames<-as.list(VarNamesInteraction)[-1]
  ListNamesCat=names(DataNonCat)
  ListNamesNew=listnames[c(1:length(b))]
  ListInterAct=listnames[c(length(b)+1:(length(listnames)-length(b)))]
  
  Orig.Var<-paste(ListNamesNew,collapse="+")
  Orig.Var<-paste(Orig.Var,"+")
  Quad.Var<-paste(ListNamesCat,"^2)")
  Quad.Var<-paste("I(",Quad.Var)
  List.Quad.InterAct<-c(ListInterAct,Quad.Var)
  formula999<-paste("treat~",Orig.Var)
  formula3<-as.formula(paste(formula999,paste(List.Quad.InterAct,collapse="+")))
  NumDim<-length(List.Quad.InterAct)
  factors2<-NULL
  storeMaxLoc2<-NULL
  dev_matrix2=matrix(0,NumDim,NumDim)
  
  for (j in 1:NumDim) {
    assign(paste0("dev", j), dev)
    assign(paste0("max_loc", j), max_loc)
    assign(paste0("dev_sub", j), dev_sub)
    assign(paste0("min_loc", j), min_loc)
    assign(paste0("factors2", j), factors)
    assign(paste0("listnames",j), listnames)
    assign(paste0("formula", j), formula)
    assign(paste0("storeMaxLoc2", j), storeMaxLoc)
    assign(paste0("indexMax",j), indexMax)
  }
  for( i in 1:NumDim) {
    factors2[i]<-List.Quad.InterAct[i]
    assign(paste0("formula",i), as.formula(paste(formula999,paste(List.Quad.InterAct[i],collapse="+"))))
    assign(paste0("Variable", i), glm(get(paste0("formula",i)),data=data4,family=binomial(link="logit"))$deviance)
    dev[i]=get(paste0("Variable",i))
  }
  for (i in 1:NumDim) {
    assign(paste0("dev",1),dev)
    assign(paste0("min_loc",1),which.min(get(paste0("dev",1))))
    dev_matrix2[,1]=get(paste0("dev",1))
    storeMaxLoc2[1]=factors2[get(paste0("min_loc",1))]
    Variable1=dev[min_loc1]
  }
  #######################################################################################################################################################
  for (j in 2:2){
    for (i in 2:NumDim){
      newlist1=List.Quad.InterAct[min_loc1]
      newlist2=unique(c(newlist1,List.Quad.InterAct))
      factors2[i]=newlist2[i]
      ###need to include the variable from the 1st step and other equation like above##############################
      assign(paste0("formula",i), as.formula(paste(formula999,paste(newlist2[c(1,i)],collapse="+"))))
      assign(paste0("Variable", i+1), glm(get(paste0("formula",i)),data=data4,family=binomial(link="logit"))$deviance)
      dev[i]=get(paste0("Variable",i+1))
      assign(paste0("dev",2),dev)
      dev_sub2=Variable1-dev2
      dev_matrix2[,2]=dev_sub2
      max_loc2=which.max(dev_sub2)
      assign(paste0("indexMax", j),get(paste0("dev_sub",j))[get(paste0("max_loc",j))])
      storeMaxLoc2[2]=factors2[get(paste0("max_loc",2))]
      assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
    }
  }
  
  ####################################################################################################################################################
  options(warn=-1)
  ##The warnings were suppressed since you may encounter fitted values from the logistic regression of either 0 or 1. The warnings is 
  ##turned on again after the logistic regression with interactions.#################################################################
  
  formula1000<-paste0(formula999,storeMaxLoc2[1])
  formula1000<-paste0(formula1000, paste("+"), paste(storeMaxLoc2[2]),paste("+"))
  factors2<-unique(c(storeMaxLoc2,factors2))
  for(t in 3:NumDim) {
    for (j in t:t) {
      for (i in 3:NumDim){
        assign(paste0("storeMaxLoc2",j),storeMaxLoc2)
        y=paste(get(paste0("storeMaxLoc2",j)),collapse="+")
        assign(paste0("formula",i), as.formula(paste(formula999,paste(y,paste("+")),paste(factors2[(i)], collapse="+"))))
        #assign(paste0("formula",i), as.formula(paste(formula1000,paste(factors2[i], collapse="+"))))
        assign(paste0("Variable", i), glm(get(paste0("formula", i)), data=data4,family=binomial(link="logit"))$deviance)
        dev[i]=get(paste0("Variable",i))
        assign(paste0("dev",j),dev)
        num=j-1
        assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
        devPrior=get(paste0("indexMax",num))
        assign(paste0("dev_sub", j),devPrior-get(paste0("dev",j)))
        dev_matrix2[,j]=get(paste0("dev_sub",j))
        assign(paste0("max_loc", j), which.max(get(paste0("dev_sub",j))))
        #assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
      }
      storeMaxLoc2[j]=factors2[get(paste0("max_loc",j))]
    }
  }
  #options(warn=0)
  ##The difference in the deviance statistic from the previous best nested model should be greater than '2.71' for Interaction terms###
  
  MatrixLen2<-length(unique(storeMaxLoc2))
  a.b<-matrix(0,MatrixLen2,MatrixLen2)
  a.b<-ifelse(dev_matrix2<=2.71,1,0)
  a.b<-as.data.frame(a.b)
  a.b<-a.b[-1,-ncol(a.b)]
  c.b<-colSums(a.b)
  d.b=length(which(c.b<ncol(a.b)))
  StoreMaxLocGreaterThan2.71=unique(storeMaxLoc2)[c(1:d.b)]
  
  interActionVars<<-as.data.frame(length(unique(StoreMaxLocGreaterThan2.71)),1)#Saved as global variables for use outside function
  linearVars<<-as.data.frame(length(unique(StoreMaxLocGreaterThanOne)),1)#Saved as global variables for use outside function
  
  interActionVars<-as.data.frame(length(unique(StoreMaxLocGreaterThan2.71)),1)#Saved as local variables for inside the function
  linearVars<-as.data.frame(length(unique(StoreMaxLocGreaterThanOne)), 1)#Saved as local variables for inside the function
  
  for(g in 1:length(unique(StoreMaxLocGreaterThanOne))){
    linearVars[g,1]<-unique(StoreMaxLocGreaterThanOne)[g]
    print(unique(StoreMaxLocGreaterThanOne)[g])}
  for(g in 1:length(unique(StoreMaxLocGreaterThan2.71))){
    interActionVars[g,1]<-unique(StoreMaxLocGreaterThan2.71)[g]
    print(unique(StoreMaxLocGreaterThan2.71)[g])}
  
  for(g in 1:length(unique(StoreMaxLocGreaterThanOne))){
    linearVars[g,1]<<-unique(StoreMaxLocGreaterThanOne)[g]}
  #print(unique(StoreMaxLocGreaterThanOne)[g])}
  for(g in 1:length(unique(StoreMaxLocGreaterThan2.71))){
    interActionVars[g,1]<<-unique(StoreMaxLocGreaterThan2.71)[g]}
  #print(unique(StoreMaxLocGreaterThan2.71)[g])}
  
  linearVars<<-linearVars[,1]#Save as global variables
  interactionVars<<-interActionVars[,1]#Save as global variables
  
  linearVars<<-linearVars[,1]#Save as global variables
  interactionVars<<-interActionVars[,1]#Save as global variables
  
  PSformulaUnWt<<-as.formula(paste("treat~",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+"))))
  formulaRegUnWt<<-as.formula(paste("Y~treat+",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+"))))
}

#######################################################################################################
#######################################################################################################

PS.CovSelection.Function.Weighted<-function(data, weights){
  
  dev=NULL
  max_loc=NULL
  min_loc=NULL
  dev_sub=NULL
  factors=NULL
  formula=NULL
  dev_matrix=matrix(0,ncol(data4),ncol(data4))
  dev_matrix1=NULL
  storeMaxLoc=NULL
  indexMax=NULL
  iter=1:ncol(data4)
  devPrior=NULL
  dev=NULL
  max_loc=NULL
  min_loc=NULL
  dev_sub=NULL
  factors=NULL
  formula=NULL
  dev_matrix=matrix(0,ncol(data4),ncol(data4))
  dev_matrix1=NULL
  storeMaxLoc=NULL
  indexMax=NULL
  iter=1:ncol(data4)
  devPrior=NULL
  listnames=NULL
  data4<-na.omit(data4)
  
  ##This will give you the deviance statistic--subtract the deviance from the nested (smaller) model with the larger model for the likelihood ratio statistic--the largest value has greatest effect
  for (j in 2:ncol(data4)) {
    assign(paste0("dev", j), dev)
    assign(paste0("max_loc", j), max_loc)
    assign(paste0("dev_sub", j), dev_sub)
    assign(paste0("min_loc", j), min_loc)
    assign(paste0("factors", j), factors)
    assign(paste0("formula", j), formula)
    assign(paste0("storeMaxLoc", j), storeMaxLoc)
    assign(paste0("indexMax",j), indexMax)
  }
  for( i in 2:ncol(data4)) {
    factors[i]<-colnames(data4)[i]
    assign(paste0("formula",1), as.formula(paste("treat~", paste(factors[i], collapse="+"))))
    assign(paste0("Variable", i), glm(formula1,data=data4,weights=weights,family="quasibinomial")$deviance)
    dev[i]=get(paste0("Variable",i))
  }
  
  for (i in 2:ncol(data4)) {
    assign(paste0("dev",1),dev)
    assign(paste0("min_loc",1),which.min(get(paste0("dev",1))))
    dev_matrix[,1]=get(paste0("dev",1))
    storeMaxLoc[1]=factors[get(paste0("min_loc",1))]
    Variable1=dev[min_loc1]
  }
  
  for (j in 2:2){
    for (i in 2:ncol(data4)){
      assign(paste0("Variable", i), glm(data4[,1]~data4[,get(paste0("min_loc",1))]+data4[,i],data=data4,weights=weights,family="quasibinomial")$deviance)
      dev[i]=get(paste0("Variable",i))
      assign(paste0("dev",2),dev)
      dev_sub2=Variable1-dev2
      dev_matrix[,2]=dev_sub2
      max_loc2=which.max(dev_sub2)
      assign(paste0("indexMax", j),get(paste0("dev_sub",j))[get(paste0("max_loc",j))])
      storeMaxLoc[2]=factors[get(paste0("max_loc",2))]
      assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
    }
  }
  #####################################################################################################
  for(t in 3:ncol(data4)) {
    for (j in t:t) {
      for (i in 2:ncol(data4)){
        assign(paste0("storeMaxLoc",j),storeMaxLoc)
        y=paste(get(paste0("storeMaxLoc",j)),collapse="+")
        assign(paste0("formula",j), as.formula(paste("treat~", paste(y,paste("+"),paste(factors[i], collapse="+")))))
        assign(paste0("Variable", i), glm(get(paste0("formula", j)), data=data4,weights=weights,family="quasibinomial")$deviance)
        dev[i]=get(paste0("Variable",i))
        assign(paste0("dev",j),dev)
        num=j-1
        devPrior=get(paste0("indexMax",num))
        assign(paste0("dev_sub", j),devPrior-get(paste0("dev",j)))
        dev_matrix[,j]=get(paste0("dev_sub",j))
        assign(paste0("max_loc", j), which.max(get(paste0("dev_sub",j))))
        assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
      }
      storeMaxLoc[j]=factors[get(paste0("max_loc",j))]
    }
  }
  
  ##The difference in the deviance statistic from the previous best nested model should be greater than '1'###
  MatrixLen<-length(unique(storeMaxLoc))
  a<-matrix(0,MatrixLen,MatrixLen)
  a<-ifelse(dev_matrix<=1,1,0)
  a<-as.data.frame(a)
  a<-a[-1,-ncol(a)]
  c<-colSums(a)
  d=length(which(c<ncol(a)))
  StoreMaxLocGreaterThanOne=unique(storeMaxLoc)[c(1:d)]
  
  ################################################Check the interaction terms. 'CL' level is 2.71 for interaction terms (see p. 288 in Imbens and Rubin (2015))#######
  ###################################################################################################################################################################
  b<-unique(storeMaxLoc)
  
  ###Need to create a data set with no categorical variables. Quadratics of categorical variables yield the same results the original variables############
  ###The best option is to create a separate data set named 'DataNonCat' in which you remove the categorical variables. This code will remove##############
  ###'integer' class variables, though if the 'integer' class is not simply binary but ordinal, it will unnecessarily remove these variables also##########
  
  #data4$age<-as.numeric(data4$age)
  #data4$educ<-as.numeric(data4$educ)
  
  DataNonCat<-data4[,b]
  DataNonCat1<-NULL
  for (i in 1:ncol(DataNonCat)){
    DataNonCat1[i]=is.integer(DataNonCat[,i])}
  DataNonCat<-DataNonCat[,c(!DataNonCat1)]
  c<-paste(b,collapse="+")
  z=paste(c,")")
  z=paste("(",z)
  z=paste(z,"^2")
  FormulaNew=as.formula(paste("treat", z, sep=" ~ "))
  VarNamesInteraction=variable.names(glm(FormulaNew,data=data4,weights=weights,family="quasibinomial"))
  listnames<-as.list(VarNamesInteraction)[-1]
  ListNamesCat=names(DataNonCat)
  ListNamesNew=listnames[c(1:length(b))]
  ListInterAct=listnames[c(length(b)+1:(length(listnames)-length(b)))]
  
  Orig.Var<-paste(ListNamesNew,collapse="+")
  Orig.Var<-paste(Orig.Var,"+")
  Quad.Var<-paste(ListNamesCat,"^2)")
  Quad.Var<-paste("I(",Quad.Var)
  List.Quad.InterAct<-c(ListInterAct,Quad.Var)
  formula999<-paste("treat~",Orig.Var)
  formula3<-as.formula(paste(formula999,paste(List.Quad.InterAct,collapse="+")))
  NumDim<-length(List.Quad.InterAct)
  factors2<-NULL
  storeMaxLoc2<-NULL
  dev_matrix2=matrix(0,NumDim,NumDim)
  
  for (j in 1:NumDim) {
    assign(paste0("dev", j), dev)
    assign(paste0("max_loc", j), max_loc)
    assign(paste0("dev_sub", j), dev_sub)
    assign(paste0("min_loc", j), min_loc)
    assign(paste0("factors2", j), factors)
    assign(paste0("listnames",j), listnames)
    assign(paste0("formula", j), formula)
    assign(paste0("storeMaxLoc2", j), storeMaxLoc)
    assign(paste0("indexMax",j), indexMax)
  }
  for( i in 1:NumDim) {
    factors2[i]<-List.Quad.InterAct[i]
    assign(paste0("formula",i), as.formula(paste(formula999,paste(List.Quad.InterAct[i],collapse="+"))))
    assign(paste0("Variable", i), glm(get(paste0("formula",i)),data=data4,weights=weights, family="quasibinomial")$deviance)
    dev[i]=get(paste0("Variable",i))
  }
  for (i in 1:NumDim) {
    assign(paste0("dev",1),dev)
    assign(paste0("min_loc",1),which.min(get(paste0("dev",1))))
    dev_matrix2[,1]=get(paste0("dev",1))
    storeMaxLoc2[1]=factors2[get(paste0("min_loc",1))]
    Variable1=dev[min_loc1]
  }
  #######################################################################################################################################################
  for (j in 2:2){
    for (i in 2:NumDim){
      newlist1=List.Quad.InterAct[min_loc1]
      newlist2=unique(c(newlist1,List.Quad.InterAct))
      factors2[i]=newlist2[i]
      ###need to include the variable from the 1st step and other equation like above##############################
      assign(paste0("formula",i), as.formula(paste(formula999,paste(newlist2[c(1,i)],collapse="+"))))
      assign(paste0("Variable", i+1), glm(get(paste0("formula",i)),data=data4,weights=weights, family="quasibinomial")$deviance)
      dev[i]=get(paste0("Variable",i+1))
      assign(paste0("dev",2),dev)
      dev_sub2=Variable1-dev2
      dev_matrix2[,2]=dev_sub2
      max_loc2=which.max(dev_sub2)
      assign(paste0("indexMax", j),get(paste0("dev_sub",j))[get(paste0("max_loc",j))])
      storeMaxLoc2[2]=factors2[get(paste0("max_loc",2))]
      assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
    }
  }
  
  ####################################################################################################################################################
  options(warn=-1)
  ##The warnings were suppressed since you may encounter fitted values from the logistic regression of either 0 or 1. The warnings is 
  ##turned on again after the logistic regression with interactions.#################################################################
  
  formula1000<-paste0(formula999,storeMaxLoc2[1])
  formula1000<-paste0(formula1000, paste("+"), paste(storeMaxLoc2[2]),paste("+"))
  factors2<-unique(c(storeMaxLoc2,factors2))
  for(t in 3:NumDim) {
    for (j in t:t) {
      for (i in 3:NumDim){
        assign(paste0("storeMaxLoc2",j),storeMaxLoc2)
        y=paste(get(paste0("storeMaxLoc2",j)),collapse="+")
        assign(paste0("formula",i), as.formula(paste(formula999,paste(y,paste("+")),paste(factors2[(i)], collapse="+"))))
        assign(paste0("Variable", i), glm(get(paste0("formula", i)), data=data4,weights=weights,family="quasibinomial")$deviance)
        dev[i]=get(paste0("Variable",i))
        assign(paste0("dev",j),dev)
        num=j-1
        assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
        devPrior=get(paste0("indexMax",num))
        assign(paste0("dev_sub", j),devPrior-get(paste0("dev",j)))
        dev_matrix2[,j]=get(paste0("dev_sub",j))
        assign(paste0("max_loc", j), which.max(get(paste0("dev_sub",j))))
        #assign(paste0("indexMax",j),get(paste0("dev",j))[which.min(get(paste0("dev",j)))])
      }
      storeMaxLoc2[j]=factors2[get(paste0("max_loc",j))]
    }
  }
  options(warn=0)
  ##The difference in the deviance statistic from the previous best nested model should be greater than '2.71' for Interaction terms###
  
  
  MatrixLen2<-length(unique(storeMaxLoc2))
  a.b<-matrix(0,MatrixLen2,MatrixLen2)
  a.b<-ifelse(dev_matrix2<=2.71,1,0)
  a.b<-as.data.frame(a.b)
  a.b<-a.b[-1,-ncol(a.b)]
  c.b<-colSums(a.b)
  d.b=length(which(c.b<ncol(a.b)))
  StoreMaxLocGreaterThan2.71=unique(storeMaxLoc2)[c(1:d.b)]
  interActionVarsWt<<-as.data.frame(length(unique(StoreMaxLocGreaterThan2.71)),1)#Saved as global variables for use outside function
  linearVarsWt<<-as.data.frame(length(unique(StoreMaxLocGreaterThanOne)),1)#Saved as global variables for use outside function
  
  for(g in 1:length(unique(StoreMaxLocGreaterThanOne))){
    linearVarsWt[g,1]<<-unique(StoreMaxLocGreaterThanOne)[g]
    print(unique(StoreMaxLocGreaterThanOne)[g])}
  for(g in 1:length(unique(StoreMaxLocGreaterThan2.71))){
    interActionVarsWt[g,1]<<-unique(StoreMaxLocGreaterThan2.71)[g]
    print(unique(StoreMaxLocGreaterThan2.71)[g])}
  linearVarsWt<<-linearVarsWt[,1]
  interActionVarsWt<<-interActionVarsWt[,1]
  PSformulaWt<<-as.formula(paste("treat~",paste(linearVarsWt,collapse="+"),paste("+",paste(interActionVarsWt,collapse="+"))))
  formulaRegWt<<-as.formula(paste("Y~treat+",paste(linearVarsWt,collapse="+"),paste("+",paste(interActionVarsWt,collapse="+"))))
}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
##This function simulates data from actual data. Need to provide a variance/covariance matrix (including the intercept), preferrably from
##a nested or hierarchical model...HLM will output a variance/covariance matrix based on a hierarchical model. This is intended for
##large datasets, as the function will create a subset of 20% of the data proportional on the treatment variable. The treatment variable
##should be the first variable in the dataset. The function will automatically label it 'treat'. Other variables you should include 
#are 'STRATUM', 'SCHID', and 'Female' (binary coded 1 for female). The variable 'Sims' is the number of simulations you would like to perform.

SimulationWithCovMat=function(CovMatrixWt,data1,Sims,mu_beta){
  library(MCMCpack)
  library(MatchIt)
  library(Zelig)
  library(caret)
  library(dplyr)
  library(plyr)
  library(twang)
  library(survey)
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
  
  ########We use the PS.CovSelection.Function to choose covariates for the PS model###############
  PS.CovSelection.Function(data4)
  formulaRegUnWt<<-as.formula(paste("Y~",paste(varNames,collapse="+")))
  formulaRegUnWt<-as.formula(paste("Y~",paste(varNames,collapse="+")))
  #View(formulaRegUnWt)
  View(test_data)
  test_data<<-read.csv("test_data_NEW.csv",header=TRUE,sep=",")
  test_data<<-as.data.frame(test_data)
  test_data$W_FSTUWT<<-data2$W_FSTUWT
  test_data<-read.csv("test_data_NEW.csv",header=TRUE,sep=",")
  test_data<-as.data.frame(test_data)
  test_data$W_FSTUWT<-data2$W_FSTUWT
  is.data.frame(test_data)
  LinModel<-lm(formulaRegUnWt,data=test_data, weights=W_FSTUWT)
  LinModel<<-lm(formulaRegUnWt,data=test_data, weights=W_FSTUWT)
  print(summary(LinModel))
  print(confint(LinModel))
  
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
  
  idVar<-seq(1,nrow(test_data))
  
  for(i in 1:Sims){
    assign(paste0("Y",i), Z %*% Beta_1[i,] + diag(Z %*% t(U))+rnorm(nrow(Z),0,sqrt(sigma_2)))
    as.data.frame(assign(paste0("test_data",i),cbind(Z,get(paste0("Y",i)), test_data[,c(v+2,v+3, v+4)])))
    yp<-paste0("test_data",i)
    tmp2<-get(yp)
    colnames(tmp2)[v+1]<-"Y"
    assign(yp,tmp2)
    #as.data.frame(assign(paste0("test_data",i),cbind(Z,get(paste0("Y",i)), data2[,c(6,9)])))
    assign(paste0("test_data",i),as.data.frame(get(paste0("test_data",i)))[,-1])
    assign(paste0("predVars",i),as.data.frame(Zt))
    assign(colnames(get(paste0("predVars",i)))[1],"Treat")
    assign(colnames(get(paste0("test_data",i)))[ncol(Z)],"Y")
    nm<-paste0("predVars",i)
    #nm<-paste0("test_data",i)
    tmp<-get(nm)
    #colnames(tmp) <- c("Treat","V1","V2","V3","Y","FSTUWT","SCHID")
    predNames<-colnames(tmp)
    #predNames<-colnames(tmp)[c(2,3,4)]
    predNames2<-c("treat",predNames)
    assign(nm,tmp)
    assign(paste0("formula",1), as.formula(paste("treat~", paste(predNames, collapse="+"))))
    assign(paste0("formula",2), as.formula(paste("Y~", paste(predNames2, collapse="+"))))
    assign(paste0("formula",3), as.formula(paste("Y~", paste(predNames, collapse="+"))))
    assign(paste0("formula",4), as.formula(paste("treat~", paste(predNames, collapse="+"))))
    assign(paste0("PSformulaUnWt",1), as.formula(paste("treat~",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+")))))
    assign(paste0("mr_orig",i), lm(formula2,data=get(paste0("test_data",i)),weights = get(paste0("test_data",i))$W_FSTUWT))
    #assign(paste0("mr_orig",i), lm(formula2,data=get(paste0("test_data",i))))
    assign(paste0("modelSumm_orig",i), summary(get(paste0("mr_orig",i))))
    #assign(paste0("m.out",i), matchit(PSformulaUnWt1,data=get(paste0("test_data",i)),replace=TRUE, method= "nearest", distance = "logit"))
    #assign(paste0("m.out",i), matchit(formula1,data=get(paste0("test_data",i)),replace=TRUE, method= "nearest", distance = "logit"))
    #assign(paste0("m.data1",i), match.data(get(paste0("m.out",i))))
    #assign(paste0("mr",i), lm(formula2, weights=get(paste0("m.data1",i))$weights, data=get(paste0("m.data1",i))))
    #assign(paste0("mr",i), lm(formula2, data=get(paste0("m.data1",i))))
    #assign(paste0("modelSumm",i), summary(get(paste0("mr",i))))
    assign(paste0("test_data",i), cbind(idVar,get(paste0("test_data",i))))
  }
  test_dataFour<<-test_data4
  test_dataFive<<-test_data5
  for (i in 1:Sims){
    #assign(paste0("m.out",i), matchit(formula1,data=get(paste0("test_data",i)), method= "subclass", suclass=5))
    #assign(paste0("m.out",i), matchit(PSformulaUnWt1,data=get(paste0("test_data",i)), method= "subclass", subclass=5))
    assign(paste0("ps.binary",i), ps(formula4,data=get(paste0("test_data",i)), estimand= 'ATE', stop.method="es.max",n.trees=500, verbose=FALSE))
    #assign(paste0("z.out",i), zelig(formula3,data=match.data(get(paste0("m.out",i)),"control"),model="ls",by = "subclass", cite=FALSE))
    #assign(paste0("x.out",i), setx(get(paste0("z.out",i)),data=match.data(get(paste0("m.out",i)),"treat"),fn=NULL, cond=TRUE))
    #assign(paste0("x1.out",i), setx(get(paste0("z.out",i)),treat=1))
    #assign(paste0("s.out",i), sim(get(paste0("z.out",i)),x=get(paste0("x.out",i))))
    #assign(paste0("x",i), unlist(get(paste0("s.out",i))[["sim.out"]][["x"]][["ev"]]))
  }
  
  for (i in 1:Sims){
    if(i >= 1 & i <= Sims) {
      counter = i
      print(counter)
    }
    assign(paste0("w",i),get.weights(get(paste0("ps.binary",i)),stop.method = "es.max"))
    assign(paste0("test_data",i), cbind(get(paste0("test_data",i)),get(paste0("w",i))))
    assign(paste0("design.ps",i),svydesign(ids=~1,weights = get(paste0("w",i)),data=get(paste0("test_data",i))))
    assign(paste0("glmOne",i),svyglm(Y~treat,design=get(paste0("design.ps",i))))
    assign(paste0("summglmOne",i),summary(get(paste0("glmOne",i))))
  }
  ##create formula5 from covariates in GBM propensity score model
  
  form3<-gsub("ESCS \\+ ","",PSformulaUnWt1)
  newvar<-gsub("ESCS:SC048Q01 \\+ ","",form3[3])
  newvar<-gsub("SC048Q01:ESCS \\+ ","",newvar)
  newvar<-gsub("ST011Q07:ESCS \\+ ","",newvar)
  newvar<-gsub("BELONG:ESCS \\+ ","",newvar)
  newvar<-gsub("ESCS:ST011Q07 \\+ ","",newvar)
  newvar<-gsub("ESCS:BELONG \\+ ","",newvar)
  newvar<-gsub("ESCS:ST118Q04 \\+ ","",newvar)
  newvar<-gsub("ESCS:MISCED \\+ ","",newvar)
  newvar<-gsub("ESCS:ST013Q01 \\+ ","",newvar)
  newvar<-gsub("\\+ ESCS:ST011Q07","",newvar)
  formula5<-as.formula(paste("ESCS ~", newvar))
  
  #######Sensitivity Analysis##############################
  
  test_dataGLM<-cbind(test_data5, ps.binary1$w)
  sen_control<<-subset(test_dataGLM,test_dataGLM$treat==0)
  sen_treat<<-subset(test_dataGLM, test_dataGLM$treat==1)
  num_cont<<-dim(sen_control)[1]
  num_treat<<-dim(sen_treat)[1]
  design.control<<-svydesign(ids = ~1, weights=sen_control$es.max.ATE,data=sen_control)
  design.treat<<-svydesign(ids = ~1, weights=sen_treat$es.max.ATE,data=sen_treat)
  lmSen_cont<<-svyglm(formula5, design=design.control)
  lmSen_treat<<-svyglm(formula5, design=design.treat)
  summ_cont<<-summary(lmSen_cont)
  summ_treat<<-summary(lmSen_treat)
  indx<<-(length(summ_cont$coefficients)/4)+1
  indx2<<-(length(summ_treat$coefficients)/4)+1
  se_1<<-summ_cont$coefficients[indx]
  se_2<<-summ_treat$coefficients[indx2]
  beta_1<<-summ_cont$coefficients[1]
  beta_2<<-summ_treat$coefficients[1]
  t_beta=(beta_1 - beta_2)/sqrt((se_1)^2 + (se_2)^2)
  nu=((((se_1)^2/num_cont)+((se_2)^2/num_treat))^2)/
    ((1/(num_cont-1))*(((se_1)^2/num_cont)^2)+(1/(num_treat-1))*(((se_2)^2/num_treat)^2))
  
  ############END of SENSITIVITY ANALYSIS#####################
  
  ps.binaryOne<<-ps.binary1
  ps.binaryTwo<<-ps.binary2
  predNamesOne<<-predNames
  predNamesTwo<<-predNames2
  formulaFour<<-formula4
  #z.outOne<<-z.out1
  #X.outOne<<-x.out1
  #s.outfour<<-s.out4
  #s.outthree<<-s.out3
  #View(z.out2$data)  
  ##Average the percent reduction in mean differences in the covariates before and after matching. Create table
  ##with nrow = to number of simulations and ncol=the numbers of variables in the model
  GeneralizedBoostedPSTable=matrix(rep(0,Sims*2),nrow=Sims, ncol=2)
  colnames(GeneralizedBoostedPSTable)=c("treat","SD")
  
  
  for (i in 1:Sims){
    GeneralizedBoostedPSTable[i,1]=get(paste0("summglmOne",i))$coefficients[2]
    GeneralizedBoostedPSTable[i,2]=get(paste0("summglmOne",i))$coefficients[4]
  }
  
  
  
  View(GeneralizedBoostedPSTable)
  
  TFrame<-as.data.frame(rep(0,3),ncol=3)
  TFrame[1,1]=mean(GeneralizedBoostedPSTable[,1])
  TFrame[1,2]=mean(GeneralizedBoostedPSTable[,2])
  TFrame[1,3]=2-TFrame[1,1]
  
  GeneralizedBoostedPSTable1<<-GeneralizedBoostedPSTable
  
  names(TFrame)<-c("MeanExpected","MeanSD","bias")
  
  View(TFrame)
  
  W<-sprintf("The mean of the treatment effect with GBM-propensiy score estimated model is %f, the standard deviation is %f, and the bias is %f", TFrame[1,1],TFrame[1,2],TFrame[1,3])
  
  print(W)
  
  W3<-sprintf("The t-test for the sensitivity analysis for matched units is t(%.0f)= %.2f",nu,t_beta)
  print(W3)
}
