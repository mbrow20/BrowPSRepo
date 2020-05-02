##This is provides the variances estimator for Imbens and Rubin (2015) for 500 simulated data sets
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
  source("Functions5VAR500.R")
  S2<<-as.matrix(CovMatrixWt,nrow=nrow(CovMatrixWt),ncol=ncol(CovMatrixWt))
  #S2<<-as.matrix(CovMatrixWt,nrow=5,ncol=5)
  set.seed(4)#ensures same results for random components (e.g., partitioning data, random variable selection, etc.)
  colnames(data1)[1]<<-"SC048Q01"
  dpart<<-createDataPartition(data1$treat,p=0.2,list=F)
  data2<<-data1[dpart,]
  #data2<<-data1
  varNames<<-colnames(data2)
  n=nrow(data2)
  Q<<-length(unique(data2$STRATUM)) #Number of Explicit strata; Should be less than 'n'
  g2<<-as.data.frame(rep(0,n),ncol=1)
  g2[,1]<<-data2$STRATUM
  #These are actual Strata indicators used in the PISA data. Each stratum is composed is composed
  #of three explicit stratification variables: Province, Language, and School size (see OECD 2017, p. 72).
  g2$SCHID<<-data2$SCHID#These school ids correspond to actual schools. Make sure schools in dataset correspond to a variable labeled 'SCHID'.
  g2<<-g2[order(g2[,1]),]
  colnames(g2)<<-c("g","SCHID")
  lis<<-as.data.frame(unique(g2$g))#There are 46 unique schools in this sample
  lis2<<-as.vector(lis[,1])
  lis3<<-as.character(lis2)
  lis4<<-as.character(seq(1,length(lis3),1))#Schools are re-labeled in sequential order. The 'mapvalues' function requires class 'character'
  g3<<-mapvalues(g2$g,c(lis3),c(lis4))#an case in the first vector that matches a case in the second vector will be replace by item in third vector
  g2<<-cbind(as.data.frame(as.numeric(g3)),g2$SCHID)
  colnames(g2)<<-c("g","SCHID")
  listSCH1<<-as.data.frame(unique(g2$SCHID))
  listSCH2<<-as.vector(listSCH1[,1])
  listSCH3<<-as.character(listSCH2)
  listSCH4<<-as.character(seq(1,length(listSCH3),1))
  g4<<-mapvalues(g2$SCHID,c(listSCH3),c(listSCH4))
  g2<<-cbind(as.data.frame(as.numeric(g2[,1])),as.data.frame(as.numeric(g4)))
  colnames(g2)<<-c("g","SCHID")
  ##Create weights for implicit stratification (corresponding to the three levels of urbanicity) within the explicit strata
  for (i in 1:nrow(g2)){
    for (j in 1:Q){
      g2[i,j+2]<<-ifelse(g2$g[i]==j,sample((1:3)/4,1,replace=TRUE),0)}}
  ExplicitWT<<-rgamma(length(unique(g2$g)),1,1)
  g2$ExplWt<<-ExplicitWT[g2$g]
  SCHWT<<-rgamma(unique(g2$SCHID),1,1)
  g2$SCHWT<<-SCHWT[g2$SCHID]
  g2$FSCHWT<<-rowSums(g2[,c(3:ncol(g2))])
  g2$scaleweights<<-g2$FSCHWT*10#increase weights by a factor of ten to correspond to actual PISA weights
  g2$Female<<-data2$Female##As stated earlier, 'female' is one of the cells in the within school non-response; therefore, it is included as weight.
  g2$FSTUWT<<-g2$Female+g2$scaleweights#the weight of the 'female' variable is 1.
  g2$CNTSTUID<<-data2$CNTSTUID
  g2<<-g2[order(g2$SCHID),]
  colnames(data2)[1]="SC048Q01"
  #mu_beta<<-c(1.1,2,-1,1,1.5,0.04,1.2,-1.3,0.8,-0.5,-1.4)##The actual population beta coefficient is '2.0'. Recall, we are using wieghts, so inferences are made to the population.
  K<<-length(mu_beta)-2
  #make sigma^2 from inverse gamma distribution with parameters 1 and 1 BayesM package
  sigma_2<<-rgamma(nrow(data2),1,1)
  #sigma_2<<-1
  sigma_2<<-sigma_2/g2$FSTUWT
  v<<-length(mu_beta)
  varNames<<-varNames[c(1:v-1)]
  Sigma_beta<<-riwish(v,S2)
  Beta<<-mvrnorm(1,mu_beta,Sigma_beta)
  U<<-mvrnorm(length(unique(g2$SCHID)),rep(0,K+2),diag(K+2))
  U<<-U[g2$SCHID,]
  length(unique(U[,1]))
  Ones<<-matrix(1,nrow(data2),1)
  Z<<-as.matrix(data2[,1:v-1])
  #Z<<-as.matrix(data2[,1:4])
  Z<<-cbind(Ones,Z)
  Y=Z %*% Beta + diag(Z %*% t(U)) + rnorm(n,0,sqrt(sigma_2))#One Simulated data of Y
  test_data<<-cbind(Z,Y,g2$FSTUWT,g2$SCHID)
  varNames[1]="SC048Q01"
  colnames(test_data)=c("V1",varNames, "Y", "FSTUWT", "SCHID")
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
  LinModel<-lm(formulaRegUnWt,data=test_data, weights=FSTUWT)
  LinModel<<-lm(formulaRegUnWt,data=test_data, weights=FSTUWT)
  print(summary(LinModel))
  print(confint(LinModel))
  
  Sims=Sims
  R<-as.data.frame(rep(0,Sims),ncol=1)
  Beta_1<-matrix(rep(0,Sims*(K+2)),Sims,v)
  Beta_1<-as.data.frame(Beta_1)
  Beta=mvrnorm(Sims,mu_beta,Sigma_beta)
  int<-NULL
  mr<-NULL
  TFrame<-NULL
  indxer2=which(colnames(Z)=="treat" )
  Zt<-Z[,-c(1,indxer2)]
  
  for(i in 1:Sims){
    assign(paste0("Y",i), Z %*% Beta[i,] + diag(Z %*% t(U))+rnorm(nrow(Z),0,sqrt(sigma_2)[i]))
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
    assign(paste0("PSformulaUnWt",1), as.formula(paste("treat~",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+")))))
    #assign(paste0("mr_orig",i), lm(formula2,data=get(paste0("test_data",i)),weights = get(paste0("test_data",i))$FSTUWT))
    assign(paste0("mr_orig",i), lm(formula2,data=get(paste0("test_data",i))))
    assign(paste0("modelSumm_orig",i), summary(get(paste0("mr_orig",i))))
    assign(paste0("modelSE_orig", i) , get(paste0("modelSumm_orig",i))[[4]][[18]])
    
    #assign(paste0("m.out",i), matchit(PSformulaUnWt1,data=get(paste0("test_data",i)),replace=TRUE, method= "nearest", distance = "logit"))
    #assign(paste0("m.out",i), matchit(formula1,data=get(paste0("test_data",i)),replace=TRUE, method= "nearest", distance = "logit"))
    #assign(paste0("m.data1",i), match.data(get(paste0("m.out",i))))
    #assign(paste0("mr",i), lm(formula2, weights=get(paste0("m.data1",i))$weights, data=get(paste0("m.data1",i))))
    #assign(paste0("mr",i), lm(formula2, data=get(paste0("m.data1",i))))
    #assign(paste0("modelSumm",i), summary(get(paste0("mr",i))))
  }
  PSformulaUnWt1<<-PSformulaUnWt1
  formula2<<-formula2
  
  for (i in 1:Sims){
    #assign(paste0("m.out",i), matchit(formula1,data=get(paste0("test_data",i)), method= "nearest"))
    assign(paste0("m.out",i), matchit(PSformulaUnWt1,data=get(paste0("test_data",i)), method= "nearest"))
    assign(paste0("z.out",i), zelig(formula2,data=match.data(get(paste0("m.out",i))),model="ls",cite = FALSE))
    assign(paste0("x.out",i), setx(get(paste0("z.out",i)),treat=0))
    assign(paste0("x1.out",i), setx(get(paste0("z.out",i)),treat=1))
    assign(paste0("s.out",i), sim(get(paste0("z.out",i)),x=get(paste0("x.out",i)),x1=get(paste0("x1.out",i))))
    assign(paste0("modelSumm",i), coef(get(paste0("s.out",i)))[2])
    assign(paste0("modelSE", i) , get_se(get(paste0("z.out",i)))[[1]][[2]])
    ###
    assign("summMatched1", summary(get(paste0("m.out",i)))$reduction$'Mean Diff.')
    assign("summMatched2", summary(get(paste0("m.out",i)))$reduction$'eQQ Med')
    assign("summMatched3", summary(get(paste0("m.out",i)))$reduction$'eQQ Mean')
    assign("summMatched4", summary(get(paste0("m.out",i)))$reduction$'eQQ Max')
  }
  View(z.out2$data)  
  s.outtwo<<-s.out2
  m.outtwo<<-m.out2
  z.outtwo<<-z.out2
  Sims<<-Sims
  
  m.outx1<<-m.out1
  m.outx2<<-m.out2
  m.outx3<<-m.out3
  m.outx4<<-m.out4
  m.outx5<<-m.out5
  m.outx6<<-m.out6
  m.outx7<<-m.out7
  m.outx8<<-m.out8
  m.outx9<<-m.out9
  m.outx10<<-m.out10
  m.outx11<<-m.out11
  m.outx12<<-m.out12
  m.outx13<<-m.out13
  m.outx14<<-m.out14
  m.outx15<<-m.out15
  m.outx16<<-m.out16
  m.outx17<<-m.out17
  m.outx18<<-m.out18
  m.outx19<<-m.out19
  m.outx20<<-m.out20
  m.outx21<<-m.out21
  m.outx22<<-m.out22
  m.outx23<<-m.out23
  m.outx24<<-m.out24
  m.outx25<<-m.out25
  m.outx26<<-m.out26
  m.outx27<<-m.out27
  m.outx28<<-m.out28
  m.outx29<<-m.out29
  m.outx30<<-m.out30
  m.outx31<<-m.out31
  m.outx32<<-m.out32
  m.outx33<<-m.out33
  m.outx34<<-m.out34
  m.outx35<<-m.out35
  m.outx36<<-m.out36
  m.outx37<<-m.out37
  m.outx38<<-m.out38
  m.outx39<<-m.out39
  m.outx40<<-m.out40
  m.outx41<<-m.out41
  m.outx42<<-m.out42
  m.outx43<<-m.out43
  m.outx44<<-m.out44
  m.outx45<<-m.out45
  m.outx46<<-m.out46
  m.outx47<<-m.out47
  m.outx48<<-m.out48
  m.outx49<<-m.out49
  m.outx50<<-m.out50
  m.outx51<<-m.out51
  m.outx52<<-m.out52
  m.outx53<<-m.out53
  m.outx54<<-m.out54
  m.outx55<<-m.out55
  m.outx56<<-m.out56
  m.outx57<<-m.out57
  m.outx58<<-m.out58
  m.outx59<<-m.out59
  m.outx60<<-m.out60
  m.outx61<<-m.out61
  m.outx62<<-m.out62
  m.outx63<<-m.out63
  m.outx64<<-m.out64
  m.outx65<<-m.out65
  m.outx66<<-m.out66
  m.outx67<<-m.out67
  m.outx68<<-m.out68
  m.outx69<<-m.out69
  m.outx70<<-m.out70
  m.outx71<<-m.out71
  m.outx72<<-m.out72
  m.outx73<<-m.out73
  m.outx74<<-m.out74
  m.outx75<<-m.out75
  m.outx76<<-m.out76
  m.outx77<<-m.out77
  m.outx78<<-m.out78
  m.outx79<<-m.out79
  m.outx80<<-m.out80
  m.outx81<<-m.out81
  m.outx82<<-m.out82
  m.outx83<<-m.out83
  m.outx84<<-m.out84
  m.outx85<<-m.out85
  m.outx86<<-m.out86
  m.outx87<<-m.out87
  m.outx88<<-m.out88
  m.outx89<<-m.out89
  m.outx90<<-m.out90
  m.outx91<<-m.out91
  m.outx92<<-m.out92
  m.outx93<<-m.out93
  m.outx94<<-m.out94
  m.outx95<<-m.out95
  m.outx96<<-m.out96
  m.outx97<<-m.out97
  m.outx98<<-m.out98
  m.outx99<<-m.out99
  m.outx100<<-m.out100
  m.outx101<<-m.out101
  m.outx102<<-m.out102
  m.outx103<<-m.out103
  m.outx104<<-m.out104
  m.outx105<<-m.out105
  m.outx106<<-m.out106
  m.outx107<<-m.out107
  m.outx108<<-m.out108
  m.outx109<<-m.out109
  m.outx110<<-m.out110
  m.outx111<<-m.out111
  m.outx112<<-m.out112
  m.outx113<<-m.out113
  m.outx114<<-m.out114
  m.outx115<<-m.out115
  m.outx116<<-m.out116
  m.outx117<<-m.out117
  m.outx118<<-m.out118
  m.outx119<<-m.out119
  m.outx120<<-m.out120
  m.outx121<<-m.out121
  m.outx122<<-m.out122
  m.outx123<<-m.out123
  m.outx124<<-m.out124
  m.outx125<<-m.out125
  m.outx126<<-m.out126
  m.outx127<<-m.out127
  m.outx128<<-m.out128
  m.outx129<<-m.out129
  m.outx130<<-m.out130
  m.outx131<<-m.out131
  m.outx132<<-m.out132
  m.outx133<<-m.out133
  m.outx134<<-m.out134
  m.outx135<<-m.out135
  m.outx136<<-m.out136
  m.outx137<<-m.out137
  m.outx138<<-m.out138
  m.outx139<<-m.out139
  m.outx140<<-m.out140
  m.outx141<<-m.out141
  m.outx142<<-m.out142
  m.outx143<<-m.out143
  m.outx144<<-m.out144
  m.outx145<<-m.out145
  m.outx146<<-m.out146
  m.outx147<<-m.out147
  m.outx148<<-m.out148
  m.outx149<<-m.out149
  m.outx150<<-m.out150
  m.outx151<<-m.out151
  m.outx152<<-m.out152
  m.outx153<<-m.out153
  m.outx154<<-m.out154
  m.outx155<<-m.out155
  m.outx156<<-m.out156
  m.outx157<<-m.out157
  m.outx158<<-m.out158
  m.outx159<<-m.out159
  m.outx160<<-m.out160
  m.outx161<<-m.out161
  m.outx162<<-m.out162
  m.outx163<<-m.out163
  m.outx164<<-m.out164
  m.outx165<<-m.out165
  m.outx166<<-m.out166
  m.outx167<<-m.out167
  m.outx168<<-m.out168
  m.outx169<<-m.out169
  m.outx170<<-m.out170
  m.outx171<<-m.out171
  m.outx172<<-m.out172
  m.outx173<<-m.out173
  m.outx174<<-m.out174
  m.outx175<<-m.out175
  m.outx176<<-m.out176
  m.outx177<<-m.out177
  m.outx178<<-m.out178
  m.outx179<<-m.out179
  m.outx180<<-m.out180
  m.outx181<<-m.out181
  m.outx182<<-m.out182
  m.outx183<<-m.out183
  m.outx184<<-m.out184
  m.outx185<<-m.out185
  m.outx186<<-m.out186
  m.outx187<<-m.out187
  m.outx188<<-m.out188
  m.outx189<<-m.out189
  m.outx190<<-m.out190
  m.outx191<<-m.out191
  m.outx192<<-m.out192
  m.outx193<<-m.out193
  m.outx194<<-m.out194
  m.outx195<<-m.out195
  m.outx196<<-m.out196
  m.outx197<<-m.out197
  m.outx198<<-m.out198
  m.outx199<<-m.out199
  m.outx200<<-m.out200
  m.outx201<<-m.out201
  m.outx202<<-m.out202
  m.outx203<<-m.out203
  m.outx204<<-m.out204
  m.outx205<<-m.out205
  m.outx206<<-m.out206
  m.outx207<<-m.out207
  m.outx208<<-m.out208
  m.outx209<<-m.out209
  m.outx210<<-m.out210
  m.outx211<<-m.out211
  m.outx212<<-m.out212
  m.outx213<<-m.out213
  m.outx214<<-m.out214
  m.outx215<<-m.out215
  m.outx216<<-m.out216
  m.outx217<<-m.out217
  m.outx218<<-m.out218
  m.outx219<<-m.out219
  m.outx220<<-m.out220
  m.outx221<<-m.out221
  m.outx222<<-m.out222
  m.outx223<<-m.out223
  m.outx224<<-m.out224
  m.outx225<<-m.out225
  m.outx226<<-m.out226
  m.outx227<<-m.out227
  m.outx228<<-m.out228
  m.outx229<<-m.out229
  m.outx230<<-m.out230
  m.outx231<<-m.out231
  m.outx232<<-m.out232
  m.outx233<<-m.out233
  m.outx234<<-m.out234
  m.outx235<<-m.out235
  m.outx236<<-m.out236
  m.outx237<<-m.out237
  m.outx238<<-m.out238
  m.outx239<<-m.out239
  m.outx240<<-m.out240
  m.outx241<<-m.out241
  m.outx242<<-m.out242
  m.outx243<<-m.out243
  m.outx244<<-m.out244
  m.outx245<<-m.out245
  m.outx246<<-m.out246
  m.outx247<<-m.out247
  m.outx248<<-m.out248
  m.outx249<<-m.out249
  m.outx250<<-m.out250
  m.outx251<<-m.out251
  m.outx252<<-m.out252
  m.outx253<<-m.out253
  m.outx254<<-m.out254
  m.outx255<<-m.out255
  m.outx256<<-m.out256
  m.outx257<<-m.out257
  m.outx258<<-m.out258
  m.outx259<<-m.out259
  m.outx260<<-m.out260
  m.outx261<<-m.out261
  m.outx262<<-m.out262
  m.outx263<<-m.out263
  m.outx264<<-m.out264
  m.outx265<<-m.out265
  m.outx266<<-m.out266
  m.outx267<<-m.out267
  m.outx268<<-m.out268
  m.outx269<<-m.out269
  m.outx270<<-m.out270
  m.outx271<<-m.out271
  m.outx272<<-m.out272
  m.outx273<<-m.out273
  m.outx274<<-m.out274
  m.outx275<<-m.out275
  m.outx276<<-m.out276
  m.outx277<<-m.out277
  m.outx278<<-m.out278
  m.outx279<<-m.out279
  m.outx280<<-m.out280
  m.outx281<<-m.out281
  m.outx282<<-m.out282
  m.outx283<<-m.out283
  m.outx284<<-m.out284
  m.outx285<<-m.out285
  m.outx286<<-m.out286
  m.outx287<<-m.out287
  m.outx288<<-m.out288
  m.outx289<<-m.out289
  m.outx290<<-m.out290
  m.outx291<<-m.out291
  m.outx292<<-m.out292
  m.outx293<<-m.out293
  m.outx294<<-m.out294
  m.outx295<<-m.out295
  m.outx296<<-m.out296
  m.outx297<<-m.out297
  m.outx298<<-m.out298
  m.outx299<<-m.out299
  m.outx300<<-m.out300
  m.outx301<<-m.out301
  m.outx302<<-m.out302
  m.outx303<<-m.out303
  m.outx304<<-m.out304
  m.outx305<<-m.out305
  m.outx306<<-m.out306
  m.outx307<<-m.out307
  m.outx308<<-m.out308
  m.outx309<<-m.out309
  m.outx310<<-m.out310
  m.outx311<<-m.out311
  m.outx312<<-m.out312
  m.outx313<<-m.out313
  m.outx314<<-m.out314
  m.outx315<<-m.out315
  m.outx316<<-m.out316
  m.outx317<<-m.out317
  m.outx318<<-m.out318
  m.outx319<<-m.out319
  m.outx320<<-m.out320
  m.outx321<<-m.out321
  m.outx322<<-m.out322
  m.outx323<<-m.out323
  m.outx324<<-m.out324
  m.outx325<<-m.out325
  m.outx326<<-m.out326
  m.outx327<<-m.out327
  m.outx328<<-m.out328
  m.outx329<<-m.out329
  m.outx330<<-m.out330
  m.outx331<<-m.out331
  m.outx332<<-m.out332
  m.outx333<<-m.out333
  m.outx334<<-m.out334
  m.outx335<<-m.out335
  m.outx336<<-m.out336
  m.outx337<<-m.out337
  m.outx338<<-m.out338
  m.outx339<<-m.out339
  m.outx340<<-m.out340
  m.outx341<<-m.out341
  m.outx342<<-m.out342
  m.outx343<<-m.out343
  m.outx344<<-m.out344
  m.outx345<<-m.out345
  m.outx346<<-m.out346
  m.outx347<<-m.out347
  m.outx348<<-m.out348
  m.outx349<<-m.out349
  m.outx350<<-m.out350
  m.outx351<<-m.out351
  m.outx352<<-m.out352
  m.outx353<<-m.out353
  m.outx354<<-m.out354
  m.outx355<<-m.out355
  m.outx356<<-m.out356
  m.outx357<<-m.out357
  m.outx358<<-m.out358
  m.outx359<<-m.out359
  m.outx360<<-m.out360
  m.outx361<<-m.out361
  m.outx362<<-m.out362
  m.outx363<<-m.out363
  m.outx364<<-m.out364
  m.outx365<<-m.out365
  m.outx366<<-m.out366
  m.outx367<<-m.out367
  m.outx368<<-m.out368
  m.outx369<<-m.out369
  m.outx370<<-m.out370
  m.outx371<<-m.out371
  m.outx372<<-m.out372
  m.outx373<<-m.out373
  m.outx374<<-m.out374
  m.outx375<<-m.out375
  m.outx376<<-m.out376
  m.outx377<<-m.out377
  m.outx378<<-m.out378
  m.outx379<<-m.out379
  m.outx380<<-m.out380
  m.outx381<<-m.out381
  m.outx382<<-m.out382
  m.outx383<<-m.out383
  m.outx384<<-m.out384
  m.outx385<<-m.out385
  m.outx386<<-m.out386
  m.outx387<<-m.out387
  m.outx388<<-m.out388
  m.outx389<<-m.out389
  m.outx390<<-m.out390
  m.outx391<<-m.out391
  m.outx392<<-m.out392
  m.outx393<<-m.out393
  m.outx394<<-m.out394
  m.outx395<<-m.out395
  m.outx396<<-m.out396
  m.outx397<<-m.out397
  m.outx398<<-m.out398
  m.outx399<<-m.out399
  m.outx400<<-m.out400
  m.outx401<<-m.out401
  m.outx402<<-m.out402
  m.outx403<<-m.out403
  m.outx404<<-m.out404
  m.outx405<<-m.out405
  m.outx406<<-m.out406
  m.outx407<<-m.out407
  m.outx408<<-m.out408
  m.outx409<<-m.out409
  m.outx410<<-m.out410
  m.outx411<<-m.out411
  m.outx412<<-m.out412
  m.outx413<<-m.out413
  m.outx414<<-m.out414
  m.outx415<<-m.out415
  m.outx416<<-m.out416
  m.outx417<<-m.out417
  m.outx418<<-m.out418
  m.outx419<<-m.out419
  m.outx420<<-m.out420
  m.outx421<<-m.out421
  m.outx422<<-m.out422
  m.outx423<<-m.out423
  m.outx424<<-m.out424
  m.outx425<<-m.out425
  m.outx426<<-m.out426
  m.outx427<<-m.out427
  m.outx428<<-m.out428
  m.outx429<<-m.out429
  m.outx430<<-m.out430
  m.outx431<<-m.out431
  m.outx432<<-m.out432
  m.outx433<<-m.out433
  m.outx434<<-m.out434
  m.outx435<<-m.out435
  m.outx436<<-m.out436
  m.outx437<<-m.out437
  m.outx438<<-m.out438
  m.outx439<<-m.out439
  m.outx440<<-m.out440
  m.outx441<<-m.out441
  m.outx442<<-m.out442
  m.outx443<<-m.out443
  m.outx444<<-m.out444
  m.outx445<<-m.out445
  m.outx446<<-m.out446
  m.outx447<<-m.out447
  m.outx448<<-m.out448
  m.outx449<<-m.out449
  m.outx450<<-m.out450
  m.outx451<<-m.out451
  m.outx452<<-m.out452
  m.outx453<<-m.out453
  m.outx454<<-m.out454
  m.outx455<<-m.out455
  m.outx456<<-m.out456
  m.outx457<<-m.out457
  m.outx458<<-m.out458
  m.outx459<<-m.out459
  m.outx460<<-m.out460
  m.outx461<<-m.out461
  m.outx462<<-m.out462
  m.outx463<<-m.out463
  m.outx464<<-m.out464
  m.outx465<<-m.out465
  m.outx466<<-m.out466
  m.outx467<<-m.out467
  m.outx468<<-m.out468
  m.outx469<<-m.out469
  m.outx470<<-m.out470
  m.outx471<<-m.out471
  m.outx472<<-m.out472
  m.outx473<<-m.out473
  m.outx474<<-m.out474
  m.outx475<<-m.out475
  m.outx476<<-m.out476
  m.outx477<<-m.out477
  m.outx478<<-m.out478
  m.outx479<<-m.out479
  m.outx480<<-m.out480
  m.outx481<<-m.out481
  m.outx482<<-m.out482
  m.outx483<<-m.out483
  m.outx484<<-m.out484
  m.outx485<<-m.out485
  m.outx486<<-m.out486
  m.outx487<<-m.out487
  m.outx488<<-m.out488
  m.outx489<<-m.out489
  m.outx490<<-m.out490
  m.outx491<<-m.out491
  m.outx492<<-m.out492
  m.outx493<<-m.out493
  m.outx494<<-m.out494
  m.outx495<<-m.out495
  m.outx496<<-m.out496
  m.outx497<<-m.out497
  m.outx498<<-m.out498
  m.outx499<<-m.out499
  m.outx500<<-m.out500
  
  
  sumMatchedOne<<-summMatched1
  sumMatchedTwo<<-summMatched2
  sumMatchedThree<<-summMatched3
  sumMatchedFour<<-summMatched4
  
  PercMeanDiffTable=matrix(rep(0,length(names(m.out2$model$coefficients))*4),nrow=4, ncol=length(names(m.out2$model$coefficients)))
  colnames(PercMeanDiffTable)=names(m.out2$model$coefficients)
  rownames(PercMeanDiffTable)=c("MeanDiff", "eQQMedian", "eQQMean", "eQQMax")
  
  
  PercMeanDiffTable[1,]=summMatched1
  PercMeanDiffTable[2,]=summMatched2
  PercMeanDiffTable[3,]=summMatched3
  PercMeanDiffTable[4,]=summMatched4
  View(PercMeanDiffTable)
  
  TFrame<-as.data.frame(rep(0,Sims),ncol=4)
  for (i in 1:Sims){
    TFrame[i,1]<-get(paste0("modelSumm",i))#$coefficients[[2]]
    TFrame[i,2]<-get(paste0("modelSE",i))#SE for matched
    TFrame[i,3]<-get(paste0("modelSumm_orig",i))$coefficients[[2]]
    TFrame[i,4]<-get(paste0("modelSE_orig",i))}
  names(TFrame)<-c("matched","matchedSE", "unmatched","unmatchedSE")
  View(TFrame)
  for(i in 1:nrow(TFrame)){
    TFrame$m.bias[i]=2-TFrame$matched[i]
    TFrame$um.bias[i]=2-TFrame$unmatched[i]}
  View(TFrame)
  range.um.bias<-range(TFrame$um.bias)
  range.m.bias<-range(TFrame$m.bias)
  W<-sprintf("The mean bias of unmatched units is %f with a mean SE of %f, and for matched units is %f with a mean SE of %f", mean(TFrame$um.bias),mean(TFrame$unmatchedSE), mean(TFrame$m.bias), mean(TFrame$matchedSE))
  W1<-sprintf("The range of the unmatched bias estimates is from %f to %f, matched is from %f to %f", range.um.bias[[1]], range.um.bias[[2]], range.m.bias[[1]], range.m.bias[[2]])
  print(W)
  print(W1)
  print(TFrame)
  PercMeanDiffTable<<-PercMeanDiffTable
  TFrameVAR<<-TFrame
  
}

