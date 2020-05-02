##Functions5.R is for the 15 Variable 2 Level Simulation##
PS.CovSelection.Function<-function(data4){
  start_time = Sys.time()
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
    assign(paste0("formula",1), as.formula(paste("treated~", paste(factors[i], collapse="+"))))
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
        assign(paste0("formula",j), as.formula(paste("treated~", paste(y,paste("+"),paste(factors[i], collapse="+")))))
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
  FormulaNew=as.formula(paste("treated", z, sep=" ~ "))
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
  formula999<-paste("treated~",Orig.Var)
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
  
  PSformulaUnWt<<-as.formula(paste("treated~",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+"))))
  formulaRegUnWt<<-as.formula(paste("Y~treated+",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+"))))
  end_time = Sys.time()
  Difference = end_time-start_time
  print(Difference)
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
    assign(paste0("formula",1), as.formula(paste("treated~", paste(factors[i], collapse="+"))))
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
        assign(paste0("formula",j), as.formula(paste("treated~", paste(y,paste("+"),paste(factors[i], collapse="+")))))
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
  FormulaNew=as.formula(paste("treated", z, sep=" ~ "))
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
  formula999<-paste("treated~",Orig.Var)
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
  PSformulaWt<<-as.formula(paste("treated~",paste(linearVarsWt,collapse="+"),paste("+",paste(interActionVarsWt,collapse="+"))))
  formulaRegWt<<-as.formula(paste("Y~treated+",paste(linearVarsWt,collapse="+"),paste("+",paste(interActionVarsWt,collapse="+"))))
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
  source("Functions5SR_2020CEM.R")
  numCores <<- detectCores()
  S2<<-as.matrix(CovMatrixWt,nrow=nrow(CovMatrixWt),ncol=ncol(CovMatrixWt))
  #S2<<-as.matrix(CovMatrixWt,nrow=5,ncol=5)
  set.seed(4)#ensures same results for random components (e.g., partitioning data, random variable selection, etc.)
  colnames(data1)[1]<<-"SC048Q01"
  #dpart<<-createDataPartition(data1$treated,p=0.6,list=F)
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
  g2$W_FSTUWT2<<-data2$W_FSTUWT
  #g2<<-g2[order(g2$SCHID),]
  colnames(data2)[1]="SC048Q01"
  #mu_beta<<-c(1.1,2,-1,1,1.5,0.04,1.2,-1.3,0.8,-0.5,-1.4)##The actual population beta coefficient is '2.0'. Recall, we are using wieghts, so inferences are made to the population.
  K<<-length(mu_beta)-2
  #make sigma^2 from inverse gamma distribution with parameters 1 and 1 BayesM package
  sigma_2<<-rgamma(nrow(data2),1,1)
  #sigma_2<<-1
  sigma_2<<-sigma_2/g2$W_FSTUWT2
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
  test_data<<-cbind(Z,Y,g2$FSTUWT,data2$SCHID)
  varNames[1]="SC048Q01"
  colnames(test_data)=c("V1",varNames, "Y", "FSTUWT", "SCHID")
  test_data<<-as.data.frame(test_data)
  write.table(test_data,file="test_data_NEW.csv",row.names=FALSE,col.names=TRUE, sep=",",append=FALSE)
  ##Examine the 'beta' of the treatment effect through a linear model
  data4<<-test_data[,c(1:length(varNames)+1)]
  data4<<-as.data.frame(data4)
  test_data<<-as.data.frame(test_data)
  indxer=which( colnames(data4)=="treated" )#--in case yourdataset has the variable 'treat' in a different column than one
  
  ###The following functions reorder 'data4' so the first variable is 'treat'##########################
  data4NotIndxer<<-data4[,-(indxer)]
  data4Indxer<<-data4[,indxer]
  data4<<-cbind(data4Indxer,data4NotIndxer)
  data4<<-as.data.frame(data4)
  colnames(data4)[1]="treated"
  data4<<-as.data.frame(data4)
  
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
  Beta=mvrnorm(Sims,mu_beta,Sigma_beta)
  int<-NULL
  mr<-NULL
  TFrame<-NULL
  indxer2=which(colnames(Z)=="treated" )
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
    assign(colnames(get(paste0("predVars",i)))[1],"Treated")
    assign(colnames(get(paste0("test_data",i)))[ncol(Z)],"Y")
    nm<-paste0("predVars",i)
    #nm<-paste0("test_data",i)
    tmp<-get(nm)
    #colnames(tmp) <- c("Treated","V1","V2","V3","Y","FSTUWT","SCHID")
    predNames<-colnames(tmp)
    #predNames<-colnames(tmp)[c(2,3,4)]
    predNames2<-c("treated",predNames)
    assign(nm,tmp)
    assign(paste0("formula",1), as.formula(paste("treated~", paste(predNames, collapse="+"))))
    assign(paste0("formula",2), as.formula(paste("Y~", paste(predNames2, collapse="+"))))
    assign(paste0("formula",3), as.formula(paste("Y~", paste(predNames2, collapse="+"),paste("+","W_FSTUWT"))))
    assign(paste0("PSformulaUnWt",1), as.formula(paste("treated~",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+")))))
    assign(paste0("PSformulaWt",1), as.formula(paste("treated~",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+")),paste("+","ate.wt"))))
    #assign(paste0("mr_orig",i), lm(formula2,data=get(paste0("test_data",i)),weights = get(paste0("test_data",i))$FSTUWT))
    assign(paste0("mr_orig",i), lm(formula2,data=get(paste0("test_data",i))))
    assign(paste0("modelSumm_orig",i), summary(get(paste0("mr_orig",i))))
    assign(paste0("modelSE_orig", i) , get(paste0("modelSumm_orig",i))[[4]][[nrow(get(paste0("modelSumm_orig",i))$coefficients)+2]])
    assign(paste0("mat",i), cem(treatment = "treated",data=get(paste0("test_data",i)), drop=c("Y","FSTUWT","W_FSTUWT","SCHID"),keep.all = TRUE))
  }
  for(i in 1:Sims){
    assign(paste0("est",i),att(get(paste0("mat",i)), Y~treated, data=get(paste0("test_data",i))))
    assign(paste0("estFull",i),att(get(paste0("mat",i)), formula2, data=get(paste0("test_data",i))))
    assign(paste0("estFullWt",i),att(get(paste0("mat",i)), formula3, data=get(paste0("test_data",i))))
    assign(paste0("estForest",i),att(get(paste0("mat",i)), formula2, model="forest",data=get(paste0("test_data",i))))
    assign(paste0("SummEst",i), get(paste0("est",i))[[1]][5])
    assign(paste0("SummEstSE",i), get(paste0("est",i))[[1]][6])
    assign(paste0("SummFull",i), get(paste0("estFull",i))[[1]][5])
    assign(paste0("SummFullSE",i), get(paste0("estFull",i))[[1]][6])
    assign(paste0("SummFullWt",i), get(paste0("estFullWt",i))[[1]][5])
    assign(paste0("SummFullWtSE",i), get(paste0("estFullWt",i))[[1]][6])
    assign(paste0("SummestForest",i), get(paste0("estForest",i))[[1]][1])
    assign(paste0("SummestForestSE",i), get(paste0("estForest",i))[[1]][2])
  }
    
    
    #assign(paste0("m.out",i), matchit(PSformulaUnWt1,data=get(paste0("test_data",i)),replace=TRUE, method= "nearest", distance = "logit"))
    #assign(paste0("m.out",i), matchit(formula1,data=get(paste0("test_data",i)),replace=TRUE, method= "nearest", distance = "logit"))
    #assign(paste0("m.data1",i), match.data(get(paste0("m.out",i))))
    #assign(paste0("mr",i), lm(formula2, weights=get(paste0("m.data1",i))$weights, data=get(paste0("m.data1",i))))
    #assign(paste0("mr",i), lm(formula2, data=get(paste0("m.data1",i))))
    #assign(paste0("modelSumm",i), summary(get(paste0("mr",i))))
  
  mattwo<<-mat2
  SummEstTwo<<-SummEst2
  SummFullWtTwo<<-SummFullWt2
  SummFullWtSETwo<<-SummFullWtSE2
  SummestForestTwo<<-SummestForest2
  SummestForestSETwo<<-SummestForestSE2
  PSformulaUnWt1<<-PSformulaUnWt1
  estTwo<<-est2
  estFullTwo<<-estFull2
  estFullWtTwo<<-estFullWt2
  estForestTwo<<-estForest2
  formula2<<-formula2
  formula3<<-formula3
  PSformulaWt1<<-PSformulaWt1
  
  for (i in 1:Sims){
    if(i >= 1 & i <= Sims) {
      counter = i
      print(counter)
    }
    
    #assign(paste0("m.out",i), matchit(formula1,data=get(paste0("test_data",i)), method= "nearest"))
    assign(paste0("t.model",i), glm(PSformulaUnWt1, data=get(paste0("test_data",i)),family="binomial",weights =get(paste0("test_data",i))$W_FSTUWT))
    assign(paste0("pscore",i), predict(get(paste0("t.model",i)), data=get(paste0("test_data",i)), type="response"))
    #tmp3$ate.wt <-ifelse(get(paste0("test_data",i))$treated==1, 1/get(paste0("pscore",i)), 1/(1-get(paste0("pscore",i))))
    assign(paste0("ate.wt",i), ifelse(get(paste0("test_data",i))$treated==1, 1/get(paste0("pscore",i)), 1/(1-get(paste0("pscore",i)))))
    assign(paste0("m.out",i), matchit(PSformulaUnWt1,data=get(paste0("test_data",i)), method= "nearest"))
    assign(paste0("z.out",i), zelig(formula2,data=match.data(get(paste0("m.out",i))),model="ls",cite = FALSE))
    assign(paste0("x.out",i), setx(get(paste0("z.out",i)),treated=0))
    assign(paste0("x1.out",i), setx(get(paste0("z.out",i)),treated=1))
    assign(paste0("s.out",i), sim(get(paste0("z.out",i)),x=get(paste0("x.out",i)),x1=get(paste0("x1.out",i))))
    assign(paste0("modelSumm",i), coef(get(paste0("s.out",i)))[2])
    assign(paste0("modelSE", i) , get_se(get(paste0("z.out",i)))[[1]][[2]])
    #assign(paste0("modelSummZ",i), coef(get(paste0("s.outz",i)))[2])
    #assign(paste0("modelSEZ", i) , get_se(get(paste0("z.outz",i)))[[1]][[2]])
    ###
    assign("summMatched1", summary(get(paste0("m.out",i)))$reduction$'Mean Diff.')
    assign("summMatched2", summary(get(paste0("m.out",i)))$reduction$'eQQ Med')
    assign("summMatched3", summary(get(paste0("m.out",i)))$reduction$'eQQ Mean')
    assign("summMatched4", summary(get(paste0("m.out",i)))$reduction$'eQQ Max')
  }
  for (i in 1:Sims){
    assign(paste0("td",i), cbind(get(paste0("test_data",i)),get(paste0("ate.wt",i))))
    tp<-paste0("td",i)
    tmp3<-get(tp)
    colnames(tmp3)[length(get(paste0("td",i)))]<-"ate.wt"
    assign(tp,tmp3)
    assign(paste0("m.outz",i), matchit(PSformulaUnWt1,data=get(paste0("td",i)), distance=get(paste0("td",i))$ate.wt))
    assign(paste0("z.outz",i), zelig(formula2,data=match.data(get(paste0("m.outz",i))),model="ls",cite = FALSE))
    assign(paste0("x.outz",i), setx(get(paste0("z.outz",i)),treated=0))
    assign(paste0("x1.outz",i), setx(get(paste0("z.outz",i)),treated=1))
    assign(paste0("s.outz",i), sim(get(paste0("z.outz",i)),x=get(paste0("x.outz",i)),x1=get(paste0("x1.outz",i))))
    assign(paste0("modelSummZ",i), coef(get(paste0("s.outz",i)))[2])
    assign(paste0("modelSEZ", i) , get_se(get(paste0("z.outz",i)))[[1]][[2]])
    assign(paste0("m.outcovadj",i), matchit(PSformulaWt1,data=get(paste0("td",i)), method= "nearest"))
    assign(paste0("z.outcovadj",i), zelig(formula3,data=match.data(get(paste0("m.outcovadj",i))),model="ls",cite = FALSE))
    assign(paste0("x.outcovadj",i), setx(get(paste0("z.outcovadj",i)),treated=0))
    assign(paste0("x1.outcovadj",i), setx(get(paste0("z.outcovadj",i)),treated=1))
    assign(paste0("s.outcovadj",i), sim(get(paste0("z.outcovadj",i)),x=get(paste0("x.outcovadj",i)),x1=get(paste0("x1.outcovadj",i))))
    assign(paste0("modelSummCovAdj",i), coef(get(paste0("s.outcovadj",i)))[2])
    assign(paste0("modelSECovAdj", i) , get_se(get(paste0("z.outcovadj",i)))[[1]][[2]])
    assign(paste0("m.out.wt.covadj",i), matchit(PSformulaWt1,data=get(paste0("td",i)), distance=get(paste0("td",i))$ate.wt))
    assign(paste0("z.out.wt.covadj",i), zelig(formula3,data=match.data(get(paste0("m.out.wt.covadj",i))),model="ls",cite = FALSE))
    assign(paste0("x.out.wt.covadj",i), setx(get(paste0("z.out.wt.covadj",i)),treated=0))
    assign(paste0("x1.out.wt.covadj",i), setx(get(paste0("z.out.wt.covadj",i)),treated=1))
    assign(paste0("s.out.wt.covadj",i), sim(get(paste0("z.out.wt.covadj",i)),x=get(paste0("x.out.wt.covadj",i)),x1=get(paste0("x1.out.wt.covadj",i))))
    assign(paste0("modelSummWtCovAdj",i), coef(get(paste0("s.out.wt.covadj",i)))[2])
    assign(paste0("modelSEWtCovAdj", i) , get_se(get(paste0("z.out.wt.covadj",i)))[[1]][[2]])
    
  }
  t.modeltwo<<-t.model2
  m.outtwo<<-m.out2
  z.outtwo<<-z.out2
  test_datatwo <<- test_data2
  ate.wttw<<-ate.wt2
  tdtwo<<-td2
  m.outztwo<<-m.outz2
  m.outcovadj<<-m.outcovadj2
  m.out.wt.covadj<<-m.out.wt.covadj2
  
  
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
  
  TFrame<-as.data.frame(rep(0,Sims),ncol=18)
  for (i in 1:Sims){
    TFrame[i,1]<-get(paste0("modelSumm",i))#$coefficients[[2]]
    TFrame[i,2]<-get(paste0("modelSE",i))#SE for matched
    TFrame[i,3]<-get(paste0("modelSumm_orig",i))$coefficients[[2]]
    TFrame[i,4]<-get(paste0("modelSE_orig",i))
    TFrame[i,5]<-get(paste0("modelSummZ",i))
    TFrame[i,6]<-get(paste0("modelSEZ",i))
    TFrame[i,7]<-get(paste0("modelSummCovAdj",i))
    TFrame[i,8]<-get(paste0("modelSECovAdj",i))
    TFrame[i,9]<-get(paste0("modelSummWtCovAdj",i))
    TFrame[i,10]<-get(paste0("modelSEWtCovAdj",i))
    TFrame[i,11]<-get(paste0("SummEst",i))
    TFrame[i,12]<-get(paste0("SummEstSE",i))
    TFrame[i,13]<-get(paste0("SummFull",i))
    TFrame[i,14]<-get(paste0("SummFullSE",i))
    TFrame[i,15]<-get(paste0("SummFullWt",i))
    TFrame[i,16]<-get(paste0("SummFullWtSE",i))
    TFrame[i,17]<-get(paste0("SummestForest",i))
    TFrame[i,18]<-get(paste0("SummestForestSE",i))
  }
  names(TFrame)<-c("matched","matchedSE", "unmatched","unmatchedSE","matchedWt", "matchedWtSE",
                   "matchedCovAdj","matchedCovAdjSE","matchedWtCovAdj","matchedWtCovAdjSE","SummEst",
                   "SummEstSE","SummFull","SummFullSE","SummFullWt","SummFullWtSE","SummestForest",
                   "SummestForestSE")
  View(TFrame)
  for(i in 1:nrow(TFrame)){
    TFrame$m.bias[i]=abs(2-TFrame$matched[i])
    TFrame$um.bias[i]=abs(2-TFrame$unmatched[i])
    TFrame$m.wt.bias[i]=abs(2-TFrame$matchedWt[i])
    TFrame$m.covadj.bias[i]=abs(2-TFrame$matchedCovAdj[i])
    TFrame$m.wt.covadj.bias[i]=abs(2-TFrame$matchedWtCovAdj[i])
    TFrame$SummEst.bias[i]=abs(2-TFrame$SummEst[i])
    TFrame$SummFull.bias[i]=abs(2-TFrame$SummFull[i])
    TFrame$SummFullWt.bias[i]=abs(2-TFrame$SummFullWt[i])
    TFrame$SummestForest.bias[i]=abs(2-TFrame$SummestForest[i])
    
  }
  View(TFrame)
  write.table(TFrame, file="TFrame2.csv", col.names = TRUE, row.names = FALSE, sep=",",append = FALSE)
  range.um.bias<-range(TFrame$um.bias)
  range.m.bias<-range(TFrame$m.bias)
  #range.m.wt.bias<-range(TFrame$m.wt.bias)
  actual.m<-as.vector(TFrame$matched)
  expected.m<-as.vector(rep(2,Sims))
  RootMean.m<-RMSE(actual.m,expected.m)
  actual.um<-as.vector(TFrame$unmatched)
  expected.um<-as.vector(rep(2,Sims))
  RootMean.um<-RMSE(actual.um,expected.um)
  actual.m.wt<-as.vector(TFrame$matchedWt)
  expected.m.wt<-as.vector(rep(2,Sims))
  RootMean.m.wt<-RMSE(actual.m.wt,expected.m.wt)
  actual.m.covadj<-as.vector(TFrame$matchedCovAdj)
  expected.m.covadj<-as.vector(rep(2,Sims))
  RootMean.m.covadj<-RMSE(actual.m.covadj,expected.m.covadj)
  actual.m.wt.covadj<-as.vector(TFrame$matchedWtCovAdj)
  expected.m.wt.covadj<-as.vector(rep(2,Sims))
  RootMean.m.wt.covadj<-RMSE(actual.m.wt.covadj,expected.m.wt.covadj)
  actual.SummEst<-as.vector(TFrame$SummEst)
  expected.SummEst<-as.vector(rep(2,Sims))
  RootMean.SummEst<-RMSE(actual.SummEst,expected.SummEst)
  actual.SummFull<-as.vector(TFrame$SummFull)
  expected.SummFull<-as.vector(rep(2,Sims))
  RootMean.SummFull<-RMSE(actual.SummFull,expected.SummFull)
  actual.SummFullWt<-as.vector(TFrame$SummFullWt)
  expected.SummFullWt<-as.vector(rep(2,Sims))
  RootMean.SummFullWt<-RMSE(actual.SummFullWt,expected.SummFullWt)
  actual.SummestForest<-as.vector(TFrame$SummestForest)
  expected.SummestForest<-as.vector(rep(2,Sims))
  RootMean.SummestForest<-RMSE(actual.SummestForest,expected.SummestForest)
  
  W<-sprintf("The mean bias of unmatched units is %f with a mean SE of %f, matched units is %f 
             with a mean SE of %f, Weighted matched units is %f with a mean SE of %f, 
             cov adjusted by weight %f with mean SE %f, weighted with cov adjusted by weight
             bias is %f with a mean SE of %f, CEM simple is %f with mean SE of %f, CEM with 
             all predictors the bias is %f with a mean SE of %f, CEM with all predictors plus
             the weight as a covariate is %f with a mean SE of %f, CEM estimated with a random
             forest algorithm is %f with a mean SE of %f", mean(TFrame$um.bias),
             mean(TFrame$unmatchedSE), mean(TFrame$m.bias), mean(TFrame$matchedSE),
             mean(TFrame$m.wt.bias), mean(TFrame$matchedWtSE),mean(TFrame$m.covadj.bias),
             mean(TFrame$matchedCovAdjSE),mean(TFrame$m.wt.covadj.bias),
             mean(TFrame$matchedWtCovAdjSE),mean(TFrame$SummEst.bias),
             mean(TFrame$SummEstSE),mean(TFrame$SummFull.bias),
             mean(TFrame$SummFullSE),mean(TFrame$SummFullWt.bias),
             mean(TFrame$SummFullWtSE),mean(TFrame$SummestForest.bias),
             mean(TFrame$SummestForestSE))
  W1<-sprintf("The range of the unmatched bias estimates is from %f to %f, matched is from %f to %f", range.um.bias[[1]], range.um.bias[[2]], range.m.bias[[1]], range.m.bias[[2]])
  W2<-sprintf("The RMSE of the matched units is %f and unmatched units is %f",RootMean.m, RootMean.um)
  W3<-sprintf("The RMSE of the matched weighted units %f, for covariate adjusted by weights %f,
              and for weighted covariate adjusted by weights %f",RootMean.m.wt,RootMean.m.covadj,
              RootMean.m.wt.covadj)
  W4<-sprintf("The RMSE of Simple CEM is  %f, for Full CEM %f, for Full Weighted CEM %f,
              and for Forest CEM %f",RootMean.SummEst,RootMean.SummFull, RootMean.SummFullWt,
              RootMean.SummestForest)
  
  fileConn1<-file("BiasSumm1.txt")
  fileConn2<-file("BiasSumm2.txt")
  fileConn3<-file("BiasSumm3.txt")
  fileConn4<-file("BiasSumm4.txt")
  fileConn5<-file("BiasBasic.txt")
  writeLines(W1, fileConn1)
  writeLines(W2, fileConn2)
  writeLines(W3, fileConn3)
  writeLines(W4, fileConn4)
  writeLines(W,  fileConn5)
  
  print(W)
  print(W1)
  print(W2)
  print(W3)
  print(W4)
  print(TFrame)
  print(PSformulaUnWt1)
  print(formula2)
  PercMeanDiffTable<<-PercMeanDiffTable
  TFrame<<-TFrame
}

