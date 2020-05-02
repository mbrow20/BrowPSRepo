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
  #options(warn=-1)
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
  
  #PSformulaUnWt<<-as.formula(paste("treat~",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+"))))
  #formulaRegUnWt<<-as.formula(paste("Y~treat+",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+"))))
  #print(PSformulaUnWt)
}