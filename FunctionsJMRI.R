##This function will perform multiple imputation using Joint Model Random Intercepts (van Buuren, 2012)##
JointModelRandomInterceptsPS=function(data,varNames){
  indxNum<-dim(data)[2]
  indxNum2<-indxNum-2
  varNames<-colnames(data[,c(2:indxNum2)])
  variablesForJMRI<-paste(varNames,collapse="+")
  variablesForJMRI<<-paste(varNames,collapse="+")
  variablesForJMRI2<-paste(variablesForJMRI,"~1+(1|SCHID)")
  variablesForJMRI2<<-paste(variablesForJMRI,"~1+(1|SCHID)")
  variablesForJMRI3<-paste("Y+",variablesForJMRI2)
  variablesForJMRI3<<-paste("Y+",variablesForJMRI2)
  fml<-variablesForJMRI3
  fml<<-variablesForJMRI3
  imp_pan<-panImpute(data,formula=as.formula(fml),n.burn=500,n.iter=100,m=5,seed=1234)
  imp_pan<<-panImpute(data,formula=as.formula(fml),n.burn=500,n.iter=100,m=5,seed=1234)
  impList<-mitmlComplete(imp_pan,print="all")
  impList<<-mitmlComplete(imp_pan,print="all")
  write.table(impList[[1]], file="USA2015REL_PSIMPUTE1.csv",col.names = T, row.names = F, sep=",", append=F)
  data4<- impList[[1]]
  data4<<-impList[[1]]
  #dpart<- createDataPartition(data4$treat,p=0.2,list=F)
  #dpart<<-createDataPartition(data4$treat,p=0.2,list=F)
  #data4<- data4[dpart,]
  #data4<<-data4[dpart,]
  indxer=  which( colnames(data4)=="treat" )
  indxer<<-which( colnames(data4)=="treat" )
  data4NotIndxer<- data4[,-(indxer)]
  data4NotIndxer<<-data4[,-(indxer)]
  data4Indxer<- data4[,indxer]
  data4Indxer<<-data4[,indxer]
  data4<-cbind(data4Indxer,data4NotIndxer)
  data4<<-cbind(data4Indxer,data4NotIndxer)
  indxY=which(colnames(data4)=="Y")
  indxY<<-which(colnames(data4)=="Y")
  data4<- data4[,-indxY]
  data4<<-data4[,-indxY]
  data4<- as.data.frame(data4)
  data4<<-as.data.frame(data4)
  colnames(data4)[1]="treat"
  colnames(data4)[1]<<-"treat"
} 
JointModelRandomInterceptsPS(data1,varNames)
  indxSCHID=which(colnames(data4)=="SCHID")
  data5<<-data4[,-indxSCHID]
  indxCNTSTUID=which(colnames(data5)=="CNTSTUID")
  data5<<-data5[,-indxCNTSTUID]
  PS.CovSelection.Function(data5)
  PSformulaUnWt<<-as.formula(paste("treat~",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+"))))
  formulaRegUnWt<<-as.formula(paste("Y~treat+",paste(linearVars,collapse="+"),paste("+",paste(interactionVars,collapse="+"))))
  print(PSformulaUnWt)


JointModelRandomInterceptsFULL=function(data){
  imp_pan<-panImpute(data,formula=as.formula(fml),n.burn=500,n.iter=100,m=5,seed=1234)
  imp_pan<<-panImpute(data,formula=as.formula(fml),n.burn=500,n.iter=100,m=5,seed=1234)
  impList<-mitmlComplete(imp_pan,print="all")
  impList<<-mitmlComplete(imp_pan,print="all")
  write.table(impList[[1]], file="USA2015REL_FULLIMPUTE1.csv",col.names = T, row.names = F, sep=",", append=F)
  write.table(impList[[2]], file="USA2015REL_FULLIMPUTE2.csv",col.names = T, row.names = F, sep=",", append=F)
  write.table(impList[[3]], file="USA2015REL_FULLIMPUTE3.csv",col.names = T, row.names = F, sep=",", append=F)
  write.table(impList[[4]], file="USA2015REL_FULLIMPUTE4.csv",col.names = T, row.names = F, sep=",", append=F)
  write.table(impList[[5]], file="USA2015REL_FULLIMPUTE5.csv",col.names = T, row.names = F, sep=",", append=F)
}
dfFULL1<- read.csv ("USA2015REL_FULLIMPUTE1.csv",header=TRUE, sep=",")
dfFULL1<<-read.csv ("USA2015REL_FULLIMPUTE1.csv",header=TRUE, sep=",")
dfFULL1<- as.data.frame(dfFULL1)
dfFULL1<<-as.data.frame(dfFULL1)
m.out1<- matchit(PSformulaUnWt, data=dfFULL1, method="nearest", distance="logit")
m.out1<<-matchit(PSformulaUnWt, data=dfFULL1, method="nearest", distance="logit")
m.out1<<-m.out1
dfFULL2<- read.csv ("USA2015REL_FULLIMPUTE2.csv",header=TRUE, sep=",")
dfFULL2<<-read.csv ("USA2015REL_FULLIMPUTE2.csv",header=TRUE, sep=",")
dfFULL2<- as.data.frame(dfFULL2)
dfFULL2<<-as.data.frame(dfFULL2)
m.out2<- matchit(PSformulaUnWt, data=dfFULL2, method="nearest", distance="logit")
m.out2<<-matchit(PSformulaUnWt, data=dfFULL2, method="nearest", distance="logit")
m.out2<<-m.out2
dfFULL3<- read.csv ("USA2015REL_FULLIMPUTE3.csv",header=TRUE, sep=",")
dfFULL3<<-read.csv ("USA2015REL_FULLIMPUTE3.csv",header=TRUE, sep=",")
dfFULL3<- as.data.frame(dfFULL3)
dfFULL3<<-as.data.frame(dfFULL3)
m.out3<- matchit(PSformulaUnWt, data=dfFULL3, method="nearest", distance="logit")
m.out3<<-matchit(PSformulaUnWt, data=dfFULL3, method="nearest", distance="logit")
m.out3<<-m.out3
dfFULL4<- read.csv ("USA2015REL_FULLIMPUTE4.csv",header=TRUE, sep=",")
dfFULL4<<-read.csv ("USA2015REL_FULLIMPUTE4.csv",header=TRUE, sep=",")
dfFULL4<- as.data.frame(dfFULL4)
dfFULL4<<-as.data.frame(dfFULL4)
m.out4<- matchit(PSformulaUnWt, data=dfFULL4, method="nearest", distance="logit")
m.out4<<-matchit(PSformulaUnWt, data=dfFULL4, method="nearest", distance="logit")
m.out4<<-m.out4
dfFULL5<- read.csv ("USA2015REL_FULLIMPUTE5.csv",header=TRUE, sep=",")
dfFULL5<<-read.csv ("USA2015REL_FULLIMPUTE5.csv",header=TRUE, sep=",")
dfFULL5<- as.data.frame(dfFULL5)
dfFULL5<<-as.data.frame(dfFULL5)
m.out5<- matchit(PSformulaUnWt, data=dfFULL5, method="nearest", distance="logit")
m.out5<<-matchit(PSformulaUnWt, data=dfFULL5, method="nearest", distance="logit")
m.out5<<-m.out5
vec_One<- as.vector(row.names(m.out$match.matrix))
vec_One<<-as.vector(row.names(m.out$match.matrix))
vec_One<- as.numeric(vec_One)
vec_One<<-as.numeric(vec_One)
vec_Two<- as.vector(m.out$match.matrix[,1])
vec_Two<<-as.vector(m.out$match.matrix[,1])
vec_Two<- as.numeric(vec_Two)
vec_Two<<-as.numeric(vec_Two)
vec_comb<- c(vec_One, vec_Two)
vec_comb<<-c(vec_One, vec_Two)

USA2015REL_FULLIMPUTE1<-read.csv("USA2015REL_FULLIMPUTE1.csv", header=T, sep=",")
USA2015REL_FULLIMPUTE2<-read.csv("USA2015REL_FULLIMPUTE2.csv", header=T, sep=",")
USA2015REL_FULLIMPUTE3<-read.csv("USA2015REL_FULLIMPUTE3.csv", header=T, sep=",")
USA2015REL_FULLIMPUTE4<-read.csv("USA2015REL_FULLIMPUTE4.csv", header=T, sep=",")
USA2015REL_FULLIMPUTE5<-read.csv("USA2015REL_FULLIMPUTE5.csv", header=T, sep=",")

write.table(USA2015REL_FULLIMPUTE1[vec_comb,], file="USA2015REL_MATCHEDIMPUTE1.csv", col.names = T, row.names = F, sep = ",", append=F)
write.table(USA2015REL_FULLIMPUTE2[vec_comb,], file="USA2015REL_MATCHEDIMPUTE2.csv", col.names = T, row.names = F, sep = ",", append=F)
write.table(USA2015REL_FULLIMPUTE3[vec_comb,], file="USA2015REL_MATCHEDIMPUTE3.csv", col.names = T, row.names = F, sep = ",", append=F)
write.table(USA2015REL_FULLIMPUTE4[vec_comb,], file="USA2015REL_MATCHEDIMPUTE4.csv", col.names = T, row.names = F, sep = ",", append=F)
write.table(USA2015REL_FULLIMPUTE5[vec_comb,], file="USA2015REL_MATCHEDIMPUTE5.csv", col.names = T, row.names = F, sep = ",", append=F)
#options(warn=0)