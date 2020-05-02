JointModelRandomInterceptsFULL=function(data){
  imp_pan<-panImpute(data,formula=as.formula(fml),n.burn=500,n.iter=100,m=5,seed=1234)
  imp_pan<<-panImpute(data,formula=as.formula(fml),n.burn=500,n.iter=100,m=5,seed=1234)
  impList<-mitmlComplete(imp_pan,print="all")
  impList<<-mitmlComplete(imp_pan,print="all")
  write.table(impList[[1]], file="CAN2015REL_FULLIMPUTE1.csv",col.names = T, row.names = F, sep=",", append=F)
  write.table(impList[[2]], file="CAN2015REL_FULLIMPUTE2.csv",col.names = T, row.names = F, sep=",", append=F)
  write.table(impList[[3]], file="CAN2015REL_FULLIMPUTE3.csv",col.names = T, row.names = F, sep=",", append=F)
  write.table(impList[[4]], file="CAN2015REL_FULLIMPUTE4.csv",col.names = T, row.names = F, sep=",", append=F)
  write.table(impList[[5]], file="CAN2015REL_FULLIMPUTE5.csv",col.names = T, row.names = F, sep=",", append=F)
}