##This is the Imbens and Rubing (2015) variance estimator functions. It needs to be run
##right after Functions5VAR500.R. The Functions5VAR500.R will create matched sets in the
##global environment that the ImbensRubin function will later utilize


ImbensRubin=function(data, Sims){
  m.treat<-as.vector(row.names(data$match.matrix))
  m.treated_df<-test_data[m.treat,]
  m.treatDist<-data$distance[m.treat]
  MTreatWithDist<-as.data.frame(cbind(m.treated_df,m.treatDist))
  MTreatSorted<-MTreatWithDist[order(m.treatDist),]
  is.even <- function(x) x %% 2 == 0 #Even though matched data by ratio = 1 are matched, other ratios may not be
  N_match<-ifelse(is.even(nrow(MTreatSorted=="TRUE")),nrow(MTreatSorted)-1,nrow(MTreatSorted)-2)
  se_df<-NULL
  listOne<-as.vector(seq(1,N_match,2))
  se_df=MTreatSorted$Y[listOne]
  se_df2<-NULL
  listTwo<-as.vector(seq(2,N_match+1,2))
  se_df2<-MTreatSorted$Y[listTwo]
  se_tot<-as.data.frame(cbind(se_df,se_df2))
  se_tot$diff<-se_tot$se_df-se_tot$se_df2
  se_tot$sq<-se_tot$diff^2
  se_tot$divtwo<-se_tot$sq/2
  m.var=mean(se_tot$divtwo)
  sd.mvar<-sqrt(m.var)
  m.se<-sd.mvar/sqrt(nrow(MTreatSorted))
  print(m.se)
  ##Extracting the numbers of matching control units from test_data and matching distances
  m.control<-as.vector(data$match.matrix[,1])
  m.control_df<-test_data[m.control,]
  m.controlDist<-data$distance[m.control]
  MControlWithDist<-as.data.frame(cbind(m.control_df,m.controlDist))
  MControlSorted<-MControlWithDist[order(m.controlDist),]
  N_control<-ifelse(is.even(nrow(MControlSorted=="TRUE")),nrow(MControlSorted)-1,nrow(MControlSorted)-2)
  se_dfC<-NULL
  listOneC<-as.vector(seq(1,N_control,2))
  se_dfC=MControlSorted$Y[listOneC]
  se_df2C<-NULL
  listTwoC<-as.vector(seq(2,N_control+1,2))
  se_df2C<-MControlSorted$Y[listTwoC]
  se_totC<-as.data.frame(cbind(se_dfC,se_df2C))
  se_totC$diff<-se_totC$se_dfC-se_totC$se_df2C
  se_totC$sq<-se_totC$diff^2
  se_totC$divtwo<-se_totC$sq/2
  m.varC=mean(se_totC$divtwo)
  sd.mvarC<-sqrt(m.varC)
  m.seC<-sd.mvarC/sqrt(nrow(MControlSorted))
  print(m.seC)
  comb_se<<-(m.se+m.seC)/2
  sprintf("The Imbens and Rubin(2015) variance estimation for a matched sample is %f", comb_se)
}
ave_se<-NULL
Sims=500
for (i in 1:Sims){
  ImbensRubin(get(paste0("m.outx",i)),Sims)
  ave_se[i]=comb_se
}

#ave_se<<-ave_se
sprintf("The Imbens and Rubin(2015) variance estimation for a matched sample is %f over %d simulations", comb_se, Sims)