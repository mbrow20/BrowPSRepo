##this is the script to run sensitivity analysis on simulated data
data8<-read.csv("CAN2015_15Vars2LevWithNormWt.csv",header=TRUE,sep=",")
names(data8)[1]="SC048Q01"
data8_control<-subset(data8, data8$treat=="0")
data8_treat<-subset(data8, data8$treat=="1")

##the intercepts mark the the value in the outcome variable (in this case, economic-social-cultural status) after controlling
##for all the remaining variables. We now perform an independent t-test to test the null hypothesis that the  if the 
##expected value for ESCS is similar for the treatment and control group, when controlling for the other
##covariates. If we fail to reject the null, we may assume there is no confounding.

formula17<-as.formula(ESCS~ST011Q07 + HOMEPOS + ST121Q01 + ST034Q01 + SC048Q01 + 
  Female + ST118Q04 + RATCMP2 + ST013Q01 + HISCED + MISCED+I(SC048Q01^2) + 
  SC048Q01:RATCMP2 + ST034Q01:ST118Q04 + ST121Q01:ST013Q01 + HISCED:MISCED + Female:ST013Q01 + 
  SC048Q01:ST013Q01 + I(ST013Q01^2) + HOMEPOS:RATCMP2 + I(RATCMP2^2) + ST034Q01:RATCMP2 + 
  ST118Q04:RATCMP2 + ST011Q07:ST118Q04 + I(ST034Q01^2) + Female:ST118Q04 + 
  Female:HISCED + ST121Q01:SC048Q01 + ST011Q07:RATCMP2)
lm17_control<-lm(formula17, data=data8_control)
lm17_treat <-lm(formula17, data = data8_treat)
summary(lm17_control)

##determine pooled standard error to conduct t-test to determine significance of
##difference between intercepts in the treated and control models

sp_square = ((nrow(data8_control)-1)*(summary(lm17_control)$coefficients[[31]])^2+
               (nrow(data8_treat)-1)*(summary(lm17_treat)$coefficients[[31]])^2)/(nrow(data8_control)+nrow(data8_treat)-2)
sp<-sqrt(sp_square)

##conduct independent t-test...difference in intercepts divided by the pooled standard error

independ_ttest<-(summary(lm17_control)$coefficients[[1]]-summary(lm17_treat)$coefficients[[1]])/sp
independ_ttest
