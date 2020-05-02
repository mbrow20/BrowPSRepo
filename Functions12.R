#This is an extension of Function11 that estimates variances from propensity-score
#adjusted model. The difference in this code is that it applies the
#finalized student weights, as recommended by Pfefferman et al. (1998)
###This code used the Zou et al. (2016) function for variance estimation, yet
#instead of using the simple linear model used in the original function,
#a propensity score model estimated by means of the Imbens and Rubin (2016)
#algorithm is used. No weights are used.Needs to be run in junction with
#ExampleCodeSim15Vars2LevNormWt with Functions5SR.R#############
#R code for variance estimates for covariate adjustment by propensity score analysis
#use simulated data for covariate adjusted variance of the treatment effect. This code
#is modified from Zou et al. (2016). For one, normalized final student weights
#should be included to account for multi-stage sampling (Pfefferman et al., 1998).
#Also, rather than estimating the propensity score from a simple linear model
#as orginally programmed in the code, the PS model from the Imbens and Rubin (2015)
# is used. Estimates are compared.The formula 'PSformulaUnWt1' needs to be generated
#as a global variable by running Functions5SR.R before running this bootstrap 
#procedure.
setwd("C:/Users/mbrow20/Desktop/BrowDiss_Fall2019")
df<-read.csv("test_data_NEW.csv",header=TRUE, sep=",")
df_new<-df[,-c(1,18:19)]
dmat<-as.data.frame(df_new[,c(16,3, 1:2, 4:15)])

# dmat: data matrix
#      1) first column is outcome
#      2) second column is treatment assignment
#      3) the rest are other observed covariates
# boot.num:   # of bootstrapping
#           <=0 indicating two-stage
# output:
#      1) trt.est: treatment effect estimate
#      2) trt.se: standard error based on the proposed estimator
#      3) naive.se: standard error based on conventional variance estimator
twoStagePS = function(dmat, boot.num=100) {
  n=dim(dmat)[1];
  y=dmat[ ,1];
  #trt=dmat[ , 2];
  treat=dmat[,2];
  #x=as.matrix(dmat[,-c(1 ,2)]);#Change--matrix does not recognize columns
  x=dmat[,-1]
  x=as.matrix(x)
  #psfit=glm(trt~x,family=binomial);#compare the ps model with linear covariates and one using the PScovariate selection model. See how different is the variance estimation
  psfit=glm(PSformulaUnWt1,data=dmat, family=binomial)
  ps=psfit$fitted;
  #yfit=lm(y~trt+ps);##will have to estimate the regression with normalized weights as recommended by Pfefferman
  yfit=lm(y~treat+ps);
  ysum=(summary(yfit))$coefficients;
  if (boot.num <= 0) {
    newx=cbind(rep(1,n), x);
    #tps=trt-ps;#R does not recognize
    #trtpscomb<-as.data.frame(cbind(trt,ps))
    trtpscomb<-as.data.frame(cbind(treat,ps))
    trtpscomb$dif<-trtpscomb$treat-trtpscomb$ps
    #trtpscomb$dif<-trtpscomb$trt-trtpscomb$ps
    tps<-trtpscomb$dif
    #11th1=newx*matrix(tps,1,n)[,];#R deos not recognize objects that begin with numbers
    Eleventh1=newx*matrix(tps,1,n)[,];
    y.res=yfit$resid;
    #sigma2=sum(y.res^2)/(n-3);#R does not recognize input
    y.resSquared<-y.res^2
    Sigma2<-sum(y.resSquared)
    Sigma2<-Sigma2/(n-3)
    #tem2=y.res^2/(2*sigma2) - 0.5;#Change according to above modification
    tem2=y.resSquared/(2*Sigma2)-0.5;
    #12th2=cbind(y.res,trt*y.res,ps*y.res,tem2)/sigma2;#See comment on line 34
    #Twelfth2=cbind(y.res,trt*y.res,ps*y.res,tem2)/Sigma2;
    Twelfth2=cbind(y.res,treat*y.res,ps*y.res,tem2)/Sigma2;
    alphaP=ysum[3,1];
    psTab<-as.data.frame(ps,ncol=1)
    psTab$oneMinus<-1-psTab$ps
    psTab$res<-y.res
    psTab$alpha<-alphaP
    psTab$pstem<-psTab$ps*psTab$oneMinus*psTab$res*psTab$alpha/Sigma2
    ppstem1=psTab$pstem
    #ppstem1=ps*(1-ps)*y.res*alphaP/Sigma2;
    #ppstem2=trt*(1-ps)-(1-trt)*ps;#Again, does not recognize this input
    #trtTab<-as.data.frame(trt,ncol=1)
    trtTab<-as.data.frame(treat,ncol=1)
    trtTab$oneMinus<-psTab$oneMinus
    #trtTab$oneMinusTr<-1-trtTab$trt
    trtTab$oneMinusTr<-1-trtTab$treat
    trtTab$ps<-psTab$ps
    #trtTab$pstem2<-trtTab$trt*trtTab$oneMinus-trtTab$oneMinusTr*trtTab$ps
    trtTab$pstem2<-trtTab$treat*trtTab$oneMinus-trtTab$oneMinusTr*trtTab$ps
    ppstem2=trtTab$pstem2
    ppstem=ppstem1+ppstem2;
    #12th1=newx*matrix(ppstem,1,n)[,];#See comments above
    Twelfth1=newx*matrix(ppstem,1,n)[,];
    #V1mat=t(11th1)%*%11th1;
    V1mat=t(Eleventh1)%*%Eleventh1;
    #V2mat=t(12th2)%*%12th2;
    V2mat=t(Twelfth2)%*%Twelfth2;
    #Cmat=t(12th2)%*%12th1;
    Cmat=t(Twelfth2)%*%Twelfth1;
    #Rmat=t(12th2)%*%11th1;
    Rmat=t(Twelfth2)%*%Eleventh1
    new.V1=solve(V1mat);
    new.V2=solve(V2mat);
    CVmat1=as.matrix(Cmat%*%new.V1%*%t(Cmat));
    CVmat2=as.matrix(Rmat%*%new.V1%*%t(Cmat));
    CVmat3=as.matrix(Cmat%*%new.V1%*%t(Rmat));
    #CVmat<- CVmat1 - CVmat2 - CVmat3;#R does not recognize formatting. Also saved CVmats to matrices
    CVmat<-CVmat1 - CVmat2 -CVmat3;
    
    Vmat=new.V2+new.V2%*%CVmat%*%new.V2;
    t.est=ysum[2,1];
    t.se=sqrt(Vmat[2,2]);
    n.se=ysum[2,2];
    return(list(trt.est=t.est,trt.se=t.se,naive.se=n.se));
  }
  #sprintf("The treatment effect is %f and variance estimate for covariate adjustment by the
          #propensity score is %f. The naive s.e. is %f",t.est,t.se,n.se)
  boot.est=rep(0,boot.num);
  index=c(1:n);
  for (i in 1:boot.num) {
    nindex=sample(index,n,replace=T);
    ny=y[nindex];
    #ntrt=trt[nindex];
    ntrt=treat[nindex];
    nx=x[nindex,];
    npsfit=glm(ntrt~nx,family=binomial);
    nps=npsfit$fitted;
    nyfit=lm(ny~ntrt+nps);
    nysum=(summary(nyfit))$coefficients;
    boot.est[i]=nysum[2,1]
    boot.est1<<-boot.est
   }
  t.est=ysum[2,1];
  t.se=sqrt(var(boot.est));
  n.se=ysum[2,2];
  print(t.est)
  print(t.se)
  print(n.se)
  return(list(trt.est=t.est,trt.se=t.se,naive.se=n.se));
}

twoStagePS(dmat = dmat, boot.num = 500)

##The following lines will save the trt estimates and trt.se
##into an Excel file. We can write formulas in Excel
##to test coverage
boot.est2<-as.data.frame(boot.est1)
boot.est2$var<-sqrt(var(boot.est2$boot.est1))
colnames(boot.est2)=c("trt","trt.se")
write.table(boot.est2, file="boot.est2.csv",col.names = TRUE, row.names = FALSE, sep=",", append=FALSE)
