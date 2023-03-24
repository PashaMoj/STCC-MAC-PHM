rm(list=ls())
##Economic Design (ED) and Economic Statistical Design (ESD) of
##""Shewhart-Type"" Control Charts Under "Multiple"" Assignable Causes
##Under Different Process Failure Mechanisms (PFM) as proportional Hazards Model (PHM) 
##and Sampling Schemes (SS)
##Speciically, for
##Xbar (with Normal or Non-normal Distributions) and T2 Control Charts
##(The Lines alpha and beta Could be Modified for any Other Control Statistics.)
##Under Weibull and Chen PFMs
##(The Corresponding Lines could be modified for Other PHM-PFMs.)
##Under Non-uniform Sampling Schemes with Constant Integrated Hazard
##(The Corresponding Lines Could be modified for any Other SSs.)

##The Main Output Quantities of This Code is the Optimal Design Parameters
##(Sample Size n, Sampling intervals based on h1, Control Limit Coefficient L)
#n=x1;  h1=x2;  L=x3
###############

## Control Statistics Distribution

INorm<-readline(prompt="Enter 1 for the Univariate Normal Case, Otherwise Enter 0:")
INorm<-as.numeric(INorm)
if (INorm!=1 & INorm!=0) {print("Error: Enter 1 or 0")}

INNorm<-readline(prompt="Enter 1 or the Univariate Non-normal(Burr) Case, Otherwise Enter 0:")
INNorm<-as.numeric(INNorm)
if (INNorm!=1 & INNorm!=0){print("Error: Enter 1 or 0")}

INB=sum(INorm,INNorm)
if (INB!=1 & INB!=0){print("Error: You Should choose one and only one of the above models.")}

IT2<-readline(prompt="Enter 1 for T2 Chart, Otherwise Enter 0:")
IT2<-as.numeric(IT2)
if (IT2!=1 & IT2!=0){print("Error: Enter 1 or 0")}

if (sum(INB,IT2)!=1){print("Error: You Should choose one and only one of the above models. For Other Control Statistics, Modify the Lines alpha and beta.")}
###############

##Shock Model/Process Failure Mechanism (PFM) and its Parameters
##with its Expected Value ETa

IW<-readline(prompt="Enter 1 for the Weibull Shock Model, Otherwise Enter 0:")
IW<-as.numeric(IW)
if (IW!=1 & IW!=0){print("Error: Enter 1 or 0")}

ICH<-readline(prompt="Enter 1 for the Chen Shock Model, Otherwise Enter 0:")
ICH<-as.numeric(ICH)
if (ICH!=1 & ICH!=0){print("Error: Enter 1 or 0")}

if (sum(IW,ICH)!=1){print("Error: You Should choose one and only one of the above models. For other PFMs, Modify the corresponding Lines.")}
###############

##Prior Distributions

INE<-readline(prompt="Enter 1 for the Negative Exponential Prior, Otherwise Enter 0:")
INE<-as.numeric(INE)
if (INE!=1 & INE!=0) {print("Error: Enter 1 or 0")}

IUN<-readline(prompt="Enter 1 or the Uniform Prior, Otherwise Enter 0:")
IUN<-as.numeric(IUN)
if (IUN!=1 & IUN!=0){print("Error: Enter 1 or 0")}

IHN<-readline(prompt="Enter 1 for Half Normal Prior, Otherwise Enter 0:")
IHN<-as.numeric(IHN)
if (IHN!=1 & IHN!=0){print("Error: Enter 1 or 0")}

if (sum(INE,IUN,IHN)!=1){print("Error: You Should choose one and only one of the above models. For Other Prior Distributions, Modify the Corresponding Lines.")}
###############

##Model Cost and Time Parameters for comparison study
##From Yang and Rahim (2005) and Pasha et. al (2023) for T2 Control Charts (IT2=1 & INB=0)
#and From Chen and Yang (2002) and Pasha et. al (2018) for Xbar Charts (IT2=0 & INB=1)
Y=IT2*500+INB*2000
D0=IT2*50+INB*210
a=IT2*20+INB*20
b=IT2*4.22+INB*20
Z0=IT2*.25+INB*1.25
Z1=IT2*.25+INB*1.25

##Model Input Quantities

##Shift Value/Shift Vector for Out-of-Control Status in Univariate/Multivariate Quality Data
#Univariate
delta<-c(1,1.5,1.8,2,2.2,2.5,3)
##Burr Parameters
cb=5;kb=6;mb=0.65513;sb=0.16103

##Multivariate
DF=2 #Chi-square Degrees of Freedom
sigmaM=matrix(c(2,1,1,2.5),nrow=2) #Known Covariance Matrix

deltaM=matrix(ncol=7,nrow = 2)
deltaM[,1]=c(0,0.25*sqrt(2))
deltaM[,2]=c(0,0.5*sqrt(2))
deltaM[,3]=c(0,0.75*sqrt(2))
deltaM[,4]=c(0,1*sqrt(2))
deltaM[,5]=c(0,1.25*sqrt(2))
deltaM[,6]=c(0,1.5*sqrt(2))
deltaM[,7]=c(0,1.75*sqrt(2))
  
NCP=c()
for (i in 1:7){
NCP[i]=as.numeric(t(deltaM[,i])%*%solve(sigmaM)%*%deltaM[,i]) #Chi-square Non-central Parameter
}
distM<-sqrt(NCP) #The column di in the table of model input parameters

D1i<-IT2*c(56.21,160.21,410.36,950,1973.08,3687.57,6222.78)+INB*c(575,1684,2901,4000,5341,7776,12602)

landa<-IT2*INE*c(.006527,.005757,.005091,.004484,.003966,.003493,.003078)+IT2*IUN*c(.003529,.003529,.003529,.003529,.003529,.003529,.003529)+IT2*IHN*c(.004678,.004572,.004394,.004158,.003875,.003556,.003213)+INB*INE*c(.004566,.003557,.003059,.002772,.002502,.002155,.001689)+INB*IUN*c(.002294,.002294,.002294,.002294,.002294,.002294,.002294)+INB*IHN*c(.004220,.003608,.003190,.002901,.002612,.002194,.001557)
landa0<-sum(landa)

Z2i<-IT2*INE*c(1.092,.963,.851,.75,.663,.584,.515)+IT2*IUN*c(.75,.75,.75,.75,.75,.75,.75)+IT2*IHN*c(.844,.825,.793,.75,.699,.641,.580)+INB*INE*c(3.293,2.565,2.206,2,1.804,1.554,1.217)+INB*IUN*c(2,2,2,2,2,2,2)+INB*IHN*c(2.909,2.488,2.198,2,1.802,1.512,1.074)

wi<-IT2*INE*c(1601,1412.2,1248.8,1100,972.9,856.8,755.1)+IT2*IUN*c(1100,1100,1100,1100,1100,1100,1100)+IT2*IHN*c(1237.5,1209.4,1162.5,1100,1025,940.6,850)+INB*INE*c(1647,1283,1103,1000,902,777,609)+INB*IUN*c(1000,1000,1000,1000,1000,1000,1000)+INB*IHN*c(1454,1244,1099,1000,901,756,537)
###############

##mu is the expected value of the process failure mechanism (ETa)
#PFM (Common) Shape Parameter
nu=2

##for Weibull PFM
mu0W=((1/landa0)^(1/nu))*gamma(1+(1/nu))
muW=c()
for (i in 1:7){
  muW[i]=((1/landa[i])^(1/nu))*gamma(1+(1/nu))}

##for Chen PFM
integral = function(t){
  t*(landa0*nu*t^(nu-1)*exp(landa0*(1-exp(t^nu))+t^nu))}
mu0CH=integrate(integral,lower=0,upper=Inf)$value
muCH=c()
for (i in 1:7){
  integral = function(t){
    t*(landa[i]*nu*t^(nu-1)*exp(landa[i]*(1-exp(t^nu))+t^nu))}
  muCH[i]=integrate(integral,lower=0,upper=Inf)$value}

mu0=IW*(mu0W)+ICH*(mu0CH)
mu=IW*(muW)+ICH*(muCH)
###############

##Convergence Criteria
TOL=1e-6

##Design Parameters for Non-Uniform Scheme with Constant IHI
##n=x1 ; h1=x2 ; L=x3

ECTPHM<-function(x){
  x1<-x[1]
  x2<-x[2]
  x3<-x[3]

  ## Type I and Type II Error Probabilities (T2 Chart)
  alpha<-INorm*(2*pnorm(-x3,0,1))+INNorm*(1+(1/((1+(mb+x3*sb)^cb)^kb))-(1/((1+(mb-x3*sb)^cb)^kb)))+IT2*(1- pchisq(x3, df=DF))
  
  beta=c()
  for (i in 1:7){
    beta[i]=INorm*(pnorm(sqrt(x1)*delta[i]+x3,0,1)-pnorm(sqrt(x1)*delta[i]-x3,0,1))+INNorm*((1/((1+(mb-x3*sb-sb*delta[i]*sqrt(x1))^cb)^kb))-(1/(1+(mb+x3*sb-sb*delta[i]*sqrt(x1))^cb)^kb))+IT2*(pchisq(x3, df=DF, ncp =x1*NCP[i] , lower.tail = TRUE))
    }
  
  ##Constant Shift Probability for each sampling interval
  ##For Weibull PFM
  p0W=(1-exp(-landa0*(x2)^nu))
  pW=c()
  for (i in 1:7){
    pW[i]=(1-exp(-landa[i]*(x2)^nu))}
  ##For Chen PFM
  p0CH=(1-exp(landa0-landa0*exp(x2^nu)))
  pCH=c()
  for (i in 1:7){
    pCH[i]=(1-exp(landa[i]-landa[i]*exp(x2^nu)))}
p0=IW*p0W+ICH*p0CH
p=IW*pW+ICH*pCH

  yW=yCH=c();yW[1]=yCH[1]=x2
  j=1;S1=S2=1000;sum1=0;c1=c2=0
  while(c1<1){
    q=j
    while(c2<1){
      ##For Weibull PFM
      yW[q]=(q^{1/nu}-(q-1)^{1/nu})*x2
      ##For Chen PFM
      yCH[q]=(log(1-q*(1-exp(x2^nu))))^{1/nu}-(log(1-(q-1)*(1-exp(x2^nu))))^{1/nu}
      wk=sum(IW*yW[1:q]+ICH*yCH[1:q])
      
      sum1=sum1+((1-p[1])^(j-1))*wk*(beta[1]^(q-j))
      if(abs(sum1-S1) == "NaN" | abs(sum1-S1) == "NA"){break}
      else if (abs(sum1-S1)<TOL){break}
      else{S1=sum1;q=q+1}}
    if (abs(sum1-S2) == "NaN" | abs(sum1-S2) == "NA"){break}
    else if(abs(sum1-S2)<TOL){break}
    else{S2<-sum1;j=j+1}}

  yW=yCH=c();yW[1]=yCH[1]=x2
  j=1;S1=S2=1000;sum2=0;c1=c2=0
  while(c1<1){
    q=j
    while(c2<1){
      ##For Weibull PFM
      yW[q]=(q^{1/nu}-(q-1)^{1/nu})*x2
      ##For Chen PFM
      yCH[q]=(log(1-q*(1-exp(x2^nu))))^{1/nu}-(log(1-(q-1)*(1-exp(x2^nu))))^{1/nu}
      wk=sum(IW*yW[1:q]+ICH*yCH[1:q])
      sum2=sum2+((1-p[2])^(j-1))*wk*(beta[2]^(q-j))
      if (abs(sum2-S1) == "NaN" | abs(sum2-S1) == "NA"){break}
      else if(abs(sum2-S1)<TOL){break}
      else{S1=sum2;q=q+1}}
    if (abs(sum2-S2) == "NaN" | abs(sum2-S2) == "NA"){break}
    else if(abs(sum2-S2)<TOL){break}
    else{S2<-sum2;j=j+1}}

  yW=yCH=c();yW[1]=yCH[1]=x2
  j=1;S1=S2=1000;sum3=0;c1=c2=0
  while(c1<1){
    q=j
    while(c2<1){
      ##For Weibull PFM
      yW[q]=(q^{1/nu}-(q-1)^{1/nu})*x2
      ##For Chen PFM
      yCH[q]=(log(1-q*(1-exp(x2^nu))))^{1/nu}-(log(1-(q-1)*(1-exp(x2^nu))))^{1/nu}
      wk=sum(IW*yW[1:q]+ICH*yCH[1:q])
      sum3=sum3+((1-p[3])^(j-1))*wk*(beta[3]^(q-j))
      if (abs(sum3-S1) == "NaN" | abs(sum3-S1) == "NA"){break}
      else if(abs(sum3-S1)<TOL){break}
      else{S1=sum3;q=q+1}}
    if (abs(sum3-S2) == "NaN" | abs(sum3-S2) == "NA"){break}
    
    else if(abs(sum3-S2)<TOL){break}
    else{S2<-sum3;j=j+1}}

  yW=yCH=c();yW[1]=yCH[1]=x2
  j=1;S1=S2=1000;sum4=0;c1=c2=0
  while(c1<1){
    q=j
    while(c2<1){
      ##For Weibull PFM
      yW[q]=(q^{1/nu}-(q-1)^{1/nu})*x2
      ##For Chen PFM
      yCH[q]=(log(1-q*(1-exp(x2^nu))))^{1/nu}-(log(1-(q-1)*(1-exp(x2^nu))))^{1/nu}
      wk=sum(IW*yW[1:q]+ICH*yCH[1:q])
      sum4=sum4+((1-p[4])^(j-1))*wk*(beta[4]^(q-j))
      if (abs(sum4-S1) == "NaN" | abs(sum4-S1) == "NA"){break}
      else if(abs(sum4-S1)<TOL){break}
      else{S1=sum4;q=q+1}}
    if (abs(sum4-S2) == "NaN" | abs(sum4-S2) == "NA"){break}
    else if(abs(sum4-S2)<TOL){break}
    else{S2<-sum4;j=j+1}}

  yW=yCH=c();yW[1]=yCH[1]=x2
  j=1;S1=S2=1000;sum5=0;c1=c2=0
  while(c1<1){
    q=j
    while(c2<1){
      ##For Weibull PFM
      yW[q]=(q^{1/nu}-(q-1)^{1/nu})*x2
      ##For Chen PFM
      yCH[q]=(log(1-q*(1-exp(x2^nu))))^{1/nu}-(log(1-(q-1)*(1-exp(x2^nu))))^{1/nu}
      wk=sum(IW*yW[1:q]+ICH*yCH[1:q])
      sum5=sum5+((1-p[5])^(j-1))*wk*(beta[5]^(q-j))
      if (abs(sum5-S1) == "NaN" | abs(sum5-S1) == "NA"){break}
      else if(abs(sum5-S1)<TOL){break}
      else{S1=sum5;q=q+1}}
    if (abs(sum5-S2) == "NaN" | abs(sum5-S2) == "NA"){break}
    else if(abs(sum5-S2)<TOL){break}
    else{S2<-sum5;j=j+1}}

  yW=yCH=c();yW[1]=yCH[1]=x2
  j=1;S1=S2=1000;sum6=0;c1=c2=0
  while(c1<1){
    q=j
    while(c2<1){
      ##For Weibull PFM
      yW[q]=(q^{1/nu}-(q-1)^{1/nu})*x2
      ##For Chen PFM
      yCH[q]=(log(1-q*(1-exp(x2^nu))))^{1/nu}-(log(1-(q-1)*(1-exp(x2^nu))))^{1/nu}
      wk=sum(IW*yW[1:q]+ICH*yCH[1:q])
      sum6=sum6+((1-p[6])^(j-1))*wk*(beta[6]^(q-j))
      if (abs(sum6-S1) == "NaN" | abs(sum6-S1) == "NA"){break}
      else if(abs(sum6-S1)<TOL){break}
      else{S1=sum6;q=q+1}}
    if (abs(sum6-S2) == "NaN" | abs(sum6-S2) == "NA"){break}
    else if(abs(sum6-S2)<TOL){break}
    else{S2<-sum6;j=j+1}}

  yW=yCH=c();yW[1]=yCH[1]=x2
  j=1;S1=S2=1000;sum7=0;c1=c2=0
  while(c1<1){
    q=j
    while(c2<1){
      ##For Weibull PFM
      yW[q]=(q^{1/nu}-(q-1)^{1/nu})*x2
      ##For Chen PFM
      yCH[q]=(log(1-q*(1-exp(x2^nu))))^{1/nu}-(log(1-(q-1)*(1-exp(x2^nu))))^{1/nu}
      wk=sum(IW*yW[1:q]+ICH*yCH[1:q])
      sum7=sum7+((1-p[7])^(j-1))*wk*(beta[7]^(q-j))
      if (abs(sum7-S1) == "NaN" | abs(sum7-S1) == "NA"){break}
      else if(abs(sum7-S1)<TOL){break}
      else{S1=sum7;q=q+1}}
    if (abs(sum7-S2) == "NaN" | abs(sum7-S2) == "NA"){break}
    else if(abs(sum7-S2)<TOL){break}
    else{S2<-sum7;j=j+1}}

  sum=c(sum1,sum2,sum3,sum4,sum5,sum6,sum7)

  AVGOOT=c()
  for (i in 1:7){
    AVGOOT[i]=p[i]*(1-beta[i])*sum[i]-mu[i]
    }
  summ1<-(landa[1]/landa0)*(Z1+Z2i[1]+AVGOOT[1])+
  (landa[2]/landa0)*(Z1+Z2i[2]+AVGOOT[2])+
    (landa[3]/landa0)*(Z1+Z2i[3]+AVGOOT[3])+
    (landa[4]/landa0)*(Z1+Z2i[4]+AVGOOT[4])+
    (landa[5]/landa0)*(Z1+Z2i[5]+AVGOOT[5])+
    (landa[6]/landa0)*(Z1+Z2i[6]+AVGOOT[6])+
    (landa[7]/landa0)*(Z1+Z2i[7]+AVGOOT[7])
  summ2<-(landa[1]/landa0)*(D1i[1]*AVGOOT[1]+wi[1])+
    (landa[2]/landa0)*(D1i[2]*AVGOOT[2]+wi[2])+
    (landa[3]/landa0)*(D1i[3]*AVGOOT[3]+wi[3])+
    (landa[4]/landa0)*(D1i[4]*AVGOOT[4]+wi[4])+
    (landa[5]/landa0)*(D1i[5]*AVGOOT[5]+wi[5])+
    (landa[6]/landa0)*(D1i[6]*AVGOOT[6]+wi[6])+
    (landa[7]/landa0)*(D1i[7]*AVGOOT[7]+wi[7])
  summ3<-(landa[1]/landa0)*(1/(1-beta[1]))+(landa[2]/landa0)*(1/(1-beta[2]))+
    (landa[3]/landa0)*(1/(1-beta[3]))+(landa[4]/landa0)*(1/(1-beta[4]))+
    (landa[5]/landa0)*(1/(1-beta[5]))+(landa[6]/landa0)*(1/(1-beta[6]))+
    (landa[7]/landa0)*(1/(1-beta[7]))

  ET<-mu0+summ1+Z0*alpha*(1-p0)/p0
  EC<-D0*mu0+summ2+(a+b*x1)*((1-p0)/p0)+(a+b*x1)*summ3+alpha*Y*((1-p0)/p0)
  EC/ET
}

hin.hs100<-function(x){
  
  alpha<-INorm*(2*pnorm(-x[3],0,1))+INNorm*(1+(1/((1+(mb+x[3]*sb)^cb)^kb))-(1/((1+(mb-x[3]*sb)^cb)^kb)))+IT2*(1- pchisq(x[3], df=DF))
  beta1=INorm*((pnorm(sqrt(x[1])*delta[1]+x[3],0,1)-pnorm(sqrt(x[1])*delta[1]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[1]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[1]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[1] , lower.tail = TRUE))
  beta2=INorm*((pnorm(sqrt(x[1])*delta[2]+x[3],0,1)-pnorm(sqrt(x[1])*delta[2]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[2]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[2]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[2] , lower.tail = TRUE))
  beta3=INorm*((pnorm(sqrt(x[1])*delta[3]+x[3],0,1)-pnorm(sqrt(x[1])*delta[3]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[3]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[3]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[3] , lower.tail = TRUE))
  beta4=INorm*((pnorm(sqrt(x[1])*delta[4]+x[3],0,1)-pnorm(sqrt(x[1])*delta[4]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[4]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[4]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[4] , lower.tail = TRUE))
  beta5=INorm*((pnorm(sqrt(x[1])*delta[5]+x[3],0,1)-pnorm(sqrt(x[1])*delta[5]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[5]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[5]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[5] , lower.tail = TRUE))
  beta6=INorm*((pnorm(sqrt(x[1])*delta[6]+x[3],0,1)-pnorm(sqrt(x[1])*delta[6]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[6]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[6]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[6] , lower.tail = TRUE))
  beta7=INorm*((pnorm(sqrt(x[1])*delta[7]+x[3],0,1)-pnorm(sqrt(x[1])*delta[7]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[7]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[7]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[7] , lower.tail = TRUE))
  
  
  landa<-IT2*INE*c(.006527,.005757,.005091,.004484,.003966,.003493,.003078)+IT2*IUN*c(.003529,.003529,.003529,.003529,.003529,.003529,.003529)+IT2*IHN*c(.004678,.004572,.004394,.004158,.003875,.003556,.003213)+INB*INE*c(.004566,.003557,.003059,.002772,.002502,.002155,.001689)+INB*IUN*c(.002294,.002294,.002294,.002294,.002294,.002294,.002294)+INB*IHN*c(.004220,.003608,.003190,.002901,.002612,.002194,.001557)
  landa0<-sum(landa)  
  beta =(landa[1]/landa0)*beta1+(landa[2]/landa0)*beta2+(landa[3]/landa0)*beta3+(landa[4]/landa0)*beta4+(landa[5]/landa0)*beta5+(landa[6]/landa0)*beta6+(landa[7]/landa0)*beta7
  h<-numeric(3)
  h[1]<- 1-alpha
  h[2]<- 1-beta
  return(h)
}
o=auglag(par=c(20,5,10.2),ECTPHM,hin=hin.hs100, control.outer = list(method="nlminb",ilack.max=10))

##Or, Use an alternative Optimization function to optimize "ECTPHM", e.g. optim
o=optim(par=c(3,1,1.5),fn=ECTPHM,method = "L-BFGS-B",
        lower=c(1,.5,.5),upper = c(100,9,20))
##########
#To See the Optimal Design Parameters:
##########
o
x1=o$par[1];x2=o$par[2];x3=o$par[3]
x=c(x1,x2,x3)

##To obtain Type I Error Probability and the Power for T2 Chart
landa<-IT2*INE*c(.006527,.005757,.005091,.004484,.003966,.003493,.003078)+IT2*IUN*c(.003529,.003529,.003529,.003529,.003529,.003529,.003529)+IT2*IHN*c(.004678,.004572,.004394,.004158,.003875,.003556,.003213)+INB*INE*c(.004566,.003557,.003059,.002772,.002502,.002155,.001689)+INB*IUN*c(.002294,.002294,.002294,.002294,.002294,.002294,.002294)+INB*IHN*c(.004220,.003608,.003190,.002901,.002612,.002194,.001557)
landa0<-sum(landa)  

alpha<-INorm*(2*pnorm(-x[3],0,1))+INNorm*(1+(1/((1+(mb+x[3]*sb)^cb)^kb))-(1/((1+(mb-x[3]*sb)^cb)^kb)))+IT2*(1- pchisq(x[3], df=DF))
beta1=INorm*((pnorm(sqrt(x[1])*delta[1]+x[3],0,1)-pnorm(sqrt(x[1])*delta[1]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[1]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[1]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[1] , lower.tail = TRUE))
beta2=INorm*((pnorm(sqrt(x[1])*delta[2]+x[3],0,1)-pnorm(sqrt(x[1])*delta[2]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[2]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[2]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[2] , lower.tail = TRUE))
beta3=INorm*((pnorm(sqrt(x[1])*delta[3]+x[3],0,1)-pnorm(sqrt(x[1])*delta[3]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[3]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[3]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[3] , lower.tail = TRUE))
beta4=INorm*((pnorm(sqrt(x[1])*delta[4]+x[3],0,1)-pnorm(sqrt(x[1])*delta[4]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[4]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[4]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[4] , lower.tail = TRUE))
beta5=INorm*((pnorm(sqrt(x[1])*delta[5]+x[3],0,1)-pnorm(sqrt(x[1])*delta[5]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[5]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[5]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[5] , lower.tail = TRUE))
beta6=INorm*((pnorm(sqrt(x[1])*delta[6]+x[3],0,1)-pnorm(sqrt(x[1])*delta[6]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[6]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[6]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[6] , lower.tail = TRUE))
beta7=INorm*((pnorm(sqrt(x[1])*delta[7]+x[3],0,1)-pnorm(sqrt(x[1])*delta[7]-x[3],0,1)))+INNorm*((1/((1+(mb-x[3]*sb-sb*delta[7]*sqrt(x[1]))^cb)^kb))-(1/(1+(mb+x[3]*sb-sb*delta[7]*sqrt(x[1]))^cb)^kb))+IT2*(pchisq(x[3], df=DF, ncp =x[1]*NCP[7] , lower.tail = TRUE))
beta =(landa[1]/landa0)*beta1+(landa[2]/landa0)*beta2+(landa[3]/landa0)*beta3+(landa[4]/landa0)*beta4+(landa[5]/landa0)*beta5+(landa[6]/landa0)*beta6+(landa[7]/landa0)*beta7

alpha
1-beta
