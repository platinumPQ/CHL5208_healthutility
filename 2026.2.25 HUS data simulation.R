


# linear fill: Ut is the vector of utility scores for one subject; it fills missing values linearly

fillU=function(Ut){
  Uf=Ut
  whichna=which(is.na(Uf))
  whichgood=which(!is.na(Uf))
  if(length(whichna)==0){
    return(Uf)
  }
  for(hh in whichna){
    if(hh==1){
      comeon=min(which(!is.na(Uf)))
      Uf[hh]=Uf[comeon]
    }else{
      if(hh==length(Uf)){
        Uf[hh]=Uf[hh-1]
      }else{
        houmiande=hh:length(Uf)
        if(mean(houmiande%in%whichna)==1){
          Uf[hh]=Uf[hh-1]
        }else{
          la1=hh-1
          la2=min(intersect(houmiande,whichgood))
          Uf[hh]=(hh-la1)/(la2-la1)*Uf[la2]+(la2-hh)/(la2-la1)*Uf[la1]
        }
      }
    }
  }
  return(Uf)
}



###################### simulate data

#sample size
n1=n2=50

tt=36   #months
tp=1:tt

set.seed(1)


#define utility functions
tcut=3 #changing point at 3 months

scenario="A1"

U1_base=U2_base=rep(NA,tt)
HR_true=1

if(scenario=="A1"){
  U1_base[1]=U2_base[1]=0.8
  U1_base[tcut]=0.5
  U2_base[tcut]=0.3
  U1_base[tt]=0.8
  U2_base[tt]=0.7
  p_missU=0.3 #probability of utility being missing
  xa=c(1,3,6,9,15,21,27,36) #timepoints when utility is recorded
}

U1_base=fillU(U1_base)
U2_base=fillU(U2_base)

#plot the baseline utility
#each individual's utility will be generated based on this, with added variations
plot(U1_base~tp,type="l",col=2,ylim=c(0,1),xlab="Month",ylab="Utility",lwd=2,main=paste("Scenario",scenario))
lines(U2_base~tp,col=4,xlab="Month",ylab="Utility",lwd=2)
legend("bottomright", legend=c("Treatment 1", "Treatment 2"),col=c(2,4), lty=1, cex=1,lwd=2)


#################parameter for covariates#######################

#cov -> treatment assignment
ba1=0.01    #older: choose more conservative (more trt2)
ba2=0.2    #male: more conservative (more trt2)

#cov -> survival
beta_age=0.04
beta_sex=0.2

#cov -> utility
theta_age1=-0.01  #older: worse
theta_sex1=-0.1  #male: worse
theta_age2=-0.01 #same for trt2
theta_sex2=-0.1

#trt2 has more males/older, and males/older tend to have worse survival and worse utility
#without adjusting, comparison will be biased

###########################################################


  #generate 10000 samples first
  set.seed(130)
  nn=10000
  age_min=35;age_max=80
  age_avg=mean(c(age_min,age_max))
  age=round(runif(nn,age_min,age_max))
  sex=rbinom(nn,1,0.5)    #1: male
  
  #generate treatments (trt2: old/control)
  prob_trt2=0.5+(age-age_avg)*ba1
  prob_trt2=prob_trt2+(sex-0.5)*ba2
  prob_trt2=pmin(prob_trt2,1)
  prob_trt2=pmax(prob_trt2,0)
  
  trt_all=rbinom(nn,1,prob_trt2)+1

  table(trt_all,sex) #trt 2 has more males
  mean(age[trt_all==1]);mean(age[trt_all==2]) #trt 2 has more older
  
  
  #select 100 for the dataset
  select1=sample(which(trt_all==1),n1)
  select2=sample(which(trt_all==2),n2)
  
  age1=age[select1]
  age2=age[select2]
  sex1=sex[select1]
  sex2=sex[select2]
  
  #simulate survival
  beta0=-4
  beta0b=beta0-log(HR_true)
  
  
  xz=beta_age*(mean(c(age1,age2)))+beta_sex*(mean(c(sex1,sex2)))
  
  rate1=exp(beta0)*exp(beta_age*(age1-age_avg)+beta_sex*(sex1-0.5))
  rate2=exp(beta0b)*exp(beta_age*(age2-age_avg)+beta_sex*(sex2-0.5))
  
  S1 <-rexp(n1, rate <-rate1)
  S2 <-rexp(n2, rate <-rate2)
  
  ## The censoring time

  zc=4.5
  
  C1<-pmin(runif(n1,min=0,max=zc*tt))
  C2<-pmin(runif(n2,min=0,max=zc*tt))

  ## The indicator
  delta1=delta10=(S1<=C1)
  delta2=delta20=(S2<=C2)
  
  ## The observed time
  S1<-pmin(S1,C1)
  S2<-pmin(S2,C2)
  
  delta1[which(S1>tt)]=FALSE
  delta2[which(S2>tt)]=FALSE
  
  S1=pmin(S1,tt)
  S2=pmin(S2,tt)
  
  S1=ceiling(S1)
  S2=ceiling(S2)
  
  S1=data.frame(y=S1,failed=delta1)
  S2=data.frame(y=S2,failed=delta2)
  
  
# Q: How to control the censor rate (say about 30% has censoring)?
  
###################


# simulate utility
U1=matrix(nrow=n1,ncol=tt)
U2=matrix(nrow=n2,ncol=tt)

ma1=theta_age1*mean(age1)+theta_sex1*mean(sex1)
ma2=theta_age2*mean(age2)+theta_sex2*mean(sex2)

for(i in 1:n1){
  uu=U1_base+theta_age1*(age1[i]-age_avg)+theta_sex1*(sex1[i]-0.5)+rnorm(tt,sd=0.1)
  uu=pmin(uu,1)
  uu=pmax(uu,0)
  U1[i,]=uu
}
for(i in 1:n2){
  uu=U2_base+theta_age2*(age2[i]-age_avg)+theta_sex2*(sex2[i]-0.5)+rnorm(tt,sd=0.1)
  uu=pmin(uu,1)
  uu=pmax(uu,0)
  U2[i,]=uu
}


# some utility scores are set to missing

missingrate=c(0)
missingrate[2:tt]=p_missU

U1b=matrix(nrow=n1,ncol=tt)
U2b=matrix(nrow=n2,ncol=tt)

for(j in 1:tt){
  hua=sample(c(1,NA),size=nrow(U1),prob=c(1-missingrate[j],missingrate[j]),replace=TRUE)
  U1b[,j]=hua*U1[,j]
  hua=sample(c(1,NA),size=nrow(U2),prob=c(1-missingrate[j],missingrate[j]),replace=TRUE)
  U2b[,j]=hua*U2[,j]
}



# S1, S2, U1b, U2b are the simulated data for HUS





# Q: We simulated the survival data using the exponential distributions, how to simulate data with competing risks?

# Q: Exponential is an easy choice, but there are other distributions like log-logistic, Weibull. What are their pros and cons?












