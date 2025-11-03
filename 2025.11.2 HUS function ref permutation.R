
# initial fill: if a timepoint has valid data, use average to impute missing
# input U1b is the whole utility matrix (each row for a subject)

fill_initial=function(U1b,add_var=1){
  U1b2=U1b
  recorded_tp1=which(colSums(!is.na(U1b))>0)
  for(j in 1:length(recorded_tp1)){
    tu=recorded_tp1[j]
    lala=which(is.na(U1b[,tu]))
    if(length(lala)>0){
      if(add_var==0){
        U1b2[lala,tu]=mean(U1b2[,tu],na.rm=TRUE)
      }else{
        # add small random disturbance
        U1b2[lala,tu]=mean(U1b2[,tu],na.rm=TRUE)+rnorm(length(lala),0,sd(U1b[,tu],na.rm=TRUE))
        U1b2[lala,tu]=pmin(U1b2[lala,tu],1)
        U1b2[lala,tu]=pmax(U1b2[lala,tu],0)
      }
    }
  }
  return(U1b2)
}


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



# basic function to calculate HUS
QALY=function(S,U,tt,wt=1,lam1=1,lam2=1){
  SS=S
  SS$failed=as.numeric(SS$failed)
  f1 <- survfit(Surv(y, failed) ~ 1, data = S)
  
  S_time=f1$time
  S_KM0=f1$surv
  
  S_KM=c()
  for(t in 1:tt){
    if(t<S_time[1]){
      S_KM[t]=1
    }else{
      now=which(S_time<=t)
      if(length(now)==0){
        S_KM[t]=0
      }else{
        now=max(now)
        S_KM[t]=S_KM0[now]
      }
    }
  }
  
  
  QA=0
  for(t in 1:tt){
    stillalive=which(SS[,1]>t)
    if(length(stillalive)<1){
      next
    }
    U_avg=mean(U[stillalive,t])
    
    if(length(wt)==1){
      QA=QA+U_avg^lam2*S_KM[t]^lam1
    }else{
      QA=QA+U_avg^lam2*S_KM[t]^lam1*wt[t]
    }
  }
  return(QA)
}



U1b=df[1:50,-c(1:4)] #utility group 1
U2b=df[51:100,-c(1:4)] #utility group 2
S1=df[1:50,c(3:4)] #survival group 1
S2=df[51:100,c(3:4)] #survival group 2
n1=nrow(U1b);n2=nrow(U2b) #sample sizes

# two steps to impute U1b, U2b (utility for group 1; utility for group 2)
U1b=fill_initial(U1b)
U1=t(apply(U1b,1,fillU))
U2b=fill_initial(U2b)
U2=t(apply(U2b,1,fillU))

#calculate QALY based on survival and imputed U
Q1a=QALY(S1,U1,tt)
Q2a=QALY(S2,U2,tt)
Q_obs=Q1a-Q2a #observed difference




#permutation test (estimate null distribution)

Q_diff=c()

set.seed(1)
for(it in 1:1000){
  #permutation
  S=rbind(S1,S2)
  U=rbind(U1,U2)
  index1=sample(1:(n1+n2),n1,replace=FALSE)
  S1_perm=S[index1,]
  S2_perm=S[-index1,]
  U1_perm=U[index1,]
  U2_perm=U[-index1,]  
  
  Q1a_perm=QALY(S1_perm,U1_perm,tt)
  Q2a_perm=QALY(S2_perm,U2_perm,tt)
  Q_diff[it]=Q1a_perm-Q2a_perm
}


#calculate p-value
pv=mean(abs(Q_diff)>=abs(Q_obs))





