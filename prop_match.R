##############################################################################
###################### Propensity Score Matching #############################
##############################################################################
library(dplyr)
library(tidyr)
library(readr)

## load data ----
df <- read_csv("df with age sex.csv") %>% 
  mutate(trt=if_else(trt.group=="A",1,0))

# propensity score matching ----
#install.packages("MatchIt")
library(MatchIt)

# Propensity score matching to balance baseline characteristics between treatment groups
# PS=predicted prob of receiving trt -> match PS between trt and remove unmatched pairs

# View initial imbalance before matching
# Q: Baseline only U_T1??
m.out0 <- matchit(trt ~ age + sex + U_T1,
                  data = df,
                  method = NULL,
                  distance = "glm")
summary(m.out0) # Very Imbalanced!

# Nearest matching method
m.out <- matchit(
  trt ~ age + sex + U_T1,
  data = df,
  method = "nearest",
  ratio = 1
)
summary(m.out)
# Simplest model: everyone was kept -> paired each treated patient with a bad match.
# All SMDs are far too large -> treated and control patients are not comparable
# Good balance/SMD should be < 0.1 (absolute value)


# Add caliper to force good match only
m.out1 <- matchit(
  trt ~ age + sex + U_T1,
  data = df,
  method = "nearest",
  ratio = 1,
  caliper = 0.2
)
summary(m.out1)

# extract matched dataset -> remove 17 pairs
df_matched <- match.data(m.out1)

# Might consider other methods to match without removing: 
# IPWT(inverse prob weighting) <- adjust weight based on prob


## Missing data imputation only within subject ----
fillU <- function(Ut) {
  Uf <- Ut
  n <- length(Uf)
  
  whichna   <- which(is.na(Uf))      # indices of missing values
  whichgood <- which(!is.na(Uf))     # indices of observed values
  
  if(length(whichna) == 0) return(Uf) 
  
  for(hh in whichna){
    if(hh == 1){
      # if first value is missing -> use first observed
      Uf[hh] <- Uf[min(whichgood)]
    } else if(hh == n){
      # if last value is missing -> carry forward previous
      Uf[hh] <- Uf[hh-1]
    } else {
      # linear interpolation between previous and next observed
      future_idx <- (hh+1):n
      next_good  <- future_idx[future_idx %in% whichgood][1]
      
      if(is.na(next_good)){
        # no future observed value, carry forward previous
        Uf[hh] <- Uf[hh-1]
      } else {
        t1 <- hh-1
        t2 <- next_good
        Uf[hh] <- ((hh-t1)/(t2-t1)) * Uf[t2] + ((t2-hh)/(t2-t1)) * Uf[t1]
      }
    }
  }
  return(Uf)
}

library(survival)

# basic function to calculate HUS
QALY=function(S,U,tt,wt=1,lam1=1,lam2=1){
  SS=S
  SS$failed=as.numeric(SS$failed)
  f1 <- survfit(Surv(y, failed) ~ 1, data = S)
  
  # Extract survival times and survival probabilities
  S_time=f1$time
  S_KM0=f1$surv
  
  S_KM=c()
  for(t in 1:tt){
    if(t<S_time[1]){
      # If t is before the first event, survival probability = 1
      S_KM[t]=1
    }else{
      now=which(S_time<=t)
      if(length(now)==0){
        S_KM[t]=0
      }else{
        # Use the last KM estimate before or at time t
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
    # Compute average utility for subjects alive at this time
    U_avg=mean(U[stillalive,t])
    
    if(length(wt)==1){
      # If a single weight, use it uniformly
      QA=QA+U_avg^lam2*S_KM[t]^lam1
    }else{
      # If time-dependent weights, multiply by wt[t]
      QA=QA+U_avg^lam2*S_KM[t]^lam1*wt[t]
    }
  }
  return(QA)
}


#utility group 1
U1b <- df_matched %>% 
  filter(trt.group=="A") %>% 
  select(starts_with("U_T")) 
 
#utility group 2
U2b <- df_matched %>% 
  filter(trt.group=="B") %>% 
  select(starts_with("U_T"))  
  
S1=df[1:33,c(3:4)] #survival group 1
S2=df[34:66,c(3:4)] #survival group 2

n1=nrow(U1b);n2=nrow(U2b) #sample sizes

# two steps to impute U1b, U2b (utility for group 1; utility for group 2)
#U1b=fill_initial(U1b)
U1=t(apply(U1b,1,fillU))
#U2b=fill_initial(U2b)
U2=t(apply(U2b,1,fillU))

#calculate QALY based on survival and imputed U
tt=36
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

#######################################################################
#### use whole sample size instead of PSM ----
U1b <- df %>% 
  filter(trt.group=="A") %>% 
  select(starts_with("U_T")) 

#utility group 2
U2b <- df %>% 
  filter(trt.group=="B") %>% 
  select(starts_with("U_T"))  

S1=df[1:50,c(3:4)] #survival group 1
S2=df[51:100,c(3:4)] #survival group 2

n1=nrow(U1b);n2=nrow(U2b) #sample sizes

# two steps to impute U1b, U2b (utility for group 1; utility for group 2)
#U1b=fill_initial(U1b)
U1=t(apply(U1b,1,fillU))
#U2b=fill_initial(U2b)
U2=t(apply(U2b,1,fillU))

#calculate QALY based on survival and imputed U
tt=36
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
pv1=mean(abs(Q_diff)>=abs(Q_obs))

