simulate_data <- function(
  n0 = 50,
  n1 = 50,
  tt = 36,
  scenario = "D",
  
  # base utility type
  type = "recovery",
  
  # treatment assignment parameters
  alpha = c(-0.5, 0.01, 0.2, 0.8),
  
  # survival parameters
  beta = c(beta0=-4, age=0.04, sex=0.2, sev=0.8),
  shape = 1.2,
  
  # utility parameters
  theta = c(age=-0.01, sex=-0.1, sev=-0.2),
  
  sigma_noise = 0.1,
  sigma_subject = 0.15,
  
  p_miss = 0.3,
  
  seed=123
  ){
  set.seed(seed)
  ##############################
  #### 1. Generate covariates
  ##############################
  
  n = 10000  # large pool to induce confounding
  
  age = runif(n, 35, 80)
  sex = rbinom(n, 1, 0.5)
  severity = rnorm(n, 0, 1)
  
  ##############################
  #### 2. Treatment assignment
  ##############################
  
  linpred = alpha[1] + alpha[2]*age + alpha[3]*sex + alpha[4]*severity
  prob_trt1 = plogis(linpred)
  trt_all = rbinom(n, 1, prob_trt1)
  
  # select balanced sample
  id0 = sample(which(trt_all==0), n0)
  id1 = sample(which(trt_all==1), n1)
  
  ##############################
  #### 3. Subset covariates
  ##############################
  
  age0 = age[id0]; age1 = age[id1]
  sex0 = sex[id0]; sex1 = sex[id1]
  sev0 = severity[id0]; sev1 = severity[id1]
  
  ##############################
  #### 4. Survival (Weibull)
  ##############################
  
  # define scenario specific survival effect
  beta_trt_scn <- 0
  
  # Null: treatment does nothing
  if(scenario=="A") beta_trt_scn <- 0
  
  # survival benefit only
  if(scenario=="B") beta_trt_scn <- -0.4
  
  # utility benefit only
  if(scenario=="C") beta_trt_scn <- 0
  
  # early harm, late benefit (moderate better survival)
  if(scenario=="D") beta_trt_scn <- -0.3
  
  # Trade-off (worse QoL, longer life)
  if(scenario=="E") beta_trt_scn <- -0.4
  
  # Trade-off (better QoL, shorter life)
  if(scenario=="F") beta_trt_scn <- 0.4
  
  # linear predictor
  lin0 = beta["beta0"] + beta["age"]*age0 + beta["sex"]*sex0 +
    beta["sev"]*sev0 + beta_trt_scn*0
  
  lin1 = beta["beta0"] + beta["age"]*age1 + beta["sex"]*sex1 +
    beta["sev"]*sev1 + beta_trt_scn*1
  
  lambda0 = exp(lin0)
  lambda1 = exp(lin1)
  
  T1 = rweibull(n1, shape=shape, scale=1/lambda1)
  T0 = rweibull(n0, shape=shape, scale=1/lambda0)
  
  # censoring
  C1 = runif(n1, 0, tt)
  C0 = runif(n0, 0, tt)
  
  time1 = pmin(T1, C1)
  time0 = pmin(T0, C0)
  
  event1 = as.numeric(T1 <= C1)
  event0 = as.numeric(T0 <= C0)
  
  time1 = ceiling(pmin(time1, tt))
  time0 = ceiling(pmin(time0, tt))
  
  S1 = data.frame(time = time1, event = event1)
  S0 = data.frame(time = time0, event = event0)
  
  ##############################
  #### 5. Utility trajectories
  ##############################

  U1 = matrix(NA, n1, tt)
  U0 = matrix(NA, n0, tt)
  
  for(i in 1:n1){
    subj_eff = rnorm(1, 0, sigma_subject)
    
    for(t in 1:tt){
      
      trt_eff = 0
      base <- 0.75
      
      if(scenario=="A"){ 
        # Null - treatment does nothing
        trt_eff = 0
        # flat base utility
        base <- 0.8
        }
      if(scenario=="B"){ 
        # survival benefit only
        trt_eff = 0
        # flat base utility
        base <- 0.75}
      if(scenario=="C"){ 
        # utility benefit only
        trt_eff = 0.1 
        # recovery base
        base <- 0.85 - 0.25*exp(-t/5)}
      if(scenario=="D"){ 
        # early harm late benefit
        if(t<=3) trt_eff = -0.1 else trt_eff = 0.08
        # recovery base utility
        base <- 0.85 - 0.25*exp(-t/5)}
      if(scenario=="E"){ 
        # trade-off: longer life vs worse QoL
        trt_eff = -0.15
        # flat base utility
        base <- 0.75}
      if(scenario=="F"){ 
        # trade-off: shorter life vs better QoL
        trt_eff = 0.15
        # recovery base
        base <- 0.85 - 0.25*exp(-t/5)}
      
      U1[i,t] = base +
        theta["age"]*(age1[i]-mean(age)) +
        theta["sex"]*(sex1[i]-0.5) +
        theta["sev"]*sev1[i] +
        trt_eff +
        subj_eff +
        rnorm(1,0,sigma_noise)
      
      U1[i,t] = min(max(U1[i,t],0),1)
      
      if(t > S1$time[i] & S1$event[i]==1) U1[i,t] = 0
    }
  }
  
  for(i in 1:n0){
    subj_eff = rnorm(1, 0, sigma_subject)
    
    for(t in 1:tt){
      
      base = 0.8 - 0.3*exp(-t/6)
      
      trt_eff = 0
      
      if(scenario=="A"){ if(t<=3) trt_eff = -0.1 }
      if(scenario=="B"){ trt_eff = 0 }
      if(scenario=="C"){ trt_eff = 0.1 }
      if(scenario=="D"){ if(t<=3) trt_eff = -0.2 else trt_eff = 0.05 }
      if(scenario=="E"){ if(t<=6) trt_eff = -0.15 else trt_eff = 0.15 }
      
      U0[i,t] = base +
        theta["age"]*(age0[i]-mean(age)) +
        theta["sex"]*(sex0[i]-0.5) +
        theta["sev"]*sev0[i] +
        trt_eff +
        subj_eff +
        rnorm(1,0,sigma_noise)
      
      U0[i,t] = min(max(U0[i,t],0),1)
      
      if(t > S0$time[i] & S0$event[i]==1) U0[i,t] = 0
    }
  }
  
  ##############################
  #### 6. Missingness
  ##############################
  
  U1b = U1
  U0b = U0
  
  for(j in 1:tt){
    miss1 = rbinom(n1,1,p_miss)
    miss0 = rbinom(n0,1,p_miss)
    
    U1b[miss1==1,j] = NA
    U0b[miss0==1,j] = NA
  }
  
  ##############################
  #### 7. Return
  ##############################
  
  return(list(
    S1 = S0,
    S2 = S1,
    U1b = U0b,
    U2b = U1b
  ))
}


#######################################################
####### Permutation test across all scenarios##########
########################################################

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
  f1 <- survfit(Surv(time, event) ~ 1, data = S)
  
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
    U_avg=mean(U[stillalive,t], na.rm=TRUE)
    
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


n0 <- 50
n1 <- 50
tt <- 36

scenarios = c("A","B","C","D","E","F")

results = data.frame(
  scenario = scenarios,
  Q_obs = NA,
  p_value = NA
)

for(s in scenarios){
  
  dat = simulate_data(scenario = s)
  
  S1 = dat$S1
  S2 = dat$S2
  U1 = dat$U1b
  U2 = dat$U2b
  
  # ---- observed ----
  U1_imp = t(apply(U1, 1, fillU))
  U2_imp = t(apply(U2, 1, fillU))
  
  Q1_obs = QALY(S1, U1_imp, tt)
  Q2_obs = QALY(S2, U2_imp, tt)
  Q_obs = Q1_obs - Q2_obs
  
  # ---- permutation ----
  Q_diff = numeric(500)  # can reduce for speed
  
  S = rbind(S1, S2)
  U = rbind(U1, U2)
  
  for(it in 1:500){
    
    idx = sample(1:(n0+n1), n1)
    
    S1p = S[idx, ]
    S2p = S[-idx, ]
    
    U1p = U[idx, ]
    U2p = U[-idx, ]
    
    U1p_imp = t(apply(U1p, 1, fillU))
    U2p_imp = t(apply(U2p, 1, fillU))
    
    Q1p = QALY(S1p, U1p_imp, tt)
    Q2p = QALY(S2p, U2p_imp, tt)
    
    Q_diff[it] = Q1p - Q2p
  }
  
  pv = mean(abs(Q_diff) >= abs(Q_obs))
  
  results[results$scenario==s,] = c(s, Q_obs, pv)
}
