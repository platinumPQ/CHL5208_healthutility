simulate_data <- function(
  n0 = 100,
  n1 = 100,
  tt = 60,
  scenario = "D",
  
  # base utility type
  type = "recovery",
  
  # treatment assignment parameters
  alpha = c(-0.5, 0.005, 0.1, 0.3),
  # reduced confounding strength
  
  # survival parameters
  beta = c(beta0=-4, age=0.04, sex=0.2, sev=0.8),
  shape = 1.2,
  
  # utility parameters
  theta = c(age=-0.01, sex=-0.1, sev=-0.2),
  
  sigma_noise = 0.05,
  sigma_subject = 0.1,
  
  p_miss = 0.1,
  
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
  
  # survival benefit only (strong signal)
  if(scenario=="B") beta_trt_scn <- -1
  
  # utility benefit only
  if(scenario=="C") beta_trt_scn <- 0.05
  
  # early harm, late benefit (moderate better survival)
  if(scenario=="D") beta_trt_scn = -0.1
    
  # Trade-off (worse QoL, longer life)
  if(scenario=="E") beta_trt_scn <- -0.25
  
  # Trade-off (better QoL, shorter life)
  if(scenario=="F") beta_trt_scn <- 0.25
  
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
  
  f_base <- function(t, scenario){
    # Null - treatment does nothing
    if(scenario=="A") return(0.8) # flat base utility
    # survival benefit only
    if(scenario=="B") return(0.75) # flat base utility
    # utility benefit only
    if(scenario=="C") return(0.85 - 0.25*exp(-t/5)) # recovery base
    # early harm late benefit
    if(scenario=="D") return(0.85 - 0.25*exp(-t/5)) # recovery base
    # trade-off: longer life vs worse QoL
    if(scenario=="E") return(0.75) # flat base
    # trade-off: shorter life vs better QoL
    if(scenario=="F") return(0.85 - 0.25*exp(-t/5)) # recovery base
  }
  
  f_trt <- function(t, scenario){
    if(scenario == "A") return(0)
    if(scenario == "B") return(0)
    if(scenario == "C") return(0.15)
    if(scenario == "D"){
      if(t <= 6){
        return(-0.10)
      } else {
        return(0.15)
      }
    }
    if(scenario == "E") return(-0.15)
    if(scenario == "F") return(0.15)
  }
  
  for(i in 1:n1){
    subj_eff = rnorm(1, 0, sigma_subject)
    
    for(t in 1:tt){
      
      trt_eff <- f_trt(t, scenario)
      base <- f_base(t, scenario)
      
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
      
      trt_eff = 0 # no trt effect for control group
      base <- f_base(t, scenario)
      
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
  cat("Mean severity trt:", mean(sev1), "\n")
  cat("Mean severity ctl:", mean(sev0), "\n")
  
  cat("Mean age trt:", mean(age1), "\n")
  cat("Mean age ctl:", mean(age0), "\n")
  
  return(list(
    S_trt = S1, #treatment
    S_ctl = S0,
    U_trt = U1b,
    U_ctl = U0b,
    age_trt = age1,
    age_ctl = age0,
    sex_trt = sex1,
    sex_ctl = sex0,
    sev_trt = sev1,
    sev_ctl = sev0
  ))
}

######################################################
####### IPTW to address confounding ##########
########################################################
# In the presence of confounding, naïve comparisons produced biased estimates of treatment effect, 
# for Scenario B treatment appeared harmful due to imbalance in disease severity <-
# Apply IPTW to recover the true casual effect
get_weights <- function(age, sex, sev, trt){
  
  df = data.frame(trt = trt, age = age, sex = sex, sev = sev)
  
  fit = glm(trt ~ age + sex + sev, data = df, family = binomial)
  ps  = predict(fit, type = "response")
  
  # stabilized weights
  p_trt = mean(trt)
  
  w = ifelse(trt == 1,
             p_trt / ps,
             (1 - p_trt) / (1 - ps))
  
  # optional trimming (recommended)
  w = pmin(w, 10)
  
  return(w)
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
# QALY=function(S,U,tt,wt=1,lam1=1,lam2=1){
#   SS=S
#   f1 <- survfit(Surv(time, event) ~ 1, data = S)
#   
#   # Extract survival times and survival probabilities
#   S_time=f1$time
#   S_KM0=f1$surv
#   
#   S_KM=c()
#   for(t in 1:tt){
#     if(t<S_time[1]){
#       # If t is before the first event, survival probability = 1
#       S_KM[t]=1
#     }else{
#       now=which(S_time<=t)
#       if(length(now)==0){
#         S_KM[t]=0
#       }else{
#         # Use the last KM estimate before or at time t
#         now=max(now)
#         S_KM[t]=S_KM0[now]
#       }
#     }
#   }
#   
#   
#   QA=0
#   for(t in 1:tt){
#     stillalive=which(SS[,1]>t)
#     if(length(stillalive)<1){
#       next
#     }
#     # Compute average utility for subjects alive at this time
#     U_avg=mean(U[stillalive,t], na.rm=TRUE)
#     
#     if(length(wt)==1){
#       # If a single weight, use it uniformly
#       QA=QA+U_avg^lam2*S_KM[t]^lam1
#     }else{
#       # If time-dependent weights, multiply by wt[t]
#       QA=QA+U_avg^lam2*S_KM[t]^lam1*wt[t]
#     }
#   }
#   return(QA)
# }

## Modify QAL to weighted version ----
QALY_IPTW=function(S,U,w,tt,wt=1,lam1=1,lam2=1){
  
  # weighted KM
  f1 <- survfit(Surv(time, event) ~ 1, data = S, weights = w)
  
  S_time=f1$time
  S_KM0=f1$surv
  
  S_KM=numeric(tt)
  
  for(t in 1:tt){
    if(length(S_time)==0 || t < S_time[1]){
      S_KM[t]=1
    }else{
      now=which(S_time<=t)
      if(length(now)==0){
        S_KM[t]=1
      }else{
        S_KM[t]=S_KM0[max(now)]
      }
    }
  }
  
  QA=0
  
  for(t in 1:tt){
    stillalive=which(S$time > t)
    if(length(stillalive)<1) next
    
    # weighted mean utility
    U_avg = sum(w[stillalive] * U[stillalive,t], na.rm=TRUE) /
      sum(w[stillalive])
    
    if(length(wt)==1){
      QA = QA + U_avg^lam2 * S_KM[t]^lam1
    }else{
      QA = QA + U_avg^lam2 * S_KM[t]^lam1 * wt[t]
    }
  }
  
  return(QA)
}

## permutation test ----

n0 <- 100
n1 <- 100
tt <- 60

scenarios = c("A","B","C","D","E","F")

results = data.frame(
  scenario = scenarios,
  Q_obs = NA,
  p_value = NA
)

for(s in scenarios){
  
  dat = simulate_data(scenario = s)
  
  S_trt = dat$S_trt
  S_ctl = dat$S_ctl
  U_trt = dat$U_trt
  U_ctl = dat$U_ctl
  
  age_trt = dat$age_trt
  age_ctl = dat$age_ctl
  sex_trt = dat$sex_trt
  sex_ctl = dat$sex_ctl
  sev_trt = dat$sev_trt
  sev_ctl = dat$sev_ctl
  
  # ---- IPTW weights ----
  age = c(age_trt, age_ctl)
  sex = c(sex_trt, sex_ctl)
  sev = c(sev_trt, sev_ctl)
  trt = c(rep(1,n1), rep(0,n0))
  
  w = get_weights(age, sex, sev, trt)
  
  w_trt = w[trt==1]
  w_ctl = w[trt==0]
  
  cat("Scenario:", s, "\n")
  cat("Naive survival diff:", mean(S_trt$time) - mean(S_ctl$time), "\n")
  
  # Q_obs  = difference under REAL grouping (treatment vs control)
  # Q_diff = differences under RANDOM grouping (null world)  
  # ---- observed ----
  U_trt_imp = t(apply(U_trt, 1, fillU))
  U_ctl_imp = t(apply(U_ctl, 1, fillU))
  
  Q_trt_obs = QALY_IPTW(S_trt, U_trt_imp, w_trt, tt)
  Q_ctl_obs = QALY_IPTW(S_ctl, U_ctl_imp, w_ctl, tt)
  Q_obs = Q_trt_obs - Q_ctl_obs
  
  cat("Weighted survival diff:",
      weighted.mean(S_trt$time, w_trt) -
        weighted.mean(S_ctl$time, w_ctl), "\n")
  cat("Before weighting:\n")
  cat("Mean severity trt (trt vs ctl):", mean(sev_trt), mean(sev_ctl), "\n")
  
  cat("After weighting:\n")
  cat("Mean severity trt (trt vs ctl):", weighted.mean(sev_trt, w_trt),
      weighted.mean(sev_ctl, w_ctl), "\n")
  
  # ---- permutation ----

  Q_diff = numeric(500)  # can reduce for speed
  
  S = rbind(S_trt, S_ctl)
  U = rbind(U_trt, U_ctl)
  
  for(it in 1:500){
    
    idx = sample(1:(n0+n1), n1)
    
    S1p = S[idx, ]
    S2p = S[-idx, ]
    
    U1p = U[idx, ]
    U2p = U[-idx, ]
    
    w1p = w[idx]
    w2p = w[-idx]
    
    U1p_imp = t(apply(U1p, 1, fillU))
    U2p_imp = t(apply(U2p, 1, fillU))
    
    Q1p = QALY_IPTW(S1p, U1p_imp, w1p, tt)
    Q2p = QALY_IPTW(S2p, U2p_imp, w2p, tt)
    
    Q_diff[it] = Q1p - Q2p
  }
  
  pv = mean(abs(Q_diff) >= abs(Q_obs))
  
  results[results$scenario==s,] = c(s, Q_obs, pv)

}



#####################################################################
################# Analysis Step ######################################
#######################################################################
surv_test <- function(S1, S2){
  S = rbind(S1, S2)
  trt = c(rep(1,nrow(S1)), rep(0,nrow(S2)))
  
  fit = survdiff(Surv(time,event) ~ trt)
  pval = 1 - pchisq(fit$chisq,1)
  
  return(pval)
}

utility_test <- function(U1, U2){
  
  U1_imp = t(apply(U1,1,fillU))
  U2_imp = t(apply(U2,1,fillU))
  
  mean1 = mean(U1_imp)
  mean2 = mean(U2_imp)
  
  diff_obs = mean1 - mean2
  
  U = rbind(U1_imp, U2_imp)
  n1 = nrow(U1_imp)
  
  diffs = numeric(500)
  
  for(i in 1:500){
    idx = sample(1:nrow(U), n1)
    diffs[i] = mean(U[idx,]) - mean(U[-idx,])
  }
  
  mean(abs(diffs) >= abs(diff_obs))
}

# final result table
results = data.frame(
  scenario = scenarios,
  surv_p = NA,
  util_p = NA,
  qaly_p = NA
)



###################################################################
# IPTW weighted
surv_test <- function(S1, S2, w1, w2){
  
  S = rbind(S1, S2)
  
  S$trt = c(rep(1, nrow(S1)),
            rep(0, nrow(S2)))
  
  w = c(w1, w2)
  
  fit = coxph(
    Surv(time, event) ~ trt,
    data = S,
    weights = w
  )
  
  pval = summary(fit)$coef[,"Pr(>|z|)"]
  
  return(pval)
}

utility_test <- function(U1, U2, w1, w2, nperm = 500){
  
  U1_imp = t(apply(U1,1,fillU))
  U2_imp = t(apply(U2,1,fillU))
  
  # observed weighted means
  mean1 = weighted.mean(as.vector(U1_imp),
                        rep(w1, ncol(U1_imp)),
                        na.rm=TRUE)
  
  mean2 = weighted.mean(as.vector(U2_imp),
                        rep(w2, ncol(U2_imp)),
                        na.rm=TRUE)
  
  diff_obs = mean1 - mean2
  
  # pooled
  U = rbind(U1_imp, U2_imp)
  w = c(w1, w2)
  
  n1 = nrow(U1_imp)
  
  diffs = numeric(nperm)
  
  for(i in 1:nperm){
    
    idx = sample(1:nrow(U), n1)
    
    diffs[i] =
      weighted.mean(as.vector(U[idx,]),
                    rep(w[idx], ncol(U)),
                    na.rm=TRUE) -
      
      weighted.mean(as.vector(U[-idx,]),
                    rep(w[-idx], ncol(U)),
                    na.rm=TRUE)
  }
  
  pval = mean(abs(diffs) >= abs(diff_obs))
  
  return(pval)
}

results = data.frame(
  scenario = scenarios,
  surv_p = NA,
  util_p = NA,
  qaly_p = NA,
  qaly_diff = NA
)

for(s in scenarios){
  
  dat = simulate_data(scenario = s)
  
  # data
  S_trt = dat$S_trt
  S_ctl = dat$S_ctl
  
  U_trt = dat$U_trt
  U_ctl = dat$U_ctl
  
  # covariates
  age = c(dat$age_trt, dat$age_ctl)
  sex = c(dat$sex_trt, dat$sex_ctl)
  sev = c(dat$sev_trt, dat$sev_ctl)
  
  trt = c(rep(1,n1), rep(0,n0))
  
  # IPTW
  w = get_weights(age, sex, sev, trt)
  
  w_trt = w[trt==1]
  w_ctl = w[trt==0]
  
  ###################################################
  # Survival test
  ###################################################
  
  surv_p =
    surv_test(S_trt, S_ctl,
              w_trt, w_ctl)
  
  ###################################################
  # Utility test
  ###################################################
  
  util_p =
    utility_test(U_trt, U_ctl,
                 w_trt, w_ctl)
  
  ###################################################
  # QALY test
  ###################################################
  
  U_trt_imp = t(apply(U_trt,1,fillU))
  U_ctl_imp = t(apply(U_ctl,1,fillU))
  
  Q1 =
    QALY_IPTW(S_trt,
              U_trt_imp,
              w_trt,
              tt)
  
  Q0 =
    QALY_IPTW(S_ctl,
              U_ctl_imp,
              w_ctl,
              tt)
  
  Q_obs = Q1 - Q0
  
  # permutation
  Q_diff = numeric(500)
  
  S = rbind(S_trt, S_ctl)
  U = rbind(U_trt_imp, U_ctl_imp)
  
  for(it in 1:500){
    
    idx = sample(1:(n0+n1), n1)
    
    Q1p =
      QALY_IPTW(S[idx,],
                U[idx,],
                w[idx],
                tt)
    
    Q0p =
      QALY_IPTW(S[-idx,],
                U[-idx,],
                w[-idx],
                tt)
    
    Q_diff[it] = Q1p - Q0p
  }
  
  qaly_p =
    mean(abs(Q_diff) >= abs(Q_obs))
  
  ###################################################
  # store
  ###################################################
  
  results[results$scenario==s,] =
    list(s,
         surv_p,
         util_p,
         qaly_p,
         Q_obs)
}


