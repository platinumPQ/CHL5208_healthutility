##############################################################################
###################### Competing Risk Model #############################
##############################################################################
library(dplyr)
library(tidyr)
library(readr)

## load data ----
library(readr)
df<- read_csv("D:/AAAAAA YUPEIQING/PhD/CHL5208/CHL5208_healthutility/df with age sex bm.csv")

library(cmprsk)
library(survival)

df$time_death <- df$`time to death or last follow up`      
df$time_brain <- df$`time to brain mets development`      

# Suppose brain mets is of interest, death is the competing risk
#############################################################################
## create event time and status ----
# event=0: censored; 
# event=1: the event of interest (brain mets); 
# event=2: a competing event (death without prior brain mets)
df$ftime <- pmin(df$time_brain, df$time_death, na.rm = TRUE)


df <- df  %>% 
  mutate(fstatus = case_when(# brain mets before or at death
                             !is.na(time_brain) & (time_brain <= time_death) ~ 1,
                             # death first or only death no brain mets
                             death & (is.na(time_brain) | time_brain > time_death) ~ 2,
                             # neither occurred
                             TRUE ~ 0
                             )
           )
table(df$fstatus)

## Cumulative incidence curves (nonparametric) ----
# probability of developing brain mets over time, accounting for competing risk death
ci <- cuminc(ftime = df$ftime, fstatus = df$fstatus, group = df$trt.group)
print(ci)
# p = 0.369: no significant difference between treatment groups in risk of brain mets.
# p = 0.088 : borderline difference in death risk between groups.

# plot
plot(ci, main = "Cumulative incidence: Brain mets (1) vs Death (2) by treatment")

## Fine-Gray model for brain mets (failcode = 1) ---
# include covariates: trt.group, age, sex
fg <- crr(ftime = df$ftime, fstatus = df$fstatus,
          cov1 = model.matrix(~ trt.group + age + sex, data = df)[,-1],
          failcode = 1)
summary(fg)


## Cause-specific Cox model for brain mets (treat competing events as censored) ----
df$event_cs_brain <- as.numeric(df$fstatus == 1)
cs_cox_brain <- coxph(Surv(ftime, event_cs_brain) ~ trt.group + age + sex, data = df)
summary(cs_cox_brain)
# Treatment B has 17% higher hazard of brain mets than A -> not significant

## QALY for competing risk ----
QALY_CIF <- function(CIF, U, tt, wt = 1, lam1 = 1, lam2 = 1) {
  
  # CIF : group-level vector of length tt for brain mets (failure type 1)
  #       1 - CIF = P(no brain mets by t)
  # U   : patient-by-time utility matrix 
  # wt  : weights (scalar or vector)
  
  if (length(wt) == 1) wt <- rep(wt, tt)
  
  # EU(t) = average utilities across subjects
  U_avg <- colMeans(U[,1:tt], na.rm = TRUE)
  
  # Use CIF instead of S_KM
  # CIF(t) = cumulative incidence for brain mets at time t
  # Prob no event by t = 1 - CIF(t)
  
  # Incorporate weights + lambda parameters
  # Marginal/group-level QALY
  Q <- sum( (U_avg^lam2) * ((1 - CIF)^lam1) * wt )
  
  return(Q)
}

# Survival replaced by CIF: S(t)=~1−CIF(t)

# predict CIF for each subject in Group A and take the mean ----
X <- model.matrix(~ trt.group + age + sex, data = df)[, -1] # model matrix as CRR
fg <- crr(ftime = df$ftime, fstatus = df$fstatus,cov1 = X,failcode = 1)
tt <- 36

# predict subject level CIF for all patients in Grp A
times <- 1:36 

## extract covariate rows belonging to Group A ----
covA_bar <- matrix(colMeans(X[df$trt.group == "A", ]), nrow = 1)

# Note!!! Fine–Gray does not naturally give CIF at every integer time
# predict.crr() returns CIF at event times, not at 1, 2, …, tt.

predA <- predict(fg, cov1 = covA_bar, times = 1:tt)
# Use right-continuous step function to impute those non-event time points:
# - before first event → 0
# - between events → stays constant
# - jumps at event times
# - after last event → stays at last value
time_ev <- predA[, 1]   # 19 event times
CIF_ev  <- predA[, 2]   # CIF at those times

CIF_fun <- stepfun(x = time_ev,y = c(0, CIF_ev))
CIF_A <- CIF_fun(1:36) #group level of length 36 (for each time points)

## extract covariate rows belonging to Group B ----
covB_bar <- matrix(colMeans(X[df$trt.group == "B", ]), nrow = 1)
predB <- predict(fg, cov1 = covB_bar, times = 1:tt)

# Use right-continuous step function to impute those non-event time points:
time_evB <- predB[, 1]   
CIF_evB  <- predB[, 2] 

CIF_funB <- stepfun(x = time_evB,y = c(0, CIF_evB))
CIF_B <- CIF_funB(1:36)

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

#utility group 1
U1b <- df %>% 
  filter(trt.group=="A") %>% 
  select(starts_with("U_T")) 

U1=t(apply(U1b,1,fillU))

#utility group 2
U2b <- df %>% 
  filter(trt.group=="B") %>% 
  select(starts_with("U_T")) 

U2=t(apply(U2b,1,fillU))

#calculate QALY based on survival and imputed U ----
tt=36

Q1a <- QALY_CIF(CIF_A, U1, tt)
Q2a <- QALY_CIF(CIF_B, U2, tt)

Q_obs <- Q1a - Q2a
Q_obs

n1=nrow(U1);n2=nrow(U2) #sample sizes

#permutation test (estimate null distribution)
Q_diff=c()
set.seed(1)
for(it in 1:1000){
  #permutation

  U=rbind(U1,U2)
  index1=sample(1:(n1+n2),n1,replace=FALSE)
  
  U1_perm=U[index1,]
  U2_perm=U[-index1,]
  
  predA <- predict(fg, cov1 = matrix(colMeans(X[index1, ]), nrow = 1), times = 1:36)
  CIF_funA <- stepfun(x = predA[, 1] ,y = c(0, predA[, 2]))
  CIF_A <- CIF_funA(1:36)
  
  predB <- predict(fg, cov1 = matrix(colMeans(X[-index1, ]), nrow = 1), times = 1:36)
  CIF_funB <- stepfun(x = predB[, 1] ,y = c(0, predB[, 2]))
  CIF_B <- CIF_funB(1:36)

  Q1a_perm=QALY_CIF(CIF_A,U1_perm,tt)
  Q2a_perm=QALY_CIF(CIF_B,U2_perm,tt)
  Q_diff[it]=Q1a_perm-Q2a_perm
}

#calculate p-value
pv=mean(abs(Q_diff)>=abs(Q_obs))
#Note: simulation might assumed a very large effect size to observe p-value=0


# Permutation tests shuffle the outcomes under the null; CIF model should remain fixed
# because recomputing it changes the expected difference and shrinks the null distribution
# -----> Not the correct way to do permutation test under null hypothesis: Please choose method 1!!
Q_diff2=c()
set.seed(1)
for(it in 1:1000){
  #permutation
  
  U=rbind(U1,U2)
  index1=sample(1:(n1+n2),n1,replace=FALSE)
  
  U1_perm=U[index1,]
  U2_perm=U[-index1,]
  
  predA <- predict(fg, cov1 = matrix(colMeans(X[df$trt.group=="A", ]), nrow = 1), times = 1:36)
  CIF_funA <- stepfun(x = predA[, 1] ,y = c(0, predA[, 2]))
  CIF_A <- CIF_funA(1:36)
  
  predB <- predict(fg, cov1 = matrix(colMeans(X[df$trt.group=="B", ]), nrow = 1), times = 1:36)
  CIF_funB <- stepfun(x = predB[, 1] ,y = c(0, predB[, 2]))
  CIF_B <- CIF_funB(1:36)
  
  Q1a_perm=QALY_CIF(CIF_A,U1_perm,tt)
  Q2a_perm=QALY_CIF(CIF_B,U2_perm,tt)
  Q_diff2[it]=Q1a_perm-Q2a_perm
}

#calculate p-value
pv2=mean(abs(Q_diff2)>=abs(Q_obs))

## This strongly rejects the null hypothesis --> QALY differs between Group A and B.


# Plot survival and CIF together ----
library(survival)
library(ggplot2)

# survival probability if brain mets free at each time -> descending over time
get_KM <- function(time, status, tt) {
  f <- survfit(Surv(time, status) ~ 1)
  s <- numeric(tt)
  for(t in 1:tt) s[t] <- if(t < f$time[1]) 1 else f$surv[max(which(f$time <= t))]
  s
}

# Build plotting data
df_plot <- data.frame(
  Time = rep(1:tt, 4),
  Probability = c(1-CIF_A, 1-CIF_B,
                  get_KM(df$ftime[df$trt.group=="A"], df$fstatus[df$trt.group=="A"]==1, tt),
                  get_KM(df$ftime[df$trt.group=="B"], df$fstatus[df$trt.group=="B"]==1, tt)),
  Group = rep(c("Treatment A","Treatment B","Treatment A","Treatment B"), each = tt),
  CurveType = rep(c("Fine-Gray","Fine-Gray","KM","KM"), each = tt)
)

# Plot
ggplot(df_plot, aes(Time, Probability, color = Group, linetype = CurveType)) +
  geom_step(size=1) +
  scale_linetype_manual(values=c("Fine-Gray"="solid","KM"="dashed")) +
  scale_color_manual(values=c("Treatment A"="#1b9e77","Treatment B"="#d95f02")) +
  labs(x="Time (months)", y="Survival Probability", color="Group", linetype="Curve",
       title = "Brain Metastases: CIF vs KM Survival") +
  theme_minimal()

