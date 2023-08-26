library(dplyr)    
library(data.table)
library(tidyr)
library(reshape2) 
library(ggplot2) 
library(ggrepel)
library(ellipse)
library(gridExtra)
library(ggthemes)   # For colorblind palettes
library(scales)     # For dollar signs and commas 
library(boot)
library(dampack)    # For CEA and PSA visualization functionality
library(darthtools) # For WCC, parameter transformation an matrix checks
library(doParallel)
library(hrbrthemes)
library(dampack)
library(tibble)
library(stats)
library(haven)

#Reisdential
rm(list = ls())   
seed=071818
set.seed(seed)

#Get Exposure Probability
# residential_pest<-read.csv("Data_Exposure/r_avg_alltime.csv")
residential_pest <- read_sas("Data_Exposure/r_avg_lag10.sas7bdat")
residential_pest<-residential_pest[which(residential_pest$PD==0),]

PD_onset<-read.csv("Data_Probability/Transition_Prob_PD.csv")
PD_onset<-PD_onset[PD_onset$Usage=="Insecticide",]
PD_onset<-PD_onset[PD_onset$Banning=="Restricted use",]

v_names_str <- c(PD_onset$Chemcode) 
residential_pest<-residential_pest[,which(colnames(residential_pest) %in% v_names_str)]

exposure_prob<-as.data.frame(apply(residential_pest,2,function(x) sum(x > 0)))
exposure_prob <- tibble::rownames_to_column(exposure_prob, "Chemcode")
exposure_prob$prob<-exposure_prob$`apply(residential_pest, 2, function(x) sum(x > 0))`/nrow(residential_pest)

both12<-sum(residential_pest$Chem383 > 0 & residential_pest$Chem216 > 0)
both13<-sum(residential_pest$Chem383 > 0 & residential_pest$Chem105 > 0)
both23<-sum(residential_pest$Chem216 > 0 & residential_pest$Chem105 > 0)
both123<-sum(residential_pest$Chem383 > 0 & residential_pest$Chem216 > 0 & residential_pest$Chem105 > 0)

# sort by Chemcode
#add multiple exposure
exposure_prob[8,1]<-"Chem383*Chem216"
exposure_prob[8,2]<-both12
exposure_prob[8,3]<-both12/nrow(residential_pest)

exposure_prob[9,1]<-"Chem383*Chem105"
exposure_prob[9,2]<-both13
exposure_prob[9,3]<-both13/nrow(residential_pest)

exposure_prob[10,1]<-"Chem216*Chem105"
exposure_prob[10,2]<-both23
exposure_prob[10,3]<-both23/nrow(residential_pest)

exposure_prob[11,1]<-"Chem383*Chem216*Chem105"
exposure_prob[11,2]<-both123
exposure_prob[11,3]<-both123/nrow(residential_pest)

exposure_prob[12,1]<-"no pest"
exposure_prob[12,2]<-0
exposure_prob[12,3]<-0

exposure_prob <- exposure_prob[order(exposure_prob$Chemcode),]


residential_para<-read.csv("Data_Probability/psa_residential.csv")
# sort by Chemcode
PD_onset[8,1]<-"Chem383*Chem216"
PD_onset[8,2]<-"Oxydemeton-Methyl*Dimethoate"
PD_onset[8,3]<-1.198*1.138
PD_onset[8,4]<-1.21100*1.208
PD_onset[8,5]<-"Insecticide"

PD_onset[9,1]<-"Chem383*Chem105"
PD_onset[9,2]<-"Oxydemeton-Methyl*Carbaryl"
PD_onset[9,3]<-1.198*1.119
PD_onset[9,4]<-1.21100*1.168
PD_onset[9,5]<-"Insecticide"

PD_onset[10,1]<-"Chem216*Chem105"
PD_onset[10,2]<-"Dimethoate*Carbaryl"
PD_onset[10,3]<-1.138*1.119
PD_onset[10,4]<-1.208*1.168
PD_onset[10,5]<-"Insecticide"

PD_onset[11,1]<-"Chem383*Chem216*Chem105"
PD_onset[11,2]<-"Oxydemeton-Methyl*Dimethoate*Carbaryl"
PD_onset[11,3]<-1.198000*1.138*1.119
PD_onset[11,4]<-1.21100*1.208*1.168
PD_onset[11,5]<-"Insecticide"

PD_onset[12,1]<-"no pest"
PD_onset[12,2]<-"no pest"
PD_onset[12,3]<-1
PD_onset[12,4]<-1
PD_onset[12,5]<-"no pest"

PD_onset <- PD_onset[order(PD_onset$Chemcode),]
v_names_str <- c(PD_onset$Chemcode) 
residential_para<-residential_para[which(residential_para$chemcode %in% v_names_str),]

residential_para <- residential_para[order(residential_para$chemcode),]
residential_para[12,1]<-"no pest"
residential_para[12,2]<-0
residential_para[12,3]<-0

## General setup
cycle_length <- 1       # cycle length equal to one year (use 1/12 for monthly)
n_age_init <- 65        # age at baseline
n_age_max  <- 85       # maximum age of follow up
n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
v_names_states <- c("H",  # the 4 health states of the model:
                    "S1", # Healthy Unexposed (H), Healthy Exposed (S1), PD Onset (S2), Dead (D)
                    "S2",
                    "D")
n_states    <- length(v_names_states)     # number of health states 

# Discounting factors
d_c <- 0.03 # annual discount rate for costs 
d_e <- 0.03 # annual discount rate for QALYs

# Strategies
v_names_str <- c(PD_onset$Chemcode) 
n_str       <- length(v_names_str)        # number of strategies

# Within-cycle correction (WCC) using Simpson's 1/3 rule
v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                             method = "Simpson1/3") # vector of wcc

n_sim=1000
# Given mean and standard deviation
mean_val <- 23041
sd_val <- 34045

# Calculate shape (k) and scale (θ) parameters
shape_param <- (mean_val / sd_val)^2
scale_param <- sd_val^2 / mean_val
set.seed(seed)
# Generate random data from the gamma distribution
gamma_data <- rgamma(n = n_sim, shape = shape_param, scale = scale_param)

# Transform the data to achieve left-skewness
c_S2_psa <- abs(gamma_data - 2 * mean_val)
c_S2_psa=as.matrix(c_S2_psa)

# Check the results
mean(c_S2_psa)
sd(c_S2_psa)
hist(c_S2_psa)


# Utilities 
set.seed(seed)
u_S2    = runif(n_sim, min = 0.61, max = 0.74)
u_S2_psa=as.matrix(u_S2)

r_HS2  <- 0.0052
df_psa_1 <- data.frame(matrix(ncol = n_str, nrow = 1000))
colnames(df_psa_1)<-v_names_str

set.seed(seed)
for(m in 1:n_str) {
  r_S1S2  = rlnorm(n_sim, meanlog = residential_para[m,2],  sdlog = residential_para[m,3])*r_HS2
  df_psa_1[,m]<-r_S1S2
}

df_cea_final<-lapply(1:1000, matrix, data= NA, nrow=n_str, ncol=2)
df_cea_finalprop<-lapply(1:1000, matrix, data= NA, nrow=n_str, ncol=4)

for(j in 1:1000) {
  ## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
  r_S1S2 <- as.numeric(df_psa_1[j,])
  r_HS2  <- 0.0052
  r_HD   <- 0.0418
  r_S1D  <- 0.0418
  r_S2D  <- 0.0418
  ## State rewards
  # Costs
  c_H_nopest       <- 0.3      # annual cost of being Healthy unexposed
  c_S1_nopest      <- 0.3         # annual cost of being Healthy exposed
  c_S2_nopest      <- 0.3+as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  
  c_H       <- 0         # annual cost of being Healthy unexposed
  c_S1      <- 0         # annual cost of being Healthy exposed
  c_S2      <- as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  c_D       <- 0         # annual cost of being dead
  
  # Utilities
  u_H    <- 1           # annual utility of being Healthy unexposed
  u_S1   <- 1            # annual utility of being Healthy exposed
  u_S2   <- as.numeric(u_S2_psa[j,])     # annual utility of being PD Onset
  u_D    <- 0            # annual utility of being Dead
  
  # Discount weight for costs and effects
  v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  
  ####################### Construct state-transition models ######################
  ## Initial state vector
  # All starting healthy
  #California Population 65+ 5976000
  v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
  
  ## Initialize cohort trace for SoC
  m_M_NoPest <- matrix(NA,
                       nrow = (n_cycles + 1), ncol = n_states, 
                       dimnames = list(0:n_cycles, v_names_states))
  # Store the initial state vector in the first row of the cohort trace
  m_M_NoPest[1, ] <- v_m_init
  ## Initialize cohort trace for strategies A, B, and AB
  # Structure and initial states are the same as for SoC
  pest_cohort <- lapply(seq_len(n_str), function(X) m_M_NoPest)
  
  for (X in 1:n_str) {
    pest_cohort[[X]][1,"S1"]<- exposure_prob$prob[X]
    pest_cohort[[X]][1,"H"]<- 1-exposure_prob$prob[X]
    
  }
  
  ## Initialize transition probability matrix for Pesticide Ban
  # all transitions to a non-death state are assumed to be conditional on survival 
  m_P_NoPest <- matrix(0, 
                       nrow = n_states, ncol = n_states, 
                       dimnames = list(v_names_states, 
                                       v_names_states)) # define row and column names
  
  ## Initialize transition probability for Pesticide
  # Update transition probabilities
  pest <- lapply(seq_len(n_str), function(X) m_P_NoPest)
  
  ## Fill in matrix
  for (x in 1:n_str) {
    pest[[x]]["H", "H"]   <- 1 -r_HS2-r_HD
    pest[[x]]["H", "S2"]  <- r_HS2
    pest[[x]]["H", "D"]   <- r_HD
    # From S1
    pest[[x]]["S1", "S1"] <- 1 -r_S1S2[x]- r_S1D
    pest[[x]]["S1", "S2"] <- r_S1S2[x]
    pest[[x]]["S1", "D"]  <- r_S1D
    # From S2
    pest[[x]]["S2", "S2"] <- 1 - r_S2D
    pest[[x]]["S2", "D"]  <- r_S2D
    # From D
    pest[[x]]["D", "D"]   <- 1
  }
  
  ####################### Run Markov model #######################
  # Iterative solution of time-independent cSTM
  for(t in 1:n_cycles){
    for (x in 1:n_str){
   pest_cohort[[x]][t + 1, ] <- pest_cohort[[x]][t, ] %*% pest[[x]]
    }
  }
  
  ## Store the cohort traces in a list
  l_m_M <- pest_cohort
  names(l_m_M) <- v_names_str
  
  #### State Rewards scaled by the cycle length ####
  ## Vector of state utilities
  v_u_noban    <- c(H  = u_H, 
                    S1   = u_S1, 
                    S2   = u_S2,
                    D    = u_D) * cycle_length
  ## Vector of state costs
  v_c_nopest    <- c(H  = c_H_nopest, 
                     S1   = c_S1_nopest,
                     S2   = c_S2_nopest, 
                     D    = c_D) * cycle_length
  v_c_noban    <- c(H  = c_H, 
                    S1   = c_S1,
                    S2   = c_S2, 
                    D    = c_D) * cycle_length
  
  ## Store the vectors of state utilities for each strategy in a list 
  l_u <- lapply(seq_len(n_str), function(X) v_u_noban)
  
  ## Store the vectors of state cost for each strategy in a list 
  l_c <- lapply(seq_len(11), function(X) v_c_noban)
  l_c <- c(l_c,list(v_c_nopest))  
  # assign strategy names to matching items in the lists
  names(l_u) <- names(l_c) <- v_names_str
  
  ## create empty vectors to store total utilities and costs 
  v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
  names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
  
  #### Loop through each strategy and calculate total utilities and costs ####
  for (i in 1:n_str) {
    v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
    v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
    
    ###* Expected QALYs and costs per cycle 
    ##* Vector of QALYs and Costs
    #* Apply state rewards 
    v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
    v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
    
    ####* Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable
    #* QALYs
    v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
    #* Costs
    v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
  }
  df_cea <-as.data.frame(cbind(v_tot_cost, v_tot_qaly))
  df_cea_final[[j]]<-df_cea
  
  final_proportion=matrix(nrow=n_str, ncol=4)
  for (i in 1:n_str) {
    final_proportion[i,]=pest_cohort[[i]][21,]
  }
  colnames(final_proportion)<-colnames(pest_cohort[[1]])
  rownames(final_proportion) <- v_names_str
  df_cea_finalprop[[j]]<-final_proportion
}


df_cea_finalprop_1<-do.call(rbind.data.frame, df_cea_finalprop)
df_cea_finalprop_1$chemical<-seq(from = 1, to = n_str)
df_cea_finalprop_1 <- df_cea_finalprop_1[order(df_cea_finalprop_1$chemical),]
df_cea_finalprop_1$time<-seq(from = 1, to = 1000)

df_cea_finalprop_2<-reshape(df_cea_finalprop_1, idvar = "time", timevar = "chemical", direction = "wide")
df_cea_finalprop_2<-df_cea_finalprop_2[,-1]
mean_psaprop<-as.data.frame(colMeans(df_cea_finalprop_2))

H_final <- mean_psaprop[seq(1, nrow(mean_psaprop), 4), ]
S1_final <- mean_psaprop[seq(2, nrow(mean_psaprop), 4), ]
S2_final <- mean_psaprop[seq(3, nrow(mean_psaprop), 4), ]
D_final <- mean_psaprop[seq(4, nrow(mean_psaprop),4), ]


percentiles_975 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

H_25 <- percentiles_25[seq(1, nrow(percentiles_25), 4), ]
H_975 <- percentiles_975[seq(1, nrow(percentiles_975), 4), ]

S1_25 <- percentiles_25[seq(2, nrow(percentiles_25), 4), ]
S1_975 <- percentiles_975[seq(2, nrow(percentiles_975), 4), ]

S2_25 <- percentiles_25[seq(3, nrow(percentiles_25), 4), ]
S2_975 <- percentiles_975[seq(3, nrow(percentiles_975), 4), ]

D_25 <- percentiles_25[seq(4, nrow(percentiles_25), 4), ]
D_975 <- percentiles_975[seq(4, nrow(percentiles_975), 4), ]


data_prop_new <- cbind(S2_final, S2_25, S2_975,
                       D_final, D_25, D_975
)


data_prop_new<-as.data.frame(round(data_prop_new * 39440,digits=0))
data_prop_new_1<-data_prop_new
data_prop_new$S2_final<-data_prop_new$S2_final-1731
data_prop_new$S2_25<-data_prop_new$S2_25-1731
data_prop_new$S2_975<-data_prop_new$S2_975-1731

data_prop_new$D_final<-data_prop_new$D_final-22650
data_prop_new$D_25<-data_prop_new$D_25-22650
data_prop_new$D_975<-data_prop_new$D_975-22650

data_prop_new<-cbind(data_prop_new,v_names_str)
data_prop_new_1<-cbind(data_prop_new_1,v_names_str)


df_cea_final_1<-do.call(rbind.data.frame, df_cea_final)
df_cea_final_1$chemical<-seq(from = 1, to = n_str)
df_cea_final_1 <- df_cea_final_1[order(df_cea_final_1$chemical),]
df_cea_final_1$time<-seq(from = 1, to = 1000)

df_cea_final_2<-reshape(df_cea_final_1, idvar = "time", timevar = "chemical", direction = "wide")

df_cea_final_2<-df_cea_final_2[,-1]
mean_psa<-as.data.frame(colMeans(df_cea_final_2))

v_tot_qaly <- mean_psa[seq(2, nrow(mean_psa), 2), ]
v_tot_cost <- mean_psa[seq(1, nrow(mean_psa), 2), ]

percentiles_975 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

v_tot_qaly_25 <- percentiles_25[seq(2, nrow(percentiles_25), 2), ]
v_tot_qaly_975 <- percentiles_975[seq(2, nrow(percentiles_975), 2), ]

v_tot_cost_25 <- percentiles_25[seq(1, nrow(percentiles_25), 2), ]
v_tot_cost_975 <- percentiles_975[seq(1, nrow(percentiles_975), 2), ]


data_new <- cbind(v_tot_qaly, v_tot_qaly_25, v_tot_qaly_975, v_tot_cost,v_tot_cost_25,v_tot_cost_975,v_names_str)
chemname<-PD_onset[,c(1:2,5)]
data_new<-merge(data_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)
data_prop_new<-merge(data_prop_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)


data_new$'Inc_Cost'<-(as.numeric(data_new$v_tot_cost)-as.numeric(data_new$v_tot_cost[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_25'<-(as.numeric(data_new$v_tot_cost_25)-as.numeric(data_new$v_tot_cost_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_975'<-(as.numeric(data_new$v_tot_cost_975)-as.numeric(data_new$v_tot_cost_975[which(data_new$v_names_str=="no pest")]))

data_new$'Inc_Cost'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Cost')

data_new$'Inc_Effect'<-(as.numeric(data_new$v_tot_qaly)-as.numeric(data_new$v_tot_qaly[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_25'<-(as.numeric(data_new$v_tot_qaly_25)-as.numeric(data_new$v_tot_qaly_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_975'<-(as.numeric(data_new$v_tot_qaly_975)-as.numeric(data_new$v_tot_qaly_975[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Effect')

data_new$'ICER'<-data_new$'Inc_Cost'/data_new$'Inc_Effect'
data_new$'ICER_25'<-data_new$'Inc_Cost_25'/data_new$'Inc_Effect_25'
data_new$'ICER_975'<-data_new$'Inc_Cost_975'/data_new$'Inc_Effect_975'




#For single
df_cea_final_2_single<-df_cea_final_2
data_new_single<-data_new
data_prop_new_1_single<-data_prop_new_1
data_prop_new_single<-data_prop_new

#For double
for(j in 1:1000) {
  ## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
  r_S1S2 <- as.numeric(df_psa_1[j,])
  r_HS2  <- 0.0052
  r_HD   <- 0.0418
  r_S1D  <- 0.0418
  r_S2D  <- 0.0418
  ## State rewards
  # Costs
  c_H_nopest       <- 0.3*2      # annual cost of being Healthy unexposed
  c_S1_nopest      <- 0.3*2         # annual cost of being Healthy exposed
  c_S2_nopest      <- 0.3*2+as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  
  c_H       <- 0         # annual cost of being Healthy unexposed
  c_S1      <- 0         # annual cost of being Healthy exposed
  c_S2      <- as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  c_D       <- 0         # annual cost of being dead
  
  # Utilities
  u_H    <- 1           # annual utility of being Healthy unexposed
  u_S1   <- 1            # annual utility of being Healthy exposed
  u_S2   <- as.numeric(u_S2_psa[j,])     # annual utility of being PD Onset
  u_D    <- 0            # annual utility of being Dead
  
  # Discount weight for costs and effects
  v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  
  ####################### Construct state-transition models ######################
  ## Initial state vector
  # All starting healthy
  #California Population 65+ 5976000
  v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
  
  ## Initialize cohort trace for SoC
  m_M_NoPest <- matrix(NA,
                       nrow = (n_cycles + 1), ncol = n_states, 
                       dimnames = list(0:n_cycles, v_names_states))
  # Store the initial state vector in the first row of the cohort trace
  m_M_NoPest[1, ] <- v_m_init
  ## Initialize cohort trace for strategies A, B, and AB
  # Structure and initial states are the same as for SoC
  pest_cohort <- lapply(seq_len(n_str), function(X) m_M_NoPest)
  
  for (X in 1:n_str) {
    pest_cohort[[X]][1,"S1"]<- exposure_prob$prob[X]
    pest_cohort[[X]][1,"H"]<- 1-exposure_prob$prob[X]
    
  }
  
  ## Initialize transition probability matrix for Pesticide Ban
  # all transitions to a non-death state are assumed to be conditional on survival 
  m_P_NoPest <- matrix(0, 
                       nrow = n_states, ncol = n_states, 
                       dimnames = list(v_names_states, 
                                       v_names_states)) # define row and column names
  
  ## Initialize transition probability for Pesticide
  # Update transition probabilities
  pest <- lapply(seq_len(n_str), function(X) m_P_NoPest)
  
  ## Fill in matrix
  for (x in 1:n_str) {
    pest[[x]]["H", "H"]   <- 1 -r_HS2-r_HD
    pest[[x]]["H", "S2"]  <- r_HS2
    pest[[x]]["H", "D"]   <- r_HD
    # From S1
    pest[[x]]["S1", "S1"] <- 1 -r_S1S2[x]- r_S1D
    pest[[x]]["S1", "S2"] <- r_S1S2[x]
    pest[[x]]["S1", "D"]  <- r_S1D
    # From S2
    pest[[x]]["S2", "S2"] <- 1 - r_S2D
    pest[[x]]["S2", "D"]  <- r_S2D
    # From D
    pest[[x]]["D", "D"]   <- 1
  }
  
  ####################### Run Markov model #######################
  # Iterative solution of time-independent cSTM
  for(t in 1:n_cycles){
    for (x in 1:n_str){
      pest_cohort[[x]][t + 1, ] <- pest_cohort[[x]][t, ] %*% pest[[x]]
    }
  }
  
  ## Store the cohort traces in a list
  l_m_M <- pest_cohort
  names(l_m_M) <- v_names_str
  
  #### State Rewards scaled by the cycle length ####
  ## Vector of state utilities
  v_u_noban    <- c(H  = u_H, 
                    S1   = u_S1, 
                    S2   = u_S2,
                    D    = u_D) * cycle_length
  ## Vector of state costs
  v_c_nopest    <- c(H  = c_H_nopest, 
                     S1   = c_S1_nopest,
                     S2   = c_S2_nopest, 
                     D    = c_D) * cycle_length
  v_c_noban    <- c(H  = c_H, 
                    S1   = c_S1,
                    S2   = c_S2, 
                    D    = c_D) * cycle_length
  
  ## Store the vectors of state utilities for each strategy in a list 
  l_u <- lapply(seq_len(n_str), function(X) v_u_noban)
  
  ## Store the vectors of state cost for each strategy in a list 
  l_c <- lapply(seq_len(11), function(X) v_c_noban)
  l_c <- c(l_c,list(v_c_nopest))  
  # assign strategy names to matching items in the lists
  names(l_u) <- names(l_c) <- v_names_str
  
  ## create empty vectors to store total utilities and costs 
  v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
  names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
  
  #### Loop through each strategy and calculate total utilities and costs ####
  for (i in 1:n_str) {
    v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
    v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
    
    ###* Expected QALYs and costs per cycle 
    ##* Vector of QALYs and Costs
    #* Apply state rewards 
    v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
    v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
    
    ####* Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable
    #* QALYs
    v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
    #* Costs
    v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
  }
  df_cea <-as.data.frame(cbind(v_tot_cost, v_tot_qaly))
  df_cea_final[[j]]<-df_cea
  
  final_proportion=matrix(nrow=n_str, ncol=4)
  for (i in 1:n_str) {
    final_proportion[i,]=pest_cohort[[i]][21,]
  }
  colnames(final_proportion)<-colnames(pest_cohort[[1]])
  rownames(final_proportion) <- v_names_str
  df_cea_finalprop[[j]]<-final_proportion
}


df_cea_finalprop_1<-do.call(rbind.data.frame, df_cea_finalprop)
df_cea_finalprop_1$chemical<-seq(from = 1, to = n_str)
df_cea_finalprop_1 <- df_cea_finalprop_1[order(df_cea_finalprop_1$chemical),]
df_cea_finalprop_1$time<-seq(from = 1, to = 1000)

df_cea_finalprop_2<-reshape(df_cea_finalprop_1, idvar = "time", timevar = "chemical", direction = "wide")
df_cea_finalprop_2<-df_cea_finalprop_2[,-1]
mean_psaprop<-as.data.frame(colMeans(df_cea_finalprop_2))

H_final <- mean_psaprop[seq(1, nrow(mean_psaprop), 4), ]
S1_final <- mean_psaprop[seq(2, nrow(mean_psaprop), 4), ]
S2_final <- mean_psaprop[seq(3, nrow(mean_psaprop), 4), ]
D_final <- mean_psaprop[seq(4, nrow(mean_psaprop),4), ]


percentiles_975 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

H_25 <- percentiles_25[seq(1, nrow(percentiles_25), 4), ]
H_975 <- percentiles_975[seq(1, nrow(percentiles_975), 4), ]

S1_25 <- percentiles_25[seq(2, nrow(percentiles_25), 4), ]
S1_975 <- percentiles_975[seq(2, nrow(percentiles_975), 4), ]

S2_25 <- percentiles_25[seq(3, nrow(percentiles_25), 4), ]
S2_975 <- percentiles_975[seq(3, nrow(percentiles_975), 4), ]

D_25 <- percentiles_25[seq(4, nrow(percentiles_25), 4), ]
D_975 <- percentiles_975[seq(4, nrow(percentiles_975), 4), ]


data_prop_new <- cbind(S2_final, S2_25, S2_975,
                       D_final, D_25, D_975
)


data_prop_new<-as.data.frame(round(data_prop_new * 39440,digits=0))
data_prop_new_1<-data_prop_new
data_prop_new$S2_final<-data_prop_new$S2_final-1731
data_prop_new$S2_25<-data_prop_new$S2_25-1731
data_prop_new$S2_975<-data_prop_new$S2_975-1731

data_prop_new$D_final<-data_prop_new$D_final-22650
data_prop_new$D_25<-data_prop_new$D_25-22650
data_prop_new$D_975<-data_prop_new$D_975-22650

data_prop_new<-cbind(data_prop_new,v_names_str)
data_prop_new_1<-cbind(data_prop_new_1,v_names_str)


df_cea_final_1<-do.call(rbind.data.frame, df_cea_final)
df_cea_final_1$chemical<-seq(from = 1, to = n_str)
df_cea_final_1 <- df_cea_final_1[order(df_cea_final_1$chemical),]
df_cea_final_1$time<-seq(from = 1, to = 1000)

df_cea_final_2<-reshape(df_cea_final_1, idvar = "time", timevar = "chemical", direction = "wide")

df_cea_final_2<-df_cea_final_2[,-1]
mean_psa<-as.data.frame(colMeans(df_cea_final_2))

v_tot_qaly <- mean_psa[seq(2, nrow(mean_psa), 2), ]
v_tot_cost <- mean_psa[seq(1, nrow(mean_psa), 2), ]

percentiles_975 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

v_tot_qaly_25 <- percentiles_25[seq(2, nrow(percentiles_25), 2), ]
v_tot_qaly_975 <- percentiles_975[seq(2, nrow(percentiles_975), 2), ]

v_tot_cost_25 <- percentiles_25[seq(1, nrow(percentiles_25), 2), ]
v_tot_cost_975 <- percentiles_975[seq(1, nrow(percentiles_975), 2), ]


data_new <- cbind(v_tot_qaly, v_tot_qaly_25, v_tot_qaly_975, v_tot_cost,v_tot_cost_25,v_tot_cost_975,v_names_str)
chemname<-PD_onset[,c(1:2,5)]
data_new<-merge(data_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)
data_prop_new<-merge(data_prop_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)


data_new$'Inc_Cost'<-(as.numeric(data_new$v_tot_cost)-as.numeric(data_new$v_tot_cost[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_25'<-(as.numeric(data_new$v_tot_cost_25)-as.numeric(data_new$v_tot_cost_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_975'<-(as.numeric(data_new$v_tot_cost_975)-as.numeric(data_new$v_tot_cost_975[which(data_new$v_names_str=="no pest")]))

data_new$'Inc_Cost'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Cost')

data_new$'Inc_Effect'<-(as.numeric(data_new$v_tot_qaly)-as.numeric(data_new$v_tot_qaly[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_25'<-(as.numeric(data_new$v_tot_qaly_25)-as.numeric(data_new$v_tot_qaly_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_975'<-(as.numeric(data_new$v_tot_qaly_975)-as.numeric(data_new$v_tot_qaly_975[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Effect')

data_new$'ICER'<-data_new$'Inc_Cost'/data_new$'Inc_Effect'
data_new$'ICER_25'<-data_new$'Inc_Cost_25'/data_new$'Inc_Effect_25'
data_new$'ICER_975'<-data_new$'Inc_Cost_975'/data_new$'Inc_Effect_975'


df_cea_final_2_double<-df_cea_final_2
data_new_double<-data_new
data_prop_new_1_double<-data_prop_new_1
data_prop_new_double<-data_prop_new

#For triple

for(j in 1:1000) {
  ## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
  r_S1S2 <- as.numeric(df_psa_1[j,])
  r_HS2  <- 0.0052
  r_HD   <- 0.0418
  r_S1D  <- 0.0418
  r_S2D  <- 0.0418
  ## State rewards
  # Costs
  c_H_nopest       <- 0.3*3      # annual cost of being Healthy unexposed
  c_S1_nopest      <- 0.3*3         # annual cost of being Healthy exposed
  c_S2_nopest      <- 0.3*3+as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  
  c_H       <- 0         # annual cost of being Healthy unexposed
  c_S1      <- 0         # annual cost of being Healthy exposed
  c_S2      <- as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  c_D       <- 0         # annual cost of being dead
  
  # Utilities
  u_H    <- 1           # annual utility of being Healthy unexposed
  u_S1   <- 1            # annual utility of being Healthy exposed
  u_S2   <- as.numeric(u_S2_psa[j,])     # annual utility of being PD Onset
  u_D    <- 0            # annual utility of being Dead
  
  # Discount weight for costs and effects
  v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  
  ####################### Construct state-transition models ######################
  ## Initial state vector
  # All starting healthy
  #California Population 65+ 5976000
  v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
  
  ## Initialize cohort trace for SoC
  m_M_NoPest <- matrix(NA,
                       nrow = (n_cycles + 1), ncol = n_states, 
                       dimnames = list(0:n_cycles, v_names_states))
  # Store the initial state vector in the first row of the cohort trace
  m_M_NoPest[1, ] <- v_m_init
  ## Initialize cohort trace for strategies A, B, and AB
  # Structure and initial states are the same as for SoC
  pest_cohort <- lapply(seq_len(n_str), function(X) m_M_NoPest)
  
  for (X in 1:n_str) {
    pest_cohort[[X]][1,"S1"]<- exposure_prob$prob[X]
    pest_cohort[[X]][1,"H"]<- 1-exposure_prob$prob[X]
    
  }
  
  ## Initialize transition probability matrix for Pesticide Ban
  # all transitions to a non-death state are assumed to be conditional on survival 
  m_P_NoPest <- matrix(0, 
                       nrow = n_states, ncol = n_states, 
                       dimnames = list(v_names_states, 
                                       v_names_states)) # define row and column names
  
  ## Initialize transition probability for Pesticide
  # Update transition probabilities
  pest <- lapply(seq_len(n_str), function(X) m_P_NoPest)
  
  ## Fill in matrix
  for (x in 1:n_str) {
    pest[[x]]["H", "H"]   <- 1 -r_HS2-r_HD
    pest[[x]]["H", "S2"]  <- r_HS2
    pest[[x]]["H", "D"]   <- r_HD
    # From S1
    pest[[x]]["S1", "S1"] <- 1 -r_S1S2[x]- r_S1D
    pest[[x]]["S1", "S2"] <- r_S1S2[x]
    pest[[x]]["S1", "D"]  <- r_S1D
    # From S2
    pest[[x]]["S2", "S2"] <- 1 - r_S2D
    pest[[x]]["S2", "D"]  <- r_S2D
    # From D
    pest[[x]]["D", "D"]   <- 1
  }
  
  ####################### Run Markov model #######################
  # Iterative solution of time-independent cSTM
  for(t in 1:n_cycles){
    for (x in 1:n_str){
      pest_cohort[[x]][t + 1, ] <- pest_cohort[[x]][t, ] %*% pest[[x]]
    }
  }
  
  ## Store the cohort traces in a list
  l_m_M <- pest_cohort
  names(l_m_M) <- v_names_str
  
  #### State Rewards scaled by the cycle length ####
  ## Vector of state utilities
  v_u_noban    <- c(H  = u_H, 
                    S1   = u_S1, 
                    S2   = u_S2,
                    D    = u_D) * cycle_length
  ## Vector of state costs
  v_c_nopest    <- c(H  = c_H_nopest, 
                     S1   = c_S1_nopest,
                     S2   = c_S2_nopest, 
                     D    = c_D) * cycle_length
  v_c_noban    <- c(H  = c_H, 
                    S1   = c_S1,
                    S2   = c_S2, 
                    D    = c_D) * cycle_length
  
  ## Store the vectors of state utilities for each strategy in a list 
  l_u <- lapply(seq_len(n_str), function(X) v_u_noban)
  
  ## Store the vectors of state cost for each strategy in a list 
  l_c <- lapply(seq_len(11), function(X) v_c_noban)
  l_c <- c(l_c,list(v_c_nopest))  
  # assign strategy names to matching items in the lists
  names(l_u) <- names(l_c) <- v_names_str
  
  ## create empty vectors to store total utilities and costs 
  v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
  names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
  
  #### Loop through each strategy and calculate total utilities and costs ####
  for (i in 1:n_str) {
    v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
    v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
    
    ###* Expected QALYs and costs per cycle 
    ##* Vector of QALYs and Costs
    #* Apply state rewards 
    v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
    v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
    
    ####* Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable
    #* QALYs
    v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
    #* Costs
    v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
  }
  df_cea <-as.data.frame(cbind(v_tot_cost, v_tot_qaly))
  df_cea_final[[j]]<-df_cea
  
  final_proportion=matrix(nrow=n_str, ncol=4)
  for (i in 1:n_str) {
    final_proportion[i,]=pest_cohort[[i]][21,]
  }
  colnames(final_proportion)<-colnames(pest_cohort[[1]])
  rownames(final_proportion) <- v_names_str
  df_cea_finalprop[[j]]<-final_proportion
}


df_cea_finalprop_1<-do.call(rbind.data.frame, df_cea_finalprop)
df_cea_finalprop_1$chemical<-seq(from = 1, to = n_str)
df_cea_finalprop_1 <- df_cea_finalprop_1[order(df_cea_finalprop_1$chemical),]
df_cea_finalprop_1$time<-seq(from = 1, to = 1000)

df_cea_finalprop_2<-reshape(df_cea_finalprop_1, idvar = "time", timevar = "chemical", direction = "wide")
df_cea_finalprop_2<-df_cea_finalprop_2[,-1]
mean_psaprop<-as.data.frame(colMeans(df_cea_finalprop_2))

H_final <- mean_psaprop[seq(1, nrow(mean_psaprop), 4), ]
S1_final <- mean_psaprop[seq(2, nrow(mean_psaprop), 4), ]
S2_final <- mean_psaprop[seq(3, nrow(mean_psaprop), 4), ]
D_final <- mean_psaprop[seq(4, nrow(mean_psaprop),4), ]


percentiles_975 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

H_25 <- percentiles_25[seq(1, nrow(percentiles_25), 4), ]
H_975 <- percentiles_975[seq(1, nrow(percentiles_975), 4), ]

S1_25 <- percentiles_25[seq(2, nrow(percentiles_25), 4), ]
S1_975 <- percentiles_975[seq(2, nrow(percentiles_975), 4), ]

S2_25 <- percentiles_25[seq(3, nrow(percentiles_25), 4), ]
S2_975 <- percentiles_975[seq(3, nrow(percentiles_975), 4), ]

D_25 <- percentiles_25[seq(4, nrow(percentiles_25), 4), ]
D_975 <- percentiles_975[seq(4, nrow(percentiles_975), 4), ]


data_prop_new <- cbind(S2_final, S2_25, S2_975,
                       D_final, D_25, D_975
)


data_prop_new<-as.data.frame(round(data_prop_new * 39440,digits=0))
data_prop_new_1<-data_prop_new
data_prop_new$S2_final<-data_prop_new$S2_final-1731
data_prop_new$S2_25<-data_prop_new$S2_25-1731
data_prop_new$S2_975<-data_prop_new$S2_975-1731

data_prop_new$D_final<-data_prop_new$D_final-22650
data_prop_new$D_25<-data_prop_new$D_25-22650
data_prop_new$D_975<-data_prop_new$D_975-22650

data_prop_new<-cbind(data_prop_new,v_names_str)
data_prop_new_1<-cbind(data_prop_new_1,v_names_str)


df_cea_final_1<-do.call(rbind.data.frame, df_cea_final)
df_cea_final_1$chemical<-seq(from = 1, to = n_str)
df_cea_final_1 <- df_cea_final_1[order(df_cea_final_1$chemical),]
df_cea_final_1$time<-seq(from = 1, to = 1000)

df_cea_final_2<-reshape(df_cea_final_1, idvar = "time", timevar = "chemical", direction = "wide")

df_cea_final_2<-df_cea_final_2[,-1]
mean_psa<-as.data.frame(colMeans(df_cea_final_2))

v_tot_qaly <- mean_psa[seq(2, nrow(mean_psa), 2), ]
v_tot_cost <- mean_psa[seq(1, nrow(mean_psa), 2), ]

percentiles_975 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

v_tot_qaly_25 <- percentiles_25[seq(2, nrow(percentiles_25), 2), ]
v_tot_qaly_975 <- percentiles_975[seq(2, nrow(percentiles_975), 2), ]

v_tot_cost_25 <- percentiles_25[seq(1, nrow(percentiles_25), 2), ]
v_tot_cost_975 <- percentiles_975[seq(1, nrow(percentiles_975), 2), ]


data_new <- cbind(v_tot_qaly, v_tot_qaly_25, v_tot_qaly_975, v_tot_cost,v_tot_cost_25,v_tot_cost_975,v_names_str)
chemname<-PD_onset[,c(1:2,5)]
data_new<-merge(data_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)
data_prop_new<-merge(data_prop_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)


data_new$'Inc_Cost'<-(as.numeric(data_new$v_tot_cost)-as.numeric(data_new$v_tot_cost[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_25'<-(as.numeric(data_new$v_tot_cost_25)-as.numeric(data_new$v_tot_cost_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_975'<-(as.numeric(data_new$v_tot_cost_975)-as.numeric(data_new$v_tot_cost_975[which(data_new$v_names_str=="no pest")]))

data_new$'Inc_Cost'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Cost')

data_new$'Inc_Effect'<-(as.numeric(data_new$v_tot_qaly)-as.numeric(data_new$v_tot_qaly[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_25'<-(as.numeric(data_new$v_tot_qaly_25)-as.numeric(data_new$v_tot_qaly_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_975'<-(as.numeric(data_new$v_tot_qaly_975)-as.numeric(data_new$v_tot_qaly_975[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Effect')

data_new$'ICER'<-data_new$'Inc_Cost'/data_new$'Inc_Effect'
data_new$'ICER_25'<-data_new$'Inc_Cost_25'/data_new$'Inc_Effect_25'
data_new$'ICER_975'<-data_new$'Inc_Cost_975'/data_new$'Inc_Effect_975'


df_cea_final_2_triple<-df_cea_final_2
data_new_triple<-data_new
data_prop_new_1_triple<-data_prop_new_1
data_prop_new_triple<-data_prop_new



#Merge all
data_new_single[4,]<-data_new_double[4,]
data_new_single[8,]<-data_new_double[8,]
data_new_single[9,]<-data_new_double[9,]

data_new_single[10,]<-data_new_triple[10,]
data_new_single[13,]<-data_new_double[12,]
data_new_single[14,]<-data_new_triple[12,]

data_prop_new_1_single[4,]<-data_prop_new_1_double[4,]
data_prop_new_1_single[8,]<-data_prop_new_1_double[8,]
data_prop_new_1_single[9,]<-data_prop_new_1_double[9,]

data_prop_new_1_single[10,]<-data_prop_new_1_triple[10,]
data_prop_new_1_single[13,]<-data_prop_new_1_double[12,]
data_prop_new_1_single[14,]<-data_prop_new_1_triple[12,]

data_prop_new_single[4,]<-data_prop_new_double[4,]
data_prop_new_single[8,]<-data_prop_new_double[8,]
data_prop_new_single[9,]<-data_prop_new_double[9,]

data_prop_new_single[10,]<-data_prop_new_triple[10,]
data_prop_new_single[13,]<-data_prop_new_double[12,]
data_prop_new_single[14,]<-data_prop_new_triple[12,]



# write.csv(df_cea_final_2,"Manuscript/df_cea_residential_new_psa.csv", row.names = FALSE)
write.csv(data_new_single,"Manuscript/df_cea_residential_new_psa_mean.csv", row.names = FALSE)
write.csv(data_prop_new_1_single,"Manuscript/df_cea_residential_new_psa_prop_mean_1.csv", row.names = FALSE)
write.csv(data_prop_new_single,"Manuscript/df_cea_residential_new_psa_prop_mean.csv", row.names = FALSE)

















#Occupational
rm(list = ls())   

seed=071818
set.seed(seed)
#Get Exposure Probability
# occupational_pest<-read.csv("Data_Exposure/c_avg_alltime.csv")
occupational_pest <- read_sas("Data_Exposure/c_avg_lag10.sas7bdat")
occupational_pest<-occupational_pest[which(occupational_pest$PD==0),]

PD_onset<-read.csv("Data_Probability/Transition_Prob_PD.csv")
PD_onset<-PD_onset[PD_onset$Usage=="Insecticide",]
PD_onset<-PD_onset[PD_onset$Banning=="Restricted use",]

v_names_str <- c(PD_onset$Chemcode) 
occupational_pest<-occupational_pest[,which(colnames(occupational_pest) %in% v_names_str)]

exposure_prob<-as.data.frame(apply(occupational_pest,2,function(x) sum(x > 0)))
exposure_prob <- tibble::rownames_to_column(exposure_prob, "Chemcode")
exposure_prob$prob<-exposure_prob$`apply(occupational_pest, 2, function(x) sum(x > 0))`/nrow(occupational_pest)

both12<-sum(occupational_pest$Chem383 > 0 & occupational_pest$Chem216 > 0)
both13<-sum(occupational_pest$Chem383 > 0 & occupational_pest$Chem105 > 0)
both23<-sum(occupational_pest$Chem216 > 0 & occupational_pest$Chem105 > 0)
both123<-sum(occupational_pest$Chem383 > 0 & occupational_pest$Chem216 > 0 & occupational_pest$Chem105 > 0)

# sort by Chemcode
#add multiple exposure
exposure_prob[8,1]<-"Chem383*Chem216"
exposure_prob[8,2]<-both12
exposure_prob[8,3]<-both12/nrow(occupational_pest)

exposure_prob[9,1]<-"Chem383*Chem105"
exposure_prob[9,2]<-both13
exposure_prob[9,3]<-both13/nrow(occupational_pest)

exposure_prob[10,1]<-"Chem216*Chem105"
exposure_prob[10,2]<-both23
exposure_prob[10,3]<-both23/nrow(occupational_pest)

exposure_prob[11,1]<-"Chem383*Chem216*Chem105"
exposure_prob[11,2]<-both123
exposure_prob[11,3]<-both123/nrow(occupational_pest)

exposure_prob[12,1]<-"no pest"
exposure_prob[12,2]<-0
exposure_prob[12,3]<-0

exposure_prob <- exposure_prob[order(exposure_prob$Chemcode),]


occupational_para<-read.csv("Data_Probability/psa_occupational.csv")
# sort by Chemcode
PD_onset[8,1]<-"Chem383*Chem216"
PD_onset[8,2]<-"Oxydemeton-Methyl*Dimethoate"
PD_onset[8,3]<-1.198*1.138
PD_onset[8,4]<-1.21100*1.208
PD_onset[8,5]<-"Insecticide"

PD_onset[9,1]<-"Chem383*Chem105"
PD_onset[9,2]<-"Oxydemeton-Methyl*Carbaryl"
PD_onset[9,3]<-1.198*1.119
PD_onset[9,4]<-1.21100*1.168
PD_onset[9,5]<-"Insecticide"

PD_onset[10,1]<-"Chem216*Chem105"
PD_onset[10,2]<-"Dimethoate*Carbaryl"
PD_onset[10,3]<-1.138*1.119
PD_onset[10,4]<-1.208*1.168
PD_onset[10,5]<-"Insecticide"

PD_onset[11,1]<-"Chem383*Chem216*Chem105"
PD_onset[11,2]<-"Oxydemeton-Methyl*Dimethoate*Carbaryl"
PD_onset[11,3]<-1.198000*1.138*1.119
PD_onset[11,4]<-1.21100*1.208*1.168
PD_onset[11,5]<-"Insecticide"

PD_onset[12,1]<-"no pest"
PD_onset[12,2]<-"no pest"
PD_onset[12,3]<-1
PD_onset[12,4]<-1
PD_onset[12,5]<-"no pest"

PD_onset <- PD_onset[order(PD_onset$Chemcode),]
v_names_str <- c(PD_onset$Chemcode) 
occupational_para<-occupational_para[which(occupational_para$chemcode %in% v_names_str),]

occupational_para <- occupational_para[order(occupational_para$chemcode),]
occupational_para[12,1]<-"no pest"
occupational_para[12,2]<-0
occupational_para[12,3]<-0

## General setup
cycle_length <- 1       # cycle length equal to one year (use 1/12 for monthly)
n_age_init <- 65        # age at baseline
n_age_max  <- 85       # maximum age of follow up
n_cycles <- (n_age_max - n_age_init)/cycle_length # time horizon, number of cycles
v_names_states <- c("H",  # the 4 health states of the model:
                    "S1", # Healthy Unexposed (H), Healthy Exposed (S1), PD Onset (S2), Dead (D)
                    "S2",
                    "D")
n_states    <- length(v_names_states)     # number of health states 

# Discounting factors
d_c <- 0.03 # annual discount rate for costs 
d_e <- 0.03 # annual discount rate for QALYs

# Strategies
v_names_str <- c(PD_onset$Chemcode) 
n_str       <- length(v_names_str)        # number of strategies

# Within-cycle correction (WCC) using Simpson's 1/3 rule
v_wcc <- darthtools::gen_wcc(n_cycles = n_cycles, 
                             method = "Simpson1/3") # vector of wcc

n_sim=1000
# Given mean and standard deviation
mean_val <- 23041
sd_val <- 34045

# Calculate shape (k) and scale (θ) parameters
shape_param <- (mean_val / sd_val)^2
scale_param <- sd_val^2 / mean_val

# Generate random data from the gamma distribution
gamma_data <- rgamma(n = n_sim, shape = shape_param, scale = scale_param)

# Transform the data to achieve left-skewness
c_S2_psa <- abs(gamma_data - 2 * mean_val)
c_S2_psa=as.matrix(c_S2_psa)

# Utilities 
u_S2    = runif(n_sim, min = 0.61, max = 0.74)
u_S2_psa=as.matrix(u_S2)

r_HS2  <- 0.0052
df_psa_1 <- data.frame(matrix(ncol = n_str, nrow = 1000))
colnames(df_psa_1)<-v_names_str

for(m in 1:n_str) {
  r_S1S2  = rlnorm(n_sim, meanlog = occupational_para[m,2],  sdlog = occupational_para[m,3])*r_HS2
  df_psa_1[,m]<-r_S1S2
}

df_cea_final<-lapply(1:1000, matrix, data= NA, nrow=n_str, ncol=2)
df_cea_finalprop<-lapply(1:1000, matrix, data= NA, nrow=n_str, ncol=4)

for(j in 1:1000) {
  ## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
  r_S1S2 <- as.numeric(df_psa_1[j,])
  r_HS2  <- 0.0052
  r_HD   <- 0.0418
  r_S1D  <- 0.0418
  r_S2D  <- 0.0418
  ## State rewards
  # Costs
  c_H_nopest       <- 0.3      # annual cost of being Healthy unexposed
  c_S1_nopest      <- 0.3         # annual cost of being Healthy exposed
  c_S2_nopest      <- 0.3+as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  
  c_H       <- 0         # annual cost of being Healthy unexposed
  c_S1      <- 0         # annual cost of being Healthy exposed
  c_S2      <- as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  c_D       <- 0         # annual cost of being dead
  
  # Utilities
  u_H    <- 1           # annual utility of being Healthy unexposed
  u_S1   <- 1            # annual utility of being Healthy exposed
  u_S2   <- as.numeric(u_S2_psa[j,])     # annual utility of being PD Onset
  u_D    <- 0            # annual utility of being Dead
  
  # Discount weight for costs and effects
  v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  
  ####################### Construct state-transition models ######################
  ## Initial state vector
  # All starting healthy
  #California Population 65+ 5976000
  v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
  
  ## Initialize cohort trace for SoC
  m_M_NoPest <- matrix(NA,
                       nrow = (n_cycles + 1), ncol = n_states, 
                       dimnames = list(0:n_cycles, v_names_states))
  # Store the initial state vector in the first row of the cohort trace
  m_M_NoPest[1, ] <- v_m_init
  ## Initialize cohort trace for strategies A, B, and AB
  # Structure and initial states are the same as for SoC
  pest_cohort <- lapply(seq_len(n_str), function(X) m_M_NoPest)
  
  for (X in 1:n_str) {
    pest_cohort[[X]][1,"S1"]<- exposure_prob$prob[X]
    pest_cohort[[X]][1,"H"]<- 1-exposure_prob$prob[X]
    
  }
  
  ## Initialize transition probability matrix for Pesticide Ban
  # all transitions to a non-death state are assumed to be conditional on survival 
  m_P_NoPest <- matrix(0, 
                       nrow = n_states, ncol = n_states, 
                       dimnames = list(v_names_states, 
                                       v_names_states)) # define row and column names
  
  ## Initialize transition probability for Pesticide
  # Update transition probabilities
  pest <- lapply(seq_len(n_str), function(X) m_P_NoPest)
  
  ## Fill in matrix
  for (x in 1:n_str) {
    pest[[x]]["H", "H"]   <- 1 -r_HS2-r_HD
    pest[[x]]["H", "S2"]  <- r_HS2
    pest[[x]]["H", "D"]   <- r_HD
    # From S1
    pest[[x]]["S1", "S1"] <- 1 -r_S1S2[x]- r_S1D
    pest[[x]]["S1", "S2"] <- r_S1S2[x]
    pest[[x]]["S1", "D"]  <- r_S1D
    # From S2
    pest[[x]]["S2", "S2"] <- 1 - r_S2D
    pest[[x]]["S2", "D"]  <- r_S2D
    # From D
    pest[[x]]["D", "D"]   <- 1
  }
  
  ####################### Run Markov model #######################
  # Iterative solution of time-independent cSTM
  for(t in 1:n_cycles){
    for (x in 1:n_str){
      pest_cohort[[x]][t + 1, ] <- pest_cohort[[x]][t, ] %*% pest[[x]]
    }
  }
  
  ## Store the cohort traces in a list
  l_m_M <- pest_cohort
  names(l_m_M) <- v_names_str
  
  #### State Rewards scaled by the cycle length ####
  ## Vector of state utilities
  v_u_noban    <- c(H  = u_H, 
                    S1   = u_S1, 
                    S2   = u_S2,
                    D    = u_D) * cycle_length
  ## Vector of state costs
  v_c_nopest    <- c(H  = c_H_nopest, 
                     S1   = c_S1_nopest,
                     S2   = c_S2_nopest, 
                     D    = c_D) * cycle_length
  v_c_noban    <- c(H  = c_H, 
                    S1   = c_S1,
                    S2   = c_S2, 
                    D    = c_D) * cycle_length
  
  ## Store the vectors of state utilities for each strategy in a list 
  l_u <- lapply(seq_len(n_str), function(X) v_u_noban)
  
  ## Store the vectors of state cost for each strategy in a list 
  l_c <- lapply(seq_len(11), function(X) v_c_noban)
  l_c <- c(l_c,list(v_c_nopest))  
  # assign strategy names to matching items in the lists
  names(l_u) <- names(l_c) <- v_names_str
  
  ## create empty vectors to store total utilities and costs 
  v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
  names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
  
  #### Loop through each strategy and calculate total utilities and costs ####
  for (i in 1:n_str) {
    v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
    v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
    
    ###* Expected QALYs and costs per cycle 
    ##* Vector of QALYs and Costs
    #* Apply state rewards 
    v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
    v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
    
    ####* Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable
    #* QALYs
    v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
    #* Costs
    v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
  }
  df_cea <-as.data.frame(cbind(v_tot_cost, v_tot_qaly))
  df_cea_final[[j]]<-df_cea
  
  final_proportion=matrix(nrow=n_str, ncol=4)
  for (i in 1:n_str) {
    final_proportion[i,]=pest_cohort[[i]][21,]
  }
  colnames(final_proportion)<-colnames(pest_cohort[[1]])
  rownames(final_proportion) <- v_names_str
  df_cea_finalprop[[j]]<-final_proportion
}


df_cea_finalprop_1<-do.call(rbind.data.frame, df_cea_finalprop)
df_cea_finalprop_1$chemical<-seq(from = 1, to = n_str)
df_cea_finalprop_1 <- df_cea_finalprop_1[order(df_cea_finalprop_1$chemical),]
df_cea_finalprop_1$time<-seq(from = 1, to = 1000)

df_cea_finalprop_2<-reshape(df_cea_finalprop_1, idvar = "time", timevar = "chemical", direction = "wide")
df_cea_finalprop_2<-df_cea_finalprop_2[,-1]
mean_psaprop<-as.data.frame(colMeans(df_cea_finalprop_2))

H_final <- mean_psaprop[seq(1, nrow(mean_psaprop), 4), ]
S1_final <- mean_psaprop[seq(2, nrow(mean_psaprop), 4), ]
S2_final <- mean_psaprop[seq(3, nrow(mean_psaprop), 4), ]
D_final <- mean_psaprop[seq(4, nrow(mean_psaprop),4), ]


percentiles_975 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

H_25 <- percentiles_25[seq(1, nrow(percentiles_25), 4), ]
H_975 <- percentiles_975[seq(1, nrow(percentiles_975), 4), ]

S1_25 <- percentiles_25[seq(2, nrow(percentiles_25), 4), ]
S1_975 <- percentiles_975[seq(2, nrow(percentiles_975), 4), ]

S2_25 <- percentiles_25[seq(3, nrow(percentiles_25), 4), ]
S2_975 <- percentiles_975[seq(3, nrow(percentiles_975), 4), ]

D_25 <- percentiles_25[seq(4, nrow(percentiles_25), 4), ]
D_975 <- percentiles_975[seq(4, nrow(percentiles_975), 4), ]


data_prop_new <- cbind(S2_final, S2_25, S2_975,
                       D_final, D_25, D_975
)


data_prop_new<-as.data.frame(round(data_prop_new * 39440,digits=0))
data_prop_new_1<-data_prop_new
data_prop_new$S2_final<-data_prop_new$S2_final-1731
data_prop_new$S2_25<-data_prop_new$S2_25-1731
data_prop_new$S2_975<-data_prop_new$S2_975-1731

data_prop_new$D_final<-data_prop_new$D_final-22650
data_prop_new$D_25<-data_prop_new$D_25-22650
data_prop_new$D_975<-data_prop_new$D_975-22650

data_prop_new<-cbind(data_prop_new,v_names_str)
data_prop_new_1<-cbind(data_prop_new_1,v_names_str)


df_cea_final_1<-do.call(rbind.data.frame, df_cea_final)
df_cea_final_1$chemical<-seq(from = 1, to = n_str)
df_cea_final_1 <- df_cea_final_1[order(df_cea_final_1$chemical),]
df_cea_final_1$time<-seq(from = 1, to = 1000)

df_cea_final_2<-reshape(df_cea_final_1, idvar = "time", timevar = "chemical", direction = "wide")

df_cea_final_2<-df_cea_final_2[,-1]
mean_psa<-as.data.frame(colMeans(df_cea_final_2))

v_tot_qaly <- mean_psa[seq(2, nrow(mean_psa), 2), ]
v_tot_cost <- mean_psa[seq(1, nrow(mean_psa), 2), ]

percentiles_975 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

v_tot_qaly_25 <- percentiles_25[seq(2, nrow(percentiles_25), 2), ]
v_tot_qaly_975 <- percentiles_975[seq(2, nrow(percentiles_975), 2), ]

v_tot_cost_25 <- percentiles_25[seq(1, nrow(percentiles_25), 2), ]
v_tot_cost_975 <- percentiles_975[seq(1, nrow(percentiles_975), 2), ]


data_new <- cbind(v_tot_qaly, v_tot_qaly_25, v_tot_qaly_975, v_tot_cost,v_tot_cost_25,v_tot_cost_975,v_names_str)
chemname<-PD_onset[,c(1:2,5)]
data_new<-merge(data_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)
data_prop_new<-merge(data_prop_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)


data_new$'Inc_Cost'<-(as.numeric(data_new$v_tot_cost)-as.numeric(data_new$v_tot_cost[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_25'<-(as.numeric(data_new$v_tot_cost_25)-as.numeric(data_new$v_tot_cost_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_975'<-(as.numeric(data_new$v_tot_cost_975)-as.numeric(data_new$v_tot_cost_975[which(data_new$v_names_str=="no pest")]))

data_new$'Inc_Cost'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Cost')

data_new$'Inc_Effect'<-(as.numeric(data_new$v_tot_qaly)-as.numeric(data_new$v_tot_qaly[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_25'<-(as.numeric(data_new$v_tot_qaly_25)-as.numeric(data_new$v_tot_qaly_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_975'<-(as.numeric(data_new$v_tot_qaly_975)-as.numeric(data_new$v_tot_qaly_975[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Effect')

data_new$'ICER'<-data_new$'Inc_Cost'/data_new$'Inc_Effect'
data_new$'ICER_25'<-data_new$'Inc_Cost_25'/data_new$'Inc_Effect_25'
data_new$'ICER_975'<-data_new$'Inc_Cost_975'/data_new$'Inc_Effect_975'

#For single
df_cea_final_2_single<-df_cea_final_2
data_new_single<-data_new
data_prop_new_1_single<-data_prop_new_1
data_prop_new_single<-data_prop_new



#For double

for(j in 1:1000) {
  ## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
  r_S1S2 <- as.numeric(df_psa_1[j,])
  r_HS2  <- 0.0052
  r_HD   <- 0.0418
  r_S1D  <- 0.0418
  r_S2D  <- 0.0418
  ## State rewards
  # Costs
  c_H_nopest       <- 0.3*2      # annual cost of being Healthy unexposed
  c_S1_nopest      <- 0.3*2         # annual cost of being Healthy exposed
  c_S2_nopest      <- 0.3*2+as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  
  c_H       <- 0         # annual cost of being Healthy unexposed
  c_S1      <- 0         # annual cost of being Healthy exposed
  c_S2      <- as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  c_D       <- 0         # annual cost of being dead
  
  # Utilities
  u_H    <- 1           # annual utility of being Healthy unexposed
  u_S1   <- 1            # annual utility of being Healthy exposed
  u_S2   <- as.numeric(u_S2_psa[j,])     # annual utility of being PD Onset
  u_D    <- 0            # annual utility of being Dead
  
  # Discount weight for costs and effects
  v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  
  ####################### Construct state-transition models ######################
  ## Initial state vector
  # All starting healthy
  #California Population 65+ 5976000
  v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
  
  ## Initialize cohort trace for SoC
  m_M_NoPest <- matrix(NA,
                       nrow = (n_cycles + 1), ncol = n_states, 
                       dimnames = list(0:n_cycles, v_names_states))
  # Store the initial state vector in the first row of the cohort trace
  m_M_NoPest[1, ] <- v_m_init
  ## Initialize cohort trace for strategies A, B, and AB
  # Structure and initial states are the same as for SoC
  pest_cohort <- lapply(seq_len(n_str), function(X) m_M_NoPest)
  
  for (X in 1:n_str) {
    pest_cohort[[X]][1,"S1"]<- exposure_prob$prob[X]
    pest_cohort[[X]][1,"H"]<- 1-exposure_prob$prob[X]
    
  }
  
  ## Initialize transition probability matrix for Pesticide Ban
  # all transitions to a non-death state are assumed to be conditional on survival 
  m_P_NoPest <- matrix(0, 
                       nrow = n_states, ncol = n_states, 
                       dimnames = list(v_names_states, 
                                       v_names_states)) # define row and column names
  
  ## Initialize transition probability for Pesticide
  # Update transition probabilities
  pest <- lapply(seq_len(n_str), function(X) m_P_NoPest)
  
  ## Fill in matrix
  for (x in 1:n_str) {
    pest[[x]]["H", "H"]   <- 1 -r_HS2-r_HD
    pest[[x]]["H", "S2"]  <- r_HS2
    pest[[x]]["H", "D"]   <- r_HD
    # From S1
    pest[[x]]["S1", "S1"] <- 1 -r_S1S2[x]- r_S1D
    pest[[x]]["S1", "S2"] <- r_S1S2[x]
    pest[[x]]["S1", "D"]  <- r_S1D
    # From S2
    pest[[x]]["S2", "S2"] <- 1 - r_S2D
    pest[[x]]["S2", "D"]  <- r_S2D
    # From D
    pest[[x]]["D", "D"]   <- 1
  }
  
  ####################### Run Markov model #######################
  # Iterative solution of time-independent cSTM
  for(t in 1:n_cycles){
    for (x in 1:n_str){
      pest_cohort[[x]][t + 1, ] <- pest_cohort[[x]][t, ] %*% pest[[x]]
    }
  }
  
  ## Store the cohort traces in a list
  l_m_M <- pest_cohort
  names(l_m_M) <- v_names_str
  
  #### State Rewards scaled by the cycle length ####
  ## Vector of state utilities
  v_u_noban    <- c(H  = u_H, 
                    S1   = u_S1, 
                    S2   = u_S2,
                    D    = u_D) * cycle_length
  ## Vector of state costs
  v_c_nopest    <- c(H  = c_H_nopest, 
                     S1   = c_S1_nopest,
                     S2   = c_S2_nopest, 
                     D    = c_D) * cycle_length
  v_c_noban    <- c(H  = c_H, 
                    S1   = c_S1,
                    S2   = c_S2, 
                    D    = c_D) * cycle_length
  
  ## Store the vectors of state utilities for each strategy in a list 
  l_u <- lapply(seq_len(n_str), function(X) v_u_noban)
  
  ## Store the vectors of state cost for each strategy in a list 
  l_c <- lapply(seq_len(11), function(X) v_c_noban)
  l_c <- c(l_c,list(v_c_nopest))  
  # assign strategy names to matching items in the lists
  names(l_u) <- names(l_c) <- v_names_str
  
  ## create empty vectors to store total utilities and costs 
  v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
  names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
  
  #### Loop through each strategy and calculate total utilities and costs ####
  for (i in 1:n_str) {
    v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
    v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
    
    ###* Expected QALYs and costs per cycle 
    ##* Vector of QALYs and Costs
    #* Apply state rewards 
    v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
    v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
    
    ####* Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable
    #* QALYs
    v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
    #* Costs
    v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
  }
  df_cea <-as.data.frame(cbind(v_tot_cost, v_tot_qaly))
  df_cea_final[[j]]<-df_cea
  
  final_proportion=matrix(nrow=n_str, ncol=4)
  for (i in 1:n_str) {
    final_proportion[i,]=pest_cohort[[i]][21,]
  }
  colnames(final_proportion)<-colnames(pest_cohort[[1]])
  rownames(final_proportion) <- v_names_str
  df_cea_finalprop[[j]]<-final_proportion
}


df_cea_finalprop_1<-do.call(rbind.data.frame, df_cea_finalprop)
df_cea_finalprop_1$chemical<-seq(from = 1, to = n_str)
df_cea_finalprop_1 <- df_cea_finalprop_1[order(df_cea_finalprop_1$chemical),]
df_cea_finalprop_1$time<-seq(from = 1, to = 1000)

df_cea_finalprop_2<-reshape(df_cea_finalprop_1, idvar = "time", timevar = "chemical", direction = "wide")
df_cea_finalprop_2<-df_cea_finalprop_2[,-1]
mean_psaprop<-as.data.frame(colMeans(df_cea_finalprop_2))

H_final <- mean_psaprop[seq(1, nrow(mean_psaprop), 4), ]
S1_final <- mean_psaprop[seq(2, nrow(mean_psaprop), 4), ]
S2_final <- mean_psaprop[seq(3, nrow(mean_psaprop), 4), ]
D_final <- mean_psaprop[seq(4, nrow(mean_psaprop),4), ]


percentiles_975 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

H_25 <- percentiles_25[seq(1, nrow(percentiles_25), 4), ]
H_975 <- percentiles_975[seq(1, nrow(percentiles_975), 4), ]

S1_25 <- percentiles_25[seq(2, nrow(percentiles_25), 4), ]
S1_975 <- percentiles_975[seq(2, nrow(percentiles_975), 4), ]

S2_25 <- percentiles_25[seq(3, nrow(percentiles_25), 4), ]
S2_975 <- percentiles_975[seq(3, nrow(percentiles_975), 4), ]

D_25 <- percentiles_25[seq(4, nrow(percentiles_25), 4), ]
D_975 <- percentiles_975[seq(4, nrow(percentiles_975), 4), ]


data_prop_new <- cbind(S2_final, S2_25, S2_975,
                       D_final, D_25, D_975
)


data_prop_new<-as.data.frame(round(data_prop_new * 39440,digits=0))
data_prop_new_1<-data_prop_new
data_prop_new$S2_final<-data_prop_new$S2_final-1731
data_prop_new$S2_25<-data_prop_new$S2_25-1731
data_prop_new$S2_975<-data_prop_new$S2_975-1731

data_prop_new$D_final<-data_prop_new$D_final-22650
data_prop_new$D_25<-data_prop_new$D_25-22650
data_prop_new$D_975<-data_prop_new$D_975-22650

data_prop_new<-cbind(data_prop_new,v_names_str)
data_prop_new_1<-cbind(data_prop_new_1,v_names_str)


df_cea_final_1<-do.call(rbind.data.frame, df_cea_final)
df_cea_final_1$chemical<-seq(from = 1, to = n_str)
df_cea_final_1 <- df_cea_final_1[order(df_cea_final_1$chemical),]
df_cea_final_1$time<-seq(from = 1, to = 1000)

df_cea_final_2<-reshape(df_cea_final_1, idvar = "time", timevar = "chemical", direction = "wide")

df_cea_final_2<-df_cea_final_2[,-1]
mean_psa<-as.data.frame(colMeans(df_cea_final_2))

v_tot_qaly <- mean_psa[seq(2, nrow(mean_psa), 2), ]
v_tot_cost <- mean_psa[seq(1, nrow(mean_psa), 2), ]

percentiles_975 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

v_tot_qaly_25 <- percentiles_25[seq(2, nrow(percentiles_25), 2), ]
v_tot_qaly_975 <- percentiles_975[seq(2, nrow(percentiles_975), 2), ]

v_tot_cost_25 <- percentiles_25[seq(1, nrow(percentiles_25), 2), ]
v_tot_cost_975 <- percentiles_975[seq(1, nrow(percentiles_975), 2), ]


data_new <- cbind(v_tot_qaly, v_tot_qaly_25, v_tot_qaly_975, v_tot_cost,v_tot_cost_25,v_tot_cost_975,v_names_str)
chemname<-PD_onset[,c(1:2,5)]
data_new<-merge(data_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)
data_prop_new<-merge(data_prop_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)


data_new$'Inc_Cost'<-(as.numeric(data_new$v_tot_cost)-as.numeric(data_new$v_tot_cost[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_25'<-(as.numeric(data_new$v_tot_cost_25)-as.numeric(data_new$v_tot_cost_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_975'<-(as.numeric(data_new$v_tot_cost_975)-as.numeric(data_new$v_tot_cost_975[which(data_new$v_names_str=="no pest")]))

data_new$'Inc_Cost'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Cost')

data_new$'Inc_Effect'<-(as.numeric(data_new$v_tot_qaly)-as.numeric(data_new$v_tot_qaly[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_25'<-(as.numeric(data_new$v_tot_qaly_25)-as.numeric(data_new$v_tot_qaly_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_975'<-(as.numeric(data_new$v_tot_qaly_975)-as.numeric(data_new$v_tot_qaly_975[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Effect')

data_new$'ICER'<-data_new$'Inc_Cost'/data_new$'Inc_Effect'
data_new$'ICER_25'<-data_new$'Inc_Cost_25'/data_new$'Inc_Effect_25'
data_new$'ICER_975'<-data_new$'Inc_Cost_975'/data_new$'Inc_Effect_975'

df_cea_final_2_double<-df_cea_final_2
data_new_double<-data_new
data_prop_new_1_double<-data_prop_new_1
data_prop_new_double<-data_prop_new


#For triple

for(j in 1:1000) {
  ## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
  r_S1S2 <- as.numeric(df_psa_1[j,])
  r_HS2  <- 0.0052
  r_HD   <- 0.0418
  r_S1D  <- 0.0418
  r_S2D  <- 0.0418
  ## State rewards
  # Costs
  c_H_nopest       <- 0.3*3      # annual cost of being Healthy unexposed
  c_S1_nopest      <- 0.3*3         # annual cost of being Healthy exposed
  c_S2_nopest      <- 0.3*3+as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  
  c_H       <- 0         # annual cost of being Healthy unexposed
  c_S1      <- 0         # annual cost of being Healthy exposed
  c_S2      <- as.numeric(c_S2_psa[j,])     # annual cost of being PD Onset
  c_D       <- 0         # annual cost of being dead
  
  # Utilities
  u_H    <- 1           # annual utility of being Healthy unexposed
  u_S1   <- 1            # annual utility of being Healthy exposed
  u_S2   <- as.numeric(u_S2_psa[j,])     # annual utility of being PD Onset
  u_D    <- 0            # annual utility of being Dead
  
  # Discount weight for costs and effects
  v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  
  ####################### Construct state-transition models ######################
  ## Initial state vector
  # All starting healthy
  #California Population 65+ 5976000
  v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
  
  ## Initialize cohort trace for SoC
  m_M_NoPest <- matrix(NA,
                       nrow = (n_cycles + 1), ncol = n_states, 
                       dimnames = list(0:n_cycles, v_names_states))
  # Store the initial state vector in the first row of the cohort trace
  m_M_NoPest[1, ] <- v_m_init
  ## Initialize cohort trace for strategies A, B, and AB
  # Structure and initial states are the same as for SoC
  pest_cohort <- lapply(seq_len(n_str), function(X) m_M_NoPest)
  
  for (X in 1:n_str) {
    pest_cohort[[X]][1,"S1"]<- exposure_prob$prob[X]
    pest_cohort[[X]][1,"H"]<- 1-exposure_prob$prob[X]
    
  }
  
  ## Initialize transition probability matrix for Pesticide Ban
  # all transitions to a non-death state are assumed to be conditional on survival 
  m_P_NoPest <- matrix(0, 
                       nrow = n_states, ncol = n_states, 
                       dimnames = list(v_names_states, 
                                       v_names_states)) # define row and column names
  
  ## Initialize transition probability for Pesticide
  # Update transition probabilities
  pest <- lapply(seq_len(n_str), function(X) m_P_NoPest)
  
  ## Fill in matrix
  for (x in 1:n_str) {
    pest[[x]]["H", "H"]   <- 1 -r_HS2-r_HD
    pest[[x]]["H", "S2"]  <- r_HS2
    pest[[x]]["H", "D"]   <- r_HD
    # From S1
    pest[[x]]["S1", "S1"] <- 1 -r_S1S2[x]- r_S1D
    pest[[x]]["S1", "S2"] <- r_S1S2[x]
    pest[[x]]["S1", "D"]  <- r_S1D
    # From S2
    pest[[x]]["S2", "S2"] <- 1 - r_S2D
    pest[[x]]["S2", "D"]  <- r_S2D
    # From D
    pest[[x]]["D", "D"]   <- 1
  }
  
  ####################### Run Markov model #######################
  # Iterative solution of time-independent cSTM
  for(t in 1:n_cycles){
    for (x in 1:n_str){
      pest_cohort[[x]][t + 1, ] <- pest_cohort[[x]][t, ] %*% pest[[x]]
    }
  }
  
  ## Store the cohort traces in a list
  l_m_M <- pest_cohort
  names(l_m_M) <- v_names_str
  
  #### State Rewards scaled by the cycle length ####
  ## Vector of state utilities
  v_u_noban    <- c(H  = u_H, 
                    S1   = u_S1, 
                    S2   = u_S2,
                    D    = u_D) * cycle_length
  ## Vector of state costs
  v_c_nopest    <- c(H  = c_H_nopest, 
                     S1   = c_S1_nopest,
                     S2   = c_S2_nopest, 
                     D    = c_D) * cycle_length
  v_c_noban    <- c(H  = c_H, 
                    S1   = c_S1,
                    S2   = c_S2, 
                    D    = c_D) * cycle_length
  
  ## Store the vectors of state utilities for each strategy in a list 
  l_u <- lapply(seq_len(n_str), function(X) v_u_noban)
  
  ## Store the vectors of state cost for each strategy in a list 
  l_c <- lapply(seq_len(11), function(X) v_c_noban)
  l_c <- c(l_c,list(v_c_nopest))  
  # assign strategy names to matching items in the lists
  names(l_u) <- names(l_c) <- v_names_str
  
  ## create empty vectors to store total utilities and costs 
  v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = n_str)
  names(v_tot_qaly) <- names(v_tot_cost) <- v_names_str
  
  #### Loop through each strategy and calculate total utilities and costs ####
  for (i in 1:n_str) {
    v_u_str <- l_u[[i]]   # select the vector of state utilities for the i-th strategy
    v_c_str <- l_c[[i]]   # select the vector of state costs for the i-th strategy
    
    ###* Expected QALYs and costs per cycle 
    ##* Vector of QALYs and Costs
    #* Apply state rewards 
    v_qaly_str <- l_m_M[[i]] %*% v_u_str # sum the utilities of all states for each cycle
    v_cost_str <- l_m_M[[i]] %*% v_c_str # sum the costs of all states for each cycle
    
    ####* Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable
    #* QALYs
    v_tot_qaly[i] <- t(v_qaly_str) %*% (v_dwe * v_wcc)
    #* Costs
    v_tot_cost[i] <- t(v_cost_str) %*% (v_dwc * v_wcc)
  }
  df_cea <-as.data.frame(cbind(v_tot_cost, v_tot_qaly))
  df_cea_final[[j]]<-df_cea
  
  final_proportion=matrix(nrow=n_str, ncol=4)
  for (i in 1:n_str) {
    final_proportion[i,]=pest_cohort[[i]][21,]
  }
  colnames(final_proportion)<-colnames(pest_cohort[[1]])
  rownames(final_proportion) <- v_names_str
  df_cea_finalprop[[j]]<-final_proportion
}


df_cea_finalprop_1<-do.call(rbind.data.frame, df_cea_finalprop)
df_cea_finalprop_1$chemical<-seq(from = 1, to = n_str)
df_cea_finalprop_1 <- df_cea_finalprop_1[order(df_cea_finalprop_1$chemical),]
df_cea_finalprop_1$time<-seq(from = 1, to = 1000)

df_cea_finalprop_2<-reshape(df_cea_finalprop_1, idvar = "time", timevar = "chemical", direction = "wide")
df_cea_finalprop_2<-df_cea_finalprop_2[,-1]
mean_psaprop<-as.data.frame(colMeans(df_cea_finalprop_2))

H_final <- mean_psaprop[seq(1, nrow(mean_psaprop), 4), ]
S1_final <- mean_psaprop[seq(2, nrow(mean_psaprop), 4), ]
S2_final <- mean_psaprop[seq(3, nrow(mean_psaprop), 4), ]
D_final <- mean_psaprop[seq(4, nrow(mean_psaprop),4), ]


percentiles_975 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_finalprop_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

H_25 <- percentiles_25[seq(1, nrow(percentiles_25), 4), ]
H_975 <- percentiles_975[seq(1, nrow(percentiles_975), 4), ]

S1_25 <- percentiles_25[seq(2, nrow(percentiles_25), 4), ]
S1_975 <- percentiles_975[seq(2, nrow(percentiles_975), 4), ]

S2_25 <- percentiles_25[seq(3, nrow(percentiles_25), 4), ]
S2_975 <- percentiles_975[seq(3, nrow(percentiles_975), 4), ]

D_25 <- percentiles_25[seq(4, nrow(percentiles_25), 4), ]
D_975 <- percentiles_975[seq(4, nrow(percentiles_975), 4), ]


data_prop_new <- cbind(S2_final, S2_25, S2_975,
                       D_final, D_25, D_975
)


data_prop_new<-as.data.frame(round(data_prop_new * 39440,digits=0))
data_prop_new_1<-data_prop_new
data_prop_new$S2_final<-data_prop_new$S2_final-1731
data_prop_new$S2_25<-data_prop_new$S2_25-1731
data_prop_new$S2_975<-data_prop_new$S2_975-1731

data_prop_new$D_final<-data_prop_new$D_final-22650
data_prop_new$D_25<-data_prop_new$D_25-22650
data_prop_new$D_975<-data_prop_new$D_975-22650

data_prop_new<-cbind(data_prop_new,v_names_str)
data_prop_new_1<-cbind(data_prop_new_1,v_names_str)


df_cea_final_1<-do.call(rbind.data.frame, df_cea_final)
df_cea_final_1$chemical<-seq(from = 1, to = n_str)
df_cea_final_1 <- df_cea_final_1[order(df_cea_final_1$chemical),]
df_cea_final_1$time<-seq(from = 1, to = 1000)

df_cea_final_2<-reshape(df_cea_final_1, idvar = "time", timevar = "chemical", direction = "wide")

df_cea_final_2<-df_cea_final_2[,-1]
mean_psa<-as.data.frame(colMeans(df_cea_final_2))

v_tot_qaly <- mean_psa[seq(2, nrow(mean_psa), 2), ]
v_tot_cost <- mean_psa[seq(1, nrow(mean_psa), 2), ]

percentiles_975 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.975))
percentiles_975<-as.data.frame(percentiles_975)

percentiles_25 <- apply(df_cea_final_2, 2, function(col) quantile(col, probs = 0.025))
percentiles_25<-as.data.frame(percentiles_25)

v_tot_qaly_25 <- percentiles_25[seq(2, nrow(percentiles_25), 2), ]
v_tot_qaly_975 <- percentiles_975[seq(2, nrow(percentiles_975), 2), ]

v_tot_cost_25 <- percentiles_25[seq(1, nrow(percentiles_25), 2), ]
v_tot_cost_975 <- percentiles_975[seq(1, nrow(percentiles_975), 2), ]


data_new <- cbind(v_tot_qaly, v_tot_qaly_25, v_tot_qaly_975, v_tot_cost,v_tot_cost_25,v_tot_cost_975,v_names_str)
chemname<-PD_onset[,c(1:2,5)]
data_new<-merge(data_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)
data_prop_new<-merge(data_prop_new,chemname,by.x="v_names_str",by.y="Chemcode",all.x=TRUE)


data_new$'Inc_Cost'<-(as.numeric(data_new$v_tot_cost)-as.numeric(data_new$v_tot_cost[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_25'<-(as.numeric(data_new$v_tot_cost_25)-as.numeric(data_new$v_tot_cost_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Cost_975'<-(as.numeric(data_new$v_tot_cost_975)-as.numeric(data_new$v_tot_cost_975[which(data_new$v_names_str=="no pest")]))

data_new$'Inc_Cost'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Cost')

data_new$'Inc_Effect'<-(as.numeric(data_new$v_tot_qaly)-as.numeric(data_new$v_tot_qaly[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_25'<-(as.numeric(data_new$v_tot_qaly_25)-as.numeric(data_new$v_tot_qaly_25[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect_975'<-(as.numeric(data_new$v_tot_qaly_975)-as.numeric(data_new$v_tot_qaly_975[which(data_new$v_names_str=="no pest")]))
data_new$'Inc_Effect'<-ifelse(data_new$'v_names_str'=="no pest",NA,data_new$'Inc_Effect')

data_new$'ICER'<-data_new$'Inc_Cost'/data_new$'Inc_Effect'
data_new$'ICER_25'<-data_new$'Inc_Cost_25'/data_new$'Inc_Effect_25'
data_new$'ICER_975'<-data_new$'Inc_Cost_975'/data_new$'Inc_Effect_975'

df_cea_final_2_triple<-df_cea_final_2
data_new_triple<-data_new
data_prop_new_1_triple<-data_prop_new_1
data_prop_new_triple<-data_prop_new

#Merge all
data_new_single[4,]<-data_new_double[4,]
data_new_single[8,]<-data_new_double[8,]
data_new_single[9,]<-data_new_double[9,]

data_new_single[10,]<-data_new_triple[10,]
data_new_single[13,]<-data_new_double[12,]
data_new_single[14,]<-data_new_triple[12,]


data_prop_new_1_single[4,]<-data_prop_new_1_double[4,]
data_prop_new_1_single[8,]<-data_prop_new_1_double[8,]
data_prop_new_1_single[9,]<-data_prop_new_1_double[9,]

data_prop_new_1_single[10,]<-data_prop_new_1_triple[10,]
data_prop_new_1_single[13,]<-data_prop_new_1_double[12,]
data_prop_new_1_single[14,]<-data_prop_new_1_triple[12,]



data_prop_new_single[4,]<-data_prop_new_double[4,]
data_prop_new_single[8,]<-data_prop_new_double[8,]
data_prop_new_single[9,]<-data_prop_new_double[9,]

data_prop_new_single[10,]<-data_prop_new_triple[10,]
data_prop_new_single[13,]<-data_prop_new_double[12,]
data_prop_new_single[14,]<-data_prop_new_triple[12,]



# write.csv(df_cea_final_2,"Manuscript/df_cea_residential_new_psa.csv", row.names = FALSE)
write.csv(data_new_single,"Manuscript/df_cea_occupational_new_psa_mean.csv", row.names = FALSE)
write.csv(data_prop_new_1_single,"Manuscript/df_cea_occupational_new_psa_prop_mean_1.csv", row.names = FALSE)
write.csv(data_prop_new_single,"Manuscript/df_cea_occupational_new_psa_prop_mean.csv", row.names = FALSE)




