### Load packages
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
library(tibble)
# devtools::install_github("DARTH-git/darthtools")
library(ggpubr)
library(readxl)
options(scipen = 999)
library(scales)
library(haven)

################################# Model input ##################################
#Residential
rm(list = ls())    # remove any variables in R's memory 

set.seed(123234)
#Get Exposure Probability
# residential_pest<-read.csv("Data_Exposure/r_avg_alltime.csv")
residential_pest <- read_sas("Data_Exposure/r_avg_lag10.sas7bdat")
residential_pest<-residential_pest[which(residential_pest$PD==0),]
mean(residential_pest$IndexYR10-1974)

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

exposure_prob <- exposure_prob[order(exposure_prob$Chemcode),]
PD_onset <- PD_onset[order(PD_onset$Chemcode),]

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

## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
r_HS2  <- 0.0052
r_HD   <- 0.0418
r_S1S2 <- r_HS2*(PD_onset$Residential)

r_S1D  <- 0.0418
r_S2D  <- 0.0418

## State rewards
# Costs
new<-matrix(seq(0, 100, by = 10), nrow = 11, ncol = 1)

c_H_nopest       <- new*1         # annual cost of being Healthy unexposed
c_S1_nopest      <- new*1         # annual cost of being Healthy exposed
c_S2_nopest      <- new*1+23041   # annual cost of being PD Onset

c_H       <- 0      # annual cost of being Healthy unexposed
c_S1      <- 0         # annual cost of being Healthy exposed
c_S2      <- 23041     # annual cost of being PD Onset
c_D       <- 0         # annual cost of being dead

# Utilities
u_H    <- 1        # annual utility of being Healthy unexposed
u_S1   <- 1        # annual utility of being Healthy exposed
u_S2   <- 0.74     # annual utility of being PD Onset
u_D    <- 0        # annual utility of being Dead

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
    pest_cohort[[x]][t + 1, ] <- pest_cohort[[x]][t, ] %*% pest[[x]]}
}

results <- vector("list", length(c_H_nopest))
# Loop through each value of c_H_nopest and run the simulation
for (k in 1:length(c_H_nopest)) {
  # Set the value of c_H_nopest for this iteration
  c_H_nopestnew <- c_H_nopest[k]
  c_S1_nopestnew <- c_S1_nopest[k]
  c_S2_nopestnew <- c_S2_nopest[k]
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
  v_c_nopest    <- c(H  = c_H_nopestnew, 
                     S1   = c_S1_nopestnew,
                     S2   = c_S2_nopestnew, 
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
  
  ########################### Cost-effectiveness analysis ########################
  ### Calculate incremental cost-effectiveness ratios (ICERs)
  df_cea <-as.data.frame(cbind(v_tot_cost, v_tot_qaly))
  df_cea$'Inc_Cost'<-(df_cea$v_tot_cost-df_cea$v_tot_cost[which(row.names(df_cea)=="no pest")])
  df_cea$'Inc_Effect'<-(df_cea$v_tot_qaly-df_cea$v_tot_qaly[which(row.names(df_cea)=="no pest")])
  results[[k]]<-df_cea
}



results_res<-do.call(cbind.data.frame, results)
results_res<-results_res[,which(colnames(results_res)=="Inc_Cost")]



###########
results_res_single<-results_res

results_res_double<-results_res

results_res_triple<-results_res

results_res_single[4,]<-results_res_double[4,]
results_res_single[8,]<-results_res_double[8,]
results_res_single[9,]<-results_res_double[9,]
results_res_single[10,]<-results_res_triple[10,]



write.csv(results_res_single,"Manuscript/uni_sen_cost_res.csv", row.names = TRUE)






################################# Model input ##################################
#Occupational
rm(list = ls())    # remove any variables in R's memory 

set.seed(123234)
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
# sort by Chemcode
exposure_prob <- exposure_prob[order(exposure_prob$Chemcode),]
PD_onset <- PD_onset[order(PD_onset$Chemcode),]

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

## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
r_HS2  <- 0.0052
r_HD   <- 0.0418
r_S1S2 <- r_HS2*(PD_onset$Occupational)

r_S1D  <- 0.0418
r_S2D  <- 0.0418

## State rewards
new<-matrix(seq(0, 100, by = 10), nrow = 11, ncol = 1)

c_H_nopest       <- new*3         # annual cost of being Healthy unexposed
c_S1_nopest      <- new*3         # annual cost of being Healthy exposed
c_S2_nopest      <- new*3+23041   # annual cost of being PD Onset

c_H       <- 0      # annual cost of being Healthy unexposed
c_S1      <- 0         # annual cost of being Healthy exposed
c_S2      <- 23041     # annual cost of being PD Onset
c_D       <- 0         # annual cost of being dead

# Utilities
u_H    <- 1        # annual utility of being Healthy unexposed
u_S1   <- 1        # annual utility of being Healthy exposed
u_S2   <- 0.74     # annual utility of being PD Onset
u_D    <- 0        # annual utility of being Dead

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
    pest_cohort[[x]][t + 1, ] <- pest_cohort[[x]][t, ] %*% pest[[x]]}
}

results <- vector("list", length(c_H_nopest))
# Loop through each value of c_H_nopest and run the simulation
for (k in 1:length(c_H_nopest)) {
  # Set the value of c_H_nopest for this iteration
  c_H_nopestnew <- c_H_nopest[k]
  c_S1_nopestnew <- c_S1_nopest[k]
  c_S2_nopestnew <- c_S2_nopest[k]
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
  v_c_nopest    <- c(H  = c_H_nopestnew, 
                     S1   = c_S1_nopestnew,
                     S2   = c_S2_nopestnew, 
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
  
  ########################### Cost-effectiveness analysis ########################
  ### Calculate incremental cost-effectiveness ratios (ICERs)
  df_cea <-as.data.frame(cbind(v_tot_cost, v_tot_qaly))
  df_cea$'Inc_Cost'<-(df_cea$v_tot_cost-df_cea$v_tot_cost[which(row.names(df_cea)=="no pest")])
  df_cea$'Inc_Effect'<-(df_cea$v_tot_qaly-df_cea$v_tot_qaly[which(row.names(df_cea)=="no pest")])
  results[[k]]<-df_cea
}



results_occ<-do.call(cbind.data.frame, results)
results_occ<-results_occ[,which(colnames(results_occ)=="Inc_Cost")]



###########
results_occ_single<-results_occ

results_occ_double<-results_occ

results_occ_triple<-results_occ

results_occ_single[4,]<-results_occ_double[4,]
results_occ_single[8,]<-results_occ_double[8,]
results_occ_single[9,]<-results_occ_double[9,]
results_occ_single[10,]<-results_occ_triple[10,]



write.csv(results_occ_single,"Manuscript/uni_sen_cost_occ.csv", row.names = TRUE)
