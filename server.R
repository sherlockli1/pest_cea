library(shiny)
library(reshape2)
library(data.table)
library(ggplot2)
library(dplyr)
library(DT)
library(kableExtra)

function(input, output) {
  ## General setup
  cycle_length <- 1       # cycle length equal to one year (use 1/12 for monthly)
  v_names_states <- c("H",  # the 4 health states of the model:
                      "S1", # Healthy Unexposed (H), Healthy Exposed (S1), PD Onset (S2), Dead (D)
                      "S2",
                      "D")
  n_states    <- length(v_names_states)     # number of health states 
  # Discounting factors
  d_c <- 0.03 # annual discount rate for costs 
  d_e <- 0.03 # annual discount rate for QALYs
  ## Transition probabilities (annual), and hazard ratios (HRs) with Pesticide Ban
  r_HS2  <- 0.0052
  r_HD   <- 0.0418
  r_S1D  <- 0.0418
  r_S2D  <- 0.0418
  # Utilities
  u_H    <- 1              # annual utility of being Healthy unexposed
  u_S1   <- 1              # annual utility of being Healthy exposed
  u_D    <- 0              # annual utility of being Dead
  u_S2   <-0.74
  # Replace the following placeholder code with your actual calculations
  v_m_init <- c(H = 1, S1 = 0, S2 = 0, D = 0) # initial state vector
  
  #######################Get No Pesticide Scenarios
  output_targeted<-reactive({
    # Within-cycle correction (WCC) using Simpson's 1/3 rule
    v_cycles <- seq(1, input$n_cycles + 1)
    v_wcc <- ((v_cycles%%2) == 0) * (2/3) + ((v_cycles%%2) != 
                                               0) * (4/3)
    v_wcc[1] <- v_wcc[input$n_cycles + 1] <- 1/3
    # Discount weight for costs and effects
    v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:input$n_cycles))
    v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:input$n_cycles))
    n_cycles <- input$n_cycles
    r_S1S2 <- r_HS2
    ## Initialize cohort trace
    pest_cohort <- matrix(NA,
                          nrow = (n_cycles + 1), ncol = n_states, 
                          dimnames = list(0:n_cycles, v_names_states))
    # Store the initial state vector in the first row of the cohort trace
    pest_cohort[1, ] <- v_m_init
    ## Initialize transition probability matrix for Pesticide Ban
    # all transitions to a non-death state are assumed to be conditional on survival 
    pest <- matrix(0, 
                   nrow = n_states, ncol = n_states, 
                   dimnames = list(v_names_states, 
                                   v_names_states)) # define row and column names
    ## Initialize transition probability for Pesticide
    # Update transition probabilities
    ## Fill in matrix
    pest["H", "H"]   <- 1 -r_HS2-r_HD
    pest["H", "S2"]  <- r_HS2
    pest["H", "D"]   <- r_HD
    # From S2
    pest["S2", "S2"] <- 1 - r_S2D
    pest["S2", "D"]  <- r_S2D
    # From D
    pest["D", "D"]   <- 1
    ####################### Run Markov model #######################
    # Iterative solution of time-independent cSTM
    for(t in 1:n_cycles){
      pest_cohort[t + 1, ] <- pest_cohort[t, ] %*% pest}  
    return(pest_cohort)
  })
  #######################Get Pesticide Scenarios
  output_case<-reactive({
    # Within-cycle correction (WCC) using Simpson's 1/3 rule
    v_cycles <- seq(1, input$n_cycles + 1)
    v_wcc <- ((v_cycles%%2) == 0) * (2/3) + ((v_cycles%%2) != 
                                               0) * (4/3)
    v_wcc[1] <- v_wcc[input$n_cycles + 1] <- 1/3
    # Discount weight for costs and effects
    v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:input$n_cycles))
    v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:input$n_cycles))
    n_cycles <- input$n_cycles
    r_S1S2 <- r_HS2*input$effect_pest
    ## Initialize cohort trace
    pest_cohort <- matrix(NA,
                          nrow = (n_cycles + 1), ncol = n_states, 
                          dimnames = list(0:n_cycles, v_names_states))
    # Store the initial state vector in the first row of the cohort trace
    pest_cohort[1, ] <- v_m_init
    pest_cohort[1,"S1"]<- input$prev_pest
    pest_cohort[1,"H"]<- 1-input$prev_pest
    ## Initialize transition probability matrix for Pesticide Ban
    # all transitions to a non-death state are assumed to be conditional on survival 
    pest <- matrix(0, 
                   nrow = n_states, ncol = n_states, 
                   dimnames = list(v_names_states, 
                                   v_names_states)) # define row and column names
    ## Initialize transition probability for Pesticide
    # Update transition probabilities
    ## Fill in matrix
    pest["H", "H"]   <- 1 -r_HS2-r_HD
    pest["H", "S2"]  <- r_HS2
    pest["H", "D"]   <- r_HD
    # From S1
    pest["S1", "S1"] <- 1 -r_S1S2- r_S1D
    pest["S1", "S2"] <- r_S1S2
    pest["S1", "D"]  <- r_S1D
    # From S2
    pest["S2", "S2"] <- 1 - r_S2D
    pest["S2", "D"]  <- r_S2D
    # From D
    pest["D", "D"]   <- 1
    ####################### Run Markov model #######################
    # Iterative solution of time-independent cSTM
    for(t in 1:n_cycles){
      pest_cohort[t + 1, ] <- pest_cohort[t, ] %*% pest}
    return(pest_cohort)
  })
  #######################Get CEA results
  output_cea<-reactive({
    ## Store the cohort
    cohort <- output_case()
    cohort_nopest <- output_targeted()
    v_dwc  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:input$n_cycles))
    v_dwe  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:input$n_cycles))
    v_cycles <- seq(1, input$n_cycles + 1)
    v_wcc <- ((v_cycles%%2) == 0) * (2/3) + ((v_cycles%%2) != 
                                               0) * (4/3)
    v_wcc[1] <- v_wcc[input$n_cycles + 1] <- 1/3    #### State Rewards scaled by the cycle length ####
    ## Vector of state utilities
    utility    <- c(H  = u_H, 
                    S1   = u_S1, 
                    S2   = u_S2,
                    D    = u_D) * cycle_length
    ## Vector of state costs
    cost    <- c(H  = 0, 
                 S1   = 0,
                 S2   = input$cost_pd, 
                 D    = 0) * cycle_length
    
    cost_nopest    <- c(H  = 0+input$cost_pest, 
                        S1   = 0+input$cost_pest,
                        S2   = input$cost_pd+input$cost_pest, 
                        D    = 0+input$cost_pest) * cycle_length
    ## create empty vectors to store total utilities and costs 
    v_tot_qaly <- v_tot_cost <- vector(mode = "numeric", length = 1)
    v_tot_qaly_nopest <- v_tot_cost_nopest <- vector(mode = "numeric", length = 1)
    
    #No pesticde    
    #* Apply state rewards 
    v_qaly_str_nopest <- cohort_nopest %*% utility # sum the utilities of all states for each cycle
    v_cost_str_nopest <- cohort_nopest %*% cost_nopest # sum the costs of all states for each cycle
    ####* Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable
    #* QALYs
    v_tot_qaly_nopest <- t(v_qaly_str_nopest) %*% (v_dwe * v_wcc)
    #* Costs
    v_tot_cost_nopest <- t(v_cost_str_nopest) %*% (v_dwc * v_wcc)
    
    #Pesticide
    #* Apply state rewards 
    v_qaly_str <- cohort %*% utility # sum the utilities of all states for each cycle
    v_cost_str <- cohort %*% cost # sum the costs of all states for each cycle
    ####* Discounted total expected QALYs and Costs per strategy and apply half-cycle correction if applicable
    #* QALYs
    v_tot_qaly <- t(v_qaly_str) %*% (v_dwe * v_wcc)
    #* Costs
    v_tot_cost <- t(v_cost_str) %*% (v_dwc * v_wcc)
    
    ICER<-(v_tot_cost-v_tot_cost_nopest)/(v_tot_qaly-v_tot_qaly_nopest)
    
    final_proportion=cohort[21,3]*input$pop
    final_proportion_nopest=cohort_nopest[21,3]*input$pop
    df_cea_1 <-as.data.frame(cbind(final_proportion, v_tot_cost,ICER))
    colnames(df_cea_1)<-c("Number of PD cases","Costs per Person ($)","Cost Reduction per Quality Adjusted Life Year")
    df_cea_2 <-as.data.frame(cbind(final_proportion_nopest, v_tot_cost_nopest,NA))
    colnames(df_cea_2)<-c("Number of PD cases","Costs per Person ($)","Cost Reduction per Quality Adjusted Life Year")
    df_cea<-rbind(df_cea_1,df_cea_2)
    df_cea <- round(df_cea, digits = 0)
    rownames(df_cea)<-c("Exposure to Pesticide","Pesticide Ban")
    return(df_cea)
  })  
  
  ########################## Generate the plot based on the calculated results
  output$output_cea <- renderTable({
    table_cea <- output_cea()  # Replace 'output_cea()' with your actual dataframe
    for (col in names(table_cea)) {
      if (is.numeric(table_cea[[col]])) {
        table_cea[[col]] <- format(table_cea[[col]], big.mark = ",", nsmall = 0)
      }
    }
    table_cea
  }, rownames = TRUE)
  output$cohort_trace <- renderPlot({
    pest_1<-as.data.frame(output_case())
    pest_1 <- tibble::rownames_to_column(pest_1, "year")
    pest_1<-melt(setDT(pest_1), id.vars = 1, variable.name = "health_state")
    pest_1$year<-as.numeric(pest_1$year)
    pest_1$health_state<-ifelse(pest_1$health_state=="H",
                                "Healthy Unexposed",
                                ifelse(pest_1$health_state=="S1",
                                       "Healthy Exposed to Pesticide",
                                       ifelse(pest_1$health_state=="S2",
                                              "PD Onset",
                                              ifelse(pest_1$health_state=="D",
                                                     "Dead",pest_1$health_state))))
    p1<-ggplot(pest_1, aes(x=year, y=value, group=health_state)) +
      geom_line(aes(color=health_state),linewidth=1)+
      scale_color_manual(values=c("#D81B60", "#0072B2", "#FFC107", "#004D40","#D55E00"))+
      theme_minimal()+ 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      labs(x ="Year", y = "Proportion of Population",color='Health Status',
           title="Change in Proportion of the Population in Each Health Status")
    p1
  })
}
