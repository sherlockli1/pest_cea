library(shiny)
library(reshape2)
library(data.table)
library(ggplot2)
library(dplyr)
library(DT)
library(kableExtra)

fluidPage(
  titlePanel("Cost-Effectiveness Analysis of Regulating Pesticide in California"),
  sidebarLayout(
    sidebarPanel(
      # Input controls for model parameters
      numericInput("prev_pest", "Prevalence of Pesticide:", value = 0.2, min=0,max=1,step=0.1),
      numericInput("n_cycles", "Years:", value = 20, min=0,step=1),
      numericInput("effect_pest", "Effect of Pesticide:", value = 1.1, min=1, max=3,step=0.1),
      numericInput("cost_pest", "Cost of Pesticide Ban ($):", value = 10, min=0,max=200,step=10),
      numericInput("pop", "Population Size:", value = 39440, min=0,max=100000,step=1000),
      numericInput("cost_pd", "Cost of Parkinson's disease ($):", value = 22553,min=0,max=100000,step=1000)
    ),
    mainPanel(
      # Output displays for results
      plotOutput("cohort_trace"),
      tableOutput("output_cea")
    )
  )
)
