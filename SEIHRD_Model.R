
## load libraries
library(readxl)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(deSolve)
library(bbmle)
library(tidyverse)
library(ggsci)

## load data
hosp <- read_csv("https://epistat.sciensano.be/Data/COVID19BE_HOSP.csv")
cases <- read_csv("https://epistat.sciensano.be/Data/COVID19BE_CASES_AGESEX.csv")
death <- read_csv("https://epistat.sciensano.be/Data/COVID19BE_MORT.csv")

## check data and structure

sapply(hosp, class)
sapply(cases, class)

## keep only data for BELGIUM

hosp <- hosp %>%
  group_by(DATE) %>%
  summarise(N_hosp = sum(TOTAL_IN),
            N_ICU = sum(TOTAL_IN_ICU)) %>%
  mutate(DATE_INT = as.integer(as.factor(DATE)))

cases <- cases %>%
  group_by(DATE) %>%
  summarise(N_cases = sum(CASES)) %>%
  mutate(DATE_INT = as.integer(as.factor(DATE)))

death <- death %>%
  group_by(DATE) %>%
  summarise(N_death = sum(DEATHS)) %>%
  mutate(DATE_INT = as.integer(as.factor(DATE)))

# First analysis dataset

cases_wave1_prev <- cases %>%
  filter(DATE < ymd("2020-09-01")) %>%
  mutate(prevalence = N_cases * 4)

cases_wave1_prev$smooth <- rollmean(cases_wave1_prev$prevalence, k = 7, fill = NA)

flu_dat1 <- cases_wave1_prev %>%
  filter(complete.cases(.)) %>%
  mutate(DATE_INT = as.integer(as.factor(DATE))-1) %>%
  select(DATE, time = DATE_INT, I = smooth, prevalence)

# Second analysis dataset

sir_start_date <- "2020-09-01"
sir_end_date <- "2021-01-31"

cases_wave2 <- cases %>%
  filter(DATE < ymd(sir_end_date),
         DATE > ymd(sir_start_date))

hosp_wave2 <- hosp %>%
  filter(DATE < ymd(sir_end_date),
         DATE > ymd(sir_start_date))

death_wave2 <- death %>%
  filter(DATE < ymd(sir_end_date),
         DATE > ymd(sir_start_date)) %>%
  mutate(cumm_death = cumsum(N_death))

cases_wave2 <- cases_wave2 %>%
  mutate(prevalence = 4 * N_cases)

wave2 <- inner_join(x = cases_wave2, y = hosp_wave2[c("DATE", "N_hosp")], by = "DATE")  

wave2 <- inner_join(x = wave2, y = death_wave2[c("DATE", "N_death", "cumm_death")], by = "DATE")

wave2$smooth <- rollmean(wave2$prevalence, k = 7, fill = NA)

flu_dat2 <- wave2 %>%
  filter(complete.cases(.)) %>%
  mutate(DATE_INT = as.integer(as.factor(DATE))-1) %>%
  select(DATE, time = DATE_INT, I = smooth, H = N_hosp,D = cumm_death, prevalence, N_death)

##  DESCRIPTIVE PLOTS
# Prevalence graph during the first wave
flu_dat1 %>%
  ggplot(aes(x=DATE, y=prevalence)) +
  geom_point(colour="red") +
  scale_x_date(date_breaks="1 month", date_labels = "%d %b") +
  scale_color_aaas() +
  labs(y="Prevalence") +
  theme_minimal()

# Prevalence graph during the second wave
flu_dat2 %>%
  select(-time, -H, -D) %>%
  ggplot(aes(x=DATE, y=I)) +
  geom_point(color = "red") +
  scale_color_aaas()+
  labs(y="Prevalence") +
  theme_minimal()

# Hospitalisation graph
hosp_wave2 %>%
  select(-time)%>%
  pivot_longer(-1) %>%
  ggplot() +
  geom_line(aes(x = DATE, y = value, color =name)) +
  scale_color_manual(name="Hospitalisations", values = c("red", "dodgerblue")) +
  labs(y="Number of Hospitalisations") +
  theme_minimal()

# Cumulative deaths graph
flu_dat2 %>%
  select(-time, -I, -H) %>%
  ggplot(aes(x=DATE, y=D)) +
  geom_point(color = "dark green") +
  scale_color_aaas()+
  labs(y="Cumulative deaths") +
  theme_minimal() 

## SIR MODEL
# SIR algorithm
SIR_fn <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    N  <- S+I+R
    
    dS <- -beta*S*I/N
    dI <- beta*S*I/N-gamma*I
    dR <- gamma*I
    
    return(list(c(dS, dI, dR)))
    
  })
  
}

# Initial values
initial_state_values <- c(S = 10000, I = 1, R = 0)

beta_start  <- 1
gamma_start <- 0.5

times <- seq(from = 0, to = 178, by = 0.01)

# Function for the SSE
SIR_SSQ <- function(parameters, dat) {  # parameters must contain beta & gamma
  
  # calculate model output using your SIR function with ode()
  
  result <- as.data.frame(ode(y = initial_state_values  # vector of initial state 
                              # values, with named elements
                              , times = times             # vector of times
                              , func = SIR_fn             # your predefined SIR function
                              , parms = parameters)       # the parameters argument
                          # entered with SIR_SSQ()
  )
  
  # SSQ calculation: needs the dat argument (the observed data you are fitting to)
  # assumes the data you are fitting to has a column "I"
  
  # select only complete cases, i.e. rows with no NAs, from the dataframe
  dat <- na.omit(dat)  
  
  # select elements where results$time is in dat$time
  deltas2 <- (result$I[result$time %in% dat$time]  
              - dat$I)^2                             
  SSQ   <- sum(deltas2)
  
  return(SSQ)
  
}

# Optimization 

optimised <- optim(par = c(beta = beta_start
                           , gamma = gamma_start)     
                   
                   , fn = SIR_SSQ
                   , dat = flu_dat1)

# Computing fitted values

opt_mod1 <- as.data.frame(ode(y = initial_state_values  
                              , times = times            
                              ,  func = SIR_fn           
                              , parms = optimised$par))

# Goodness-of-fit Assessment
opt_plot <- ggplot()
opt_plot <- opt_plot + geom_point(aes(x = time, y = prevalence)
                                  , colour = "red"
                                  , data = flu_dat1)
opt_plot <- opt_plot + geom_line(aes(x = time, y = I)
                                 , colour = "blue"
                                 , data   = opt_mod1)+
  labs(y="Prevalence",
       x= "Time (Days)")+
  theme_minimal() 
opt_plot

## SEIHRD MODEL
# SEIHRD algorithm
SEIR_fn <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    N  <- S+E+I+H+R+M
    
    dS <- -beta*S*I/N
    dE <- beta*S*I/N-delta*E
    dI <- (delta * E) - (gamma * I) 
    dH <- rho*(gamma * I) - (mu*H) - (lambda*H) 
    dM <- (mu*H)
    dR <- (1-rho)*(gamma*I) + (lambda*H) 
    
    return(list(c(dS, dE, dI, dH, dM, dR)))
    
  })
  
}

# Model initiation
initial_state_values <- c(S = 500000, E = 1, I = 0, H = 0, M = 0, R = 0)

beta_start  <- 1
gamma_start <- 0.5
delta  <- 1/5.2
mu_start <- 1/50
lambda <- 1/12.4
rho <- 0.05

times <- seq(from = 0, to = 145, by = 0.01)

# SSE function
SEIR_SSQ <- function(parameters, dat) { 
  
  result <- as.data.frame(ode(y = initial_state_values 
                              , times = times             
                              , func = SEIR_fn             
                              , parms = parameters)       
  )
  

  dat <- na.omit(dat)  
  
 
  deltas2 <- (result$I[result$time %in% dat$time]  
              - dat$I)^2                             
  SSQ1   <- sum(deltas2)
  
  deltas3 <- (result$H[result$time %in% dat$time]  
              - dat$H)^2                             
  SSQ2   <- sum(deltas3)
  
  deltas4 <- (result$M[result$time %in% dat$time]  
              - dat$D)^2                             
  SSQ3   <- sum(deltas4)
  
  SSQ <- SSQ1 + SSQ2 + SSQ3  
  return(SSQ)
  
}

# Model optimization

optimised <- optim(par = c(beta = beta_start, 
                           gamma = gamma_start,
                           mu = mu_start)      
                   
                   , fn = SEIR_SSQ
                   , dat = flu_dat2  
)

# Computing fitted values
opt_mod <- as.data.frame(ode(y = initial_state_values  
                             , times = times            
                             ,  func = SEIR_fn            
                             , parms = optimised$par))

## Goodness-of-fit Assessment
# Infectious state
opt_plot <- ggplot()
opt_plot <- opt_plot + geom_point(aes(x = time, y = prevalence)
                                  , colour = "red"
                                  , data = flu_dat2)
opt_plot <- opt_plot + geom_line(aes(x = time, y = I)
                                 , colour = "blue"
                                 , data   = opt_mod) +
  labs(y="Prevalence",
       x= "Time (Days)") +
  theme_minimal() 
opt_plot

# Hospitalisation state
opt_plot <- ggplot()
opt_plot <- opt_plot + geom_point(aes(x = time, y = H)
                                  , colour = "red"
                                  , data = flu_dat2)
opt_plot <- opt_plot + geom_line(aes(x = time, y = H)
                                 , colour = "blue"
                                 , data   = opt_mod) +
  labs(y="Number of Hospitalisations (Prevalence)",
       x= "Time (Days)") +
  theme_minimal() 
opt_plot

# Death state

opt_plot <- ggplot()
opt_plot <- opt_plot + geom_point(aes(x = time, y = D)
                                  , colour = "red"
                                  , data = flu_dat2)
opt_plot <- opt_plot + geom_line(aes(x = time, y = M)
                                 , colour = "blue"
                                 , data   = opt_mod) +
  labs(y="Cumulative deaths",
       x= "Time (Days)") +
  theme_minimal() 
opt_plot

## Modelling Intervention measures during first wave
# Model functions for different intervention measures
SIR_lokdown1 <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    beta = ifelse(
      (time <= start_lockdown || time >= end_lockdown),
      0.484, 0.05
    )
    N  <- S+I+R
    
    dS <- -beta*S*I/N
    dI <- beta*S*I/N-gamma*I
    dR <- gamma*I
    
    return(list(c(dS, dI, dR)))
    
  })
  
}


SIR_lokdown2 <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    beta = ifelse(
      (time <= start_lockdown || time >= end_lockdown),
      0.484, 0.1
    )
    N  <- S+I+R
    
    dS <- -beta*S*I/N
    dI <- beta*S*I/N-gamma*I
    dR <- gamma*I
    
    return(list(c(dS, dI, dR)))
    
  })
  
}

SIR_lokdown3 <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    beta = ifelse(
      (time <= start_lockdown || time >= end_lockdown),
      0.484, 0.2
    )
    N  <- S+I+R
    
    dS <- -beta*S*I/N
    dI <- beta*S*I/N-gamma*I
    dR <- gamma*I
    
    return(list(c(dS, dI, dR)))
    
  })
  
}

# Initial values
initial_state_values <- c(S = 10000, I = 1, R = 0)

beta_start  <- 1
gamma_start <- 0.5

times <- seq(from = 0, to = 178, by = 0.01)

# Computing fitted values for different interventionmeasures
params <- c(
  gamma=0.04,
  start_lockdown=18,
  end_lockdown=72
)
predictions1 <- as.data.frame(ode(y = initial_state_values  
                                  
                                  , times = times           
                                  ,  func = SIR_lokdown1          
                                  , parms = params))

predictions1 <- predictions1 %>%
  select(time, "β = 0.05" = I)


predictions2 <- as.data.frame(ode(y = initial_state_values  
                                  
                                  , times = times           
                                  ,  func = SIR_lokdown2         
                                  , parms = params))

predictions2 <- predictions2 %>%
  select(time, "β = 0.1" = I)

predictions3 <- as.data.frame(ode(y = initial_state_values  
                                  
                                  , times = times           
                                  ,  func = SIR_lokdown3          
                                  , parms = params))

predictions3 <- predictions3 %>%
  select(time, "β = 0.2" = I)

# Keeping only required variables
opt_mod11 <- opt_mod1 %>%
  select(time, Baseline  = I)

wave1_lockdown <- inner_join(x = opt_mod11, y = predictions1[c("time", "β = 0.05")], by = "time")  

wave1_lockdown <- inner_join(x = wave1_lockdown, y = predictions2[c("time", "β = 0.1")], by = "time")

wave1_lockdown <- inner_join(x = wave1_lockdown, y = predictions3[c("time", "β = 0.2")], by = "time")

# Graph with intervention measures
wave1_lockdown %>%
  pivot_longer(-1) %>%
  ggplot() +
  geom_line(aes(x = time, y = value, color =name)) +
  scale_color_manual(name="Degree of Intervention measures",
                     values = c("red","green", "dodgerblue","purple")) +
  labs(y="Prevalence",
       x="time (Days)") +
  theme_minimal() 

## Forecasting of Hospitalisations and deaths
# Processing the data
initial_state_values <- c(S = 500000, E = 1, I = 0, H = 0, M = 0, R = 0)
times <- seq(from = 0, to = 200, by = 0.01)

opt_mod <- as.data.frame(ode(y = initial_state_values  
                             , times = times            
                             ,  func = SEIR_fn            
                             , parms = optimised$par))

opt_mod_1 <- as.data.frame(ode(y = initial_state_values  
                               , times = times            
                               ,  func = SEIR_fn           
                               , parms = c(beta = 0.704, gamma = 0.37, delta = 1/5.2, mu = 0.04, lambda = 0.089)))

opt_mod_2 <- as.data.frame(ode(y = initial_state_values  
                               , times = times            
                               ,  func = SEIR_fn           
                               , parms = c(beta = 1.76, gamma = 0.518, delta = 1/5.2, mu = 0.04, lambda = 0.097)))

opt_mod_3 <- as.data.frame(ode(y = initial_state_values  
                               , times = times            
                               ,  func = SEIR_fn           
                               , parms = c(beta = 0.53, gamma = 0.481, delta = 1/5.2, mu = 0.04, lambda = 1/12.4)))

opt_mod_4 <- as.data.frame(ode(y = initial_state_values  
                               , times = times            
                               ,  func = SEIR_fn           
                               , parms = c(beta = 1.41, gamma = 0.592, delta = 1/5.2, mu = 0.04, lambda = 0.113)))

# Keeping only the required variables for hospitalisations
opt_mod_1H <- opt_mod_1 %>%
  select(time, "Scenario 1: 60% reduction in β; 10% increase in ʎ" = H)

opt_mod_2H <- opt_mod_2 %>%
  select(time, "Scenario 2: 20% increase in ʎ;40% increase in ϒ" = H)

opt_mod_3H <- opt_mod_3 %>%
  select(time,"Scenario 3: 70% reduction in β;30% increase in ϒ"  = H)

opt_mod_4H <- opt_mod_4 %>%
  select(time, "Scenario 4: 20% reduction in β;60% increase in ϒ; 40% increase in ʎ" = H)

# Keeping only the required variables for deaths
opt_mod_1 <- opt_mod_1 %>%
  select(time, "Scenario 1: 60% reduction in β; 10% increase in ʎ" = M)

opt_mod_2 <- opt_mod_2 %>%
  select(time, "Scenario 2: 20% increase in ʎ;40% increase in ϒ" = M)

opt_mod_3 <- opt_mod_3 %>%
  select(time,"Scenario 3: 70% reduction in β;30% increase in ϒ"  = M)

opt_mod_4 <- opt_mod_4 %>%
  select(time, "Scenario 4: 20% reduction in β;60% increase in ϒ; 40% increase in ʎ" = M)

# Creating the analysis dataset for hospitalisations
opt_mod1 <- opt_mod %>%
  select(time, baseline  = H)

wave2_Hforecast <- inner_join(x = opt_mod1, y = opt_mod_1H[c("time","Scenario 1: 60% reduction in β; 10% increase in ʎ" )], by = "time")  

wave2_Hforecast <- inner_join(x = wave2_Hforecast, y = opt_mod_2H[c("time","Scenario 2: 20% increase in ʎ;40% increase in ϒ" )], by = "time")

wave2_Hforecast <- inner_join(x = wave2_Hforecast, y = opt_mod_3H[c("time","Scenario 3: 70% reduction in β;30% increase in ϒ" )], by = "time")

wave2_Hforecast <- inner_join(x = wave2_Hforecast, y = opt_mod_4H[c("time","Scenario 4: 20% reduction in β;60% increase in ϒ; 40% increase in ʎ")], by = "time")

# Creating the analysis dataset for deaths
opt_mod11 <- opt_mod %>%
  select(time, baseline  = M)

wave2_Mforecast <- inner_join(x = opt_mod11, y = opt_mod_1[c("time","Scenario 1: 60% reduction in β; 10% increase in ʎ" )], by = "time")  

wave2_Mforecast <- inner_join(x = wave2_Mforecast, y = opt_mod_2[c("time","Scenario 2: 20% increase in ʎ;40% increase in ϒ" )], by = "time")

wave2_Mforecast <- inner_join(x = wave2_Mforecast, y = opt_mod_3[c("time","Scenario 3: 70% reduction in β;30% increase in ϒ" )], by = "time")

wave2_Mforecast <- inner_join(x = wave2_Mforecast, y = opt_mod_4[c("time","Scenario 4: 20% reduction in β;60% increase in ϒ; 40% increase in ʎ")], by = "time")

# Graph forecasting hospitalisations in the short-term
wave2_Hforecast %>%
  pivot_longer(-1) %>%
  ggplot() +
  geom_line(aes(x = time, y = value, color =name)) +
  scale_color_manual(name="Scenario", values = c("red","brown", "dodgerblue","green","purple")) +
  labs(y="Number of Hospitalisations",
       x="Time (Days)") +
  theme(legend.position = c(0.7,0.8)) +
  theme_minimal()

# Graph forecasting deaths in the short-term
wave2_Mforecast %>%
  pivot_longer(-1) %>%
  ggplot() +
  geom_line(aes(x = time, y = value, color =name)) +
  scale_color_manual(name="Scenario", values = c("red","brown", "dodgerblue","green","purple")) +
  labs(y="Number of Hospitalisations",
       x="Time (Days)") +
  theme(legend.position = c(0.7,0.8)) +
  theme_minimal()
