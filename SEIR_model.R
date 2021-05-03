
## Importing the dataset

library(readxl)
library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(lubridate)
library(deSolve)

hosp <- read_excel("C:/Users/JEAN-CALVIN/Downloads/thesis_data/COVID19BE.xlsx",
                   sheet=4)

## Processing the data

attach(hosp)
hosp$DATE<-as.Date(hosp$DATE)

# create a data frame for prevalence 

df1 <- hosp %>%
  filter(DATE >= "2020-09-01" & DATE <= "2021-01-29") %>%
  group_by(DATE) %>%
  summarise(cases = sum(TOTAL_IN, na.rm = TRUE),
            recovered = sum(NEW_OUT)) %>%
  arrange(DATE) %>%
  ungroup()%>%
  mutate(cases_cum = cumsum(cases))

## daily prevalence graph
df1 %>%
  ggplot(aes(x = DATE)) +
  geom_line(aes(y = cases, colour = "red")) +
  labs(
    y = "Prevalence",
    title = "Daily Prevalence, Belgium"
  ) +
  theme_minimal()

## Writing the function of the SEIR model
SEIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dE <- (beta * I * S / N) - (alpha * E)
    dI <- (alpha * E) - (gamma * I)
    dR <- gamma * I
    list(c(dS, dE, dI, dR))
  })
}


## defining the period for model prediction  

## phase 1: sept to oct 30  

sir_start_date <- "2020-09-01"
sir_end_date <- "2020-10-30"

## putting the daily prevalence in a vector

Infected <- subset(df1, DATE >= ymd(sir_start_date) & DATE <= ymd(sir_end_date))$cases

# Create an incrementing Day vector the same length as our
# cases vector

Day <- 1:(length(Infected))
N <- 11589623


## Parameters values
# phase 1
beta_value = 0.5277758
gamma_value = 0.4722227
alpha_value = 1/5.0

parameter_list <- c (beta = beta_value, gamma = gamma_value, alpha = alpha_value)

## Initial values for the differential equation
I_0 <- Infected[1]
E_0 <- I_0*37
R_0 <- 0

# initial values for phase 1 (1 sep to 30 oct)
initial_values <- c(S = N-E_0-I_0-R_0, E=E_0, I=I_0,  R=R_0 )

# time in days for predictions
t <- 1:as.integer(ymd(sir_end_date) + 1 - ymd(sir_start_date))

## get the fitted values from our SEIR model
predictions <- data.frame(ode(
  y = initial_values, times = t,
  func = SEIR, parms = parameter_list
))


# add a Date column and the observed prevalence data for phase 1
fitted_prevalence1 <- predictions %>%
  mutate(
    Date = ymd(sir_start_date) + days(t- 1),
    Country = "Belgium",
    prevalence = Infected
  )

## plot the data to check for fit

fitted_prevalence1 %>%
  filter(Date <= ymd("2021-01-29")) %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_point(aes(y = prevalence), colour = "blue") +
  labs(
    y = "Prevalence",
    title = "COVID-19 fitted vs prevalence-phase 1, Belgium",
    subtitle = "(Red = fitted from SEIR model, blue = observed)"
  ) +
  theme_minimal()


# phase 2: oct 31 to jan 29 
sir_start_date <- "2020-10-31"
sir_end_date <- "2021-01-29"

## putting the daily prevalence in a vector

Infected <- subset(df1, DATE >= ymd(sir_start_date) & DATE <= ymd(sir_end_date))$cases
recovered <- subset(df1, DATE >= ymd(sir_start_date) & DATE <= ymd(sir_end_date))$recovered

# Create an incrementing Day vector the same length as our
# cases vector

Day <- 1:(length(Infected))
N <- 11589623


## Parameters values
# phase 2
beta_value = 0.8205509
gamma_value = 0.8283531
alpha_value = 1/5.0

parameter_list <- c (beta = beta_value, gamma = gamma_value, alpha = alpha_value)

## Initial values for the differential equation
I_0 <- Infected[1]
E_0 <- I_0*37
r_0 <- recovered[1]

# initial values for phase 2 (31/10/2020 to 29/01/2021)
initial_values <- c(I=I_0, E=E_0, R=r_0, S = N-E_0-I_0-r_0)

# time in days for predictions
t <- 1:as.integer(ymd(sir_end_date) + 1 - ymd(sir_start_date))

## get the fitted values from our SEIR model
predictions <- data.frame(ode(
  y = initial_values, times = t,
  func = SEIR, parms = parameter_list
))


# add a Date column and the observed prevalence data for phase 2
fitted_prevalence2 <- predictions %>%
  mutate(
    Date = ymd(sir_start_date) + days(t- 1),
    Country = "Belgium",
    prevalence = Infected
  )

## plot the data to check for fit

fitted_prevalence2 %>%
  filter(Date <= ymd("2021-01-29")) %>%
  ggplot(aes(x = Date)) +
  geom_line(aes(y = I), colour = "red") +
  geom_point(aes(y = prevalence), colour = "blue") +
  labs(
    y = "Prevalence",
    title = "COVID-19 fitted vs prevalence-phase 2, Belgium",
    subtitle = "(Red = fitted from SEIR model, blue = observed)"
  ) +
  theme_minimal()

# merging fitted values for both time-points
fitted_values <- rbind(fitted_prevalence1,fitted_prevalence2)

# fitting graph for both phases

fitted_values %>%
  ggplot(aes(x = Date)) +
  geom_point(aes(y = I), colour = "red") +
  geom_line(aes(y = prevalence), colour = "blue") +
  labs(
    y = "Prevalence",
    title = "COVID-19 fitted vs prevalence, Belgium",
    subtitle = "(Red = fitted from SEIR model, blue = observed)"
  ) +
  theme_minimal()