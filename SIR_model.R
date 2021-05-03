

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


## Writing the function of the SIR model
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
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

# define a function to calculate the residual sum of squares
# (RSS), passing in parameters beta and gamma that are to be
# optimised for the best fit to the incidence data

RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = initial, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((Infected - fit)^2)
}


# now find the values of beta and gamma that give the
# smallest RSS, which represents the best fit to the data.
# Start with values of 0.5 for each, and constrain them to
# the interval 0 to 1.0

Opt <- optim(c(0.5, 0.5),
             RSS,
             method = "L-BFGS-B",
             lower = c(0, 0),
             upper = c(1, 1)
)

## Now we can examine the fitted values for β and γ:
Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par

# initial values for phase 1 (1 sep to 30 oct)
initial <- c(S = N-Infected[1], I =Infected[1] , R =0)

# time in days for predictions
t <- 1:as.integer(ymd(sir_end_date) + 1 - ymd(sir_start_date))

## get the fitted values from our SIR model
predictions <- data.frame(ode(
  y = initial, times = t,
  func = SIR, parms = Opt_par
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
    subtitle = "(Red = fitted from SIR model, blue = observed)"
  ) +
  theme_minimal()


# phase 2: oct 31 to jan 29 (gives convergence of RSS)
sir_start_date <- "2020-10-31"
sir_end_date <- "2021-01-29"

## putting the daily prevalence in a vector

Infected <- subset(df1, DATE >= ymd(sir_start_date) & DATE <= ymd(sir_end_date))$cases
recovered <- subset(df1, DATE >= ymd(sir_start_date) & DATE <= ymd(sir_end_date))$recovered

# Create an incrementing Day vector the same length as our
# cases vector

Day <- 1:(length(Infected))
N <- 11589623

# initial values for phase 2 (31/10/2020 to 29/01/2021)
initial <- c(S = N-Infected[1]-recovered[1], I =Infected[1] , R =recovered[1])

# define a function to calculate the residual sum of squares
# (RSS), passing in parameters beta and gamma that are to be
# optimised for the best fit to the incidence data

RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = initial, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((Infected - fit)^2)
}


# now find the values of beta and gamma that give the
# smallest RSS, which represents the best fit to the data.
# Start with values of 0.5 for each, and constrain them to
# the interval 0 to 1.0

Opt <- optim(c(0.5, 0.5),
             RSS,
             method = "L-BFGS-B",
             lower = c(0, 0),
             upper = c(1, 1)
)

## Now we can compute the fitted values for β and γ:
Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par

# time in days for predictions
t <- 1:as.integer(ymd(sir_end_date) + 1 - ymd(sir_start_date))

## get the fitted values from our SIR model
predictions <- data.frame(ode(
  y = initial, times = t,
  func = SIR, parms = Opt_par
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
    subtitle = "(Red = fitted from SIR model, blue = observed)"
  ) +
  theme_minimal()