#' ##########################
#' JMSim_Single.R   Author: James Murray
#' --------------------------
#' Single run of simulation for quadratic random effects.
#' Longitudinal part is exact same as done in previous work.
#' Survival times are done using grid-step method.
#' --------------------------
#' Functionisation will be done in another script.
#' ##########################

# Prerequisites ------------------------------------------------------------
dev.off()
rm(list = ls())
library(tidyverse)
theme_set(theme_light())
library(lme4)
library(survival)

# Setting-out the scenario ------------------------------------------------
# Some diabetes trial, higher value of outcome is worse.
# Binary covariate is receiving treatment (yes = good)
# Factor covariate is BMI category at baseline (fat = bad)
# Continuous covariate is age at baseline (older = bad)
# Six treatment times (t) - Could choose RV Unif[10,15] for this.

# 'Global' parameters //
n_i <- 6 # Number of measurements per participant
m <- 250 # Number of participants
N <- n_i * m
sigma.e <- 1.5 # Epsilon 
gamma <- c(1, 1, 1) # Latent association

# Generate Random effects //
sigma.0 <- 3   # Intercept RE
sigma.1 <- 2   # Slope RE
sigma.2 <- 1   # Quadratic RE 
# Covariance matrix --
sigma <- matrix(
  c(sigma.0 ^ 2,       sigma.0 * sigma.1, sigma.0 * sigma.2,
    sigma.1 * sigma.0, sigma.1 ^ 2,       sigma.1 * sigma.2,
    sigma.2 * sigma.0, sigma.2 * sigma.1, sigma.2 ^ 2), nrow = 3, byrow = T
)

RE <- MASS::mvrnorm(m, c(0, 0, 0), sigma)
U.0 <- RE[, 1]; U.1 <- RE[, 2]; U.2 <- RE[, 3] # Intercept; Slope; Quadratic

# Covariates - baseline //
id <- 1:m
x1 <- rbinom(m, 1, 0.5) # Treatment received
x2 <- gl(3, 1, m) # Factor
x3 <- floor(rnorm(m, 65, 10)) # Age

# Longitudinal part -----
# Coefficients
b0 <- 40; b1 <- -10; b22 <- 5; b23 <- 15; b3 <- 0.1
Bl <- matrix(c(b0, b1, b22, b23, b3), nrow = 1) # Cast to matrix
# Baseline covariates
x1_l <- rep(x1, each = n_i)
x2_l <- rep(x2, each = n_i)
x3_l <- rep(x3, each = n_i)
Xl <- model.matrix(~x1_l+x2_l+x3_l)
time <- rep(0:(n_i-1), m)
# REs
U0l <- rep(U.0, each = n_i)
U1l <- rep(U.1, each = n_i)
epsilon <- rnorm(N, 0, sigma.e)
# Response
Y <- Xl %*% t(Bl) + U0l + U1l * time + epsilon
# Data and quick model
long_data <- data.frame(id = rep(id, each = n_i), x1_l, x2_l, x3_l, time, Y)
summary(lmer(Y ~ x1_l + x2_l + x3_l + time + (1+time|id), data = long_data)) # Cool!

# Survival part ----
theta.0 <- exp(-3)
theta.1 <- 1
XsBs <- cbind(x1,x3) %*% c(-0.1, 0.05)


# Define grid steps
dt <- 0.001
gridt <- seq(0, (n_i-1), dt)[-1] # Removing time = 0
# Cast to data.frame then matrix (Definitely not the most efficient way!)
hazards <- data.frame(
  id = rep(id, each = length(gridt)),
  XsBs = rep(XsBs, each = length(gridt)),
  time = rep(gridt, m),
  U.0 = rep(U.0, each = length(gridt)),
  U.1 = rep(U.1, each = length(gridt)),
  U.2 = rep(U.2, each = length(gridt)),
  Unif = runif(m * length(gridt))
)
hazards <- as.matrix(hazards)

# Calculate hazards at time t + dt, scale by dt to obtain Pr() and compare against Unif,
# If this probability doesn't exceed Unif, then survival time set as Tau, else t
# ID is then allocated the minimum of its possible times, tt.

surv.times <- c()
for(i in id){
  temp <- hazards[hazards[,"id"] == i, ]
  lambda.t <- exp(theta.0 + theta.1 * temp[, 3]) * exp(temp[,2] + gamma[1] * temp[,4] + gamma[2] * temp[,5] * temp[,3] + 
                                                       gamma[3] * temp[,6] * temp[,3] ^ 2)
  lambda.tdt <- lambda.t * dt
  
  candidate.times <- hazards[lambda.tdt > hazards[,7], 3]
  candidate.times <- c(candidate.times, (n_i-1)) # Adding truncation time in for safety.
  surv.times[i] <- min(candidate.times)
  message(i, " done----\n")
}
hist(surv.times)

# NB: This is how it's done in joineR code: A lot more Efficiently(!) ----
bl.haz <- exp(theta.0 + theta.1 * gridt)
fixed.haz <- dt * exp(XsBs) %*% bl.haz
gamma.t <- gamma * joineR:::getD(3, gridt)
gamma.RE <- exp(cbind(U.0,U.1,U.2) %*% gamma.t) 
lambda.t.pete <- fixed.haz * gamma.RE
# Simulate survival times
# The baseline hazard function is specified Gompertz distribution with
# theta1 shape and scale exp(theta0)


summary(coxph(Surv(survtime, status) ~ x1 + x3, data = surv_data)) # Way further off than just R.I!

# Single-run of jointdata() and joint()
library(joineR)
jd <- jointdata(
  longitudinal = long_data,
  survival = surv_data,
  baseline = surv_data[,c("id", "x1", "x3")],
  id.col = "id", time.col = "time"
)

jfit <- joint(jd,
              long.formula = Y ~ x1_l + x2_l + x3_l + time,
              surv.formula = Surv(survtime, status) ~ x1 + x3,
              sepassoc = T)


# Functionise the data simulation -----------------------------------------

joint_sim <- function(m = 200, n_i = 5, 
                      Bl = c(40, -10, 5, 15, 0.1), # Longit: Intercept, binary, factor2-3, continuous
                      Bs = c(-0.3, 0.05), # Survival: log-odds binary and continuous,
                      sigma.i = 3, sigma.s = 2, rho = 0.3, # RE parameters
                      sigma.e = 1.5, # Error parameter
                      theta0 = log(0.1), theta1 = 0.1){ # Scale and shape of Gompertz.
  # 'Global' parameters //
  N <- n_i * m
  id <- 1:m
  tau <- n_i-1 # Upper limit of time variable.
  rho <- 0.3
  sigma.i <- 3   # Random intercept SD
  sigma.s <- 2   # Random slope SD
  sigma.e <- 1.5 # Epsilon SD
  
  # Covariance matrix //
  sigma <- matrix(
    c(sigma.i ^ 2, rho * sigma.i * sigma.s,
      rho * sigma.i * sigma.s, sigma.s^2), nrow = 2, byrow = T
  )
  
  # Generate Random effects //
  RE <- MASS::mvrnorm(m, c(0,0), sigma)
  Ui <- RE[, 1]; Us <- RE[, 2] # Intercept; Slope
  
  # Covariates - baseline //
  x1 <- rbinom(m, 1, 0.5) # Treatment received
  x2 <- gl(3, 1, m) # Factor
  x3 <- floor(rnorm(m, 65, 10)) # Age
  
  # Longitudinal part //
  # Coefficients
  Bl <- matrix(Bl, nrow = 1) # Cast to matrix
  # Baseline covariates
  x1l <- rep(x1, each = n_i)
  x2l <- rep(x2, each = n_i)
  x3l <- rep(x3, each = n_i)
  Xl <- model.matrix(~x1l+x2l+x3l)
  time <- rep(0:tau, m)
  # REs
  U1l <- rep(Ui, each = n_i)
  U2l <- rep(Us, each = n_i)
  epsilon <- rnorm(N, 0, sigma.e)
  # Response
  Y <- Xl %*% t(Bl) + U1l + U2l * time + epsilon
  # Data and quick model
  long_data <- data.frame(id = rep(id, each = n_i), x1l, x2l, x3l, time, Y)
  
  # Survival part //
  Bs <- matrix(Bs, nrow = 1)
  Xs <- model.matrix(~ x1 + x3 - 1)
  
  # Simulate survival times
  uu <- runif(m)
  uu <- log(uu)
  suppressWarnings(
    tt <- log(1-(uu * (theta1 + Us))/exp(theta0 + Xs %*% t(Bs) + Ui))/(theta1 + Us)
  )
  tt[is.nan(tt)] <- Inf
  # Censoring
  censor <- rexp(m, 0.01)
  survtime <- pmin(tt, censor, 5)
  status <- ifelse(survtime == tt, 1, 0)
  
  surv_data <- data.frame(id, x1, x3, survtime, status)
  
  # Extra output - number of events
  pc_events <- length(which(survtime < tau))/m * 100
  
  return(list(long_data, surv_data, pc_events))
  
}
joint_sim()

# Separate investigations -------------------------------------------------

separate_fits <- function(df){
  lmm_fit <- lmer(Y ~ x1l + x2l + x3l + time + (1+time|id), data = df[[1]])
  surv_fit <- coxph(Surv(survtime, status) ~ x1 + x3, data = df[[2]])
  return(
    list(lmm_fit, surv_fit)
  )
}

pb <- progress::progress_bar$new(total = 1000)
longit_beta <- data.frame(beta0 = NA, beta1 = NA, beta22 = NA, beta23 = NA, beta3 = NA, 
                          sigma.e = NA, sigma.ui = NA, sigma_us = NA)
surv_beta <- data.frame(beta1s = NA, beta3s = NA)
pc_events <- c()

for(i in 1:1000){
  dat <- joint_sim()
  pc_events[i] <- dat[[3]]
  fits <- separate_fits(dat)
  long_coefs <- fits[[1]]@beta[1:5]
  long_sigma.e <- sigma(fits[[1]])
  long_sigma.u <- as.numeric(attr(VarCorr(fits[[1]])$id, "stddev"))
  longit_beta[i,] <- c(long_coefs, long_sigma.e, long_sigma.u)
  surv_beta[i, ] <- as.numeric(fits[[2]]$coefficients)
  pb$tick()
}

ex <- expression
to_plot <- cbind(longit_beta, surv_beta, pc_events) %>% tibble %>% 
  gather("parameter", "estimate") %>% 
  mutate(param = factor(parameter, levels = c("beta0", "beta1", "beta22", "beta23", "beta3", 
                                              "sigma.e", "sigma.ui", 'sigma_us',
                                              "beta1s", "beta3s", "pc_events"),
                        labels = c(ex(beta[0]), ex(beta[1]), ex(beta[22]), ex(beta[23]), ex(beta[3]),
                                   ex(sigma[e]), ex(sigma[u*"i"]), ex(sigma[u*"s"]),
                                   ex(beta[1*"S"]), ex(beta[3*"S"]), ex("Percent"*" experiencing"*" events")))
  )

plot_lines <- to_plot %>% distinct(param)
plot_lines$xint <- c(40, -10, 5, 15, 0.1, 1.5, 3, 2, -0.3, 0.05, NA)

to_plot %>% 
  ggplot(aes(x = estimate)) + 
  geom_density(fill = "grey20", alpha = .2) + 
  geom_vline(data = plot_lines, aes(xintercept = xint), colour = "blue", alpha = .5, lty = 3) + 
  facet_wrap(~param, scales = "free", nrow = 6, ncol = 2, labeller = label_parsed) + 
  labs(title = "Separate investigation", x = "Estimate")
ggsave("./JM-sims-plots/Separate_Investigation_Redux.png")

