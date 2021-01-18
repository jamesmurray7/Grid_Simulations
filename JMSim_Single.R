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
gamma <- c(1, 0.1, 0.25) # Latent association

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

# Find nearest positive-definite matrix to the one created 
sigma.corrected <- as.matrix(Matrix::nearPD(sigma)$mat)

RE <- MASS::mvrnorm(m, c(0, 0, 0), Sigma = sigma.corrected)
U.0 <- RE[, 1]; U.1 <- RE[, 2]; U.2 <- RE[, 3] # Intercept; Slope; Quadratic

# Covariates - baseline //
id <- 1:m
x1 <- rbinom(m, 1, 0.5) # Treatment received
x2 <- gl(3, 1, m) # Factor
x3 <- floor(rnorm(m, 65, 10)) # Age
x3s <- scale(x3, scale = F)

# Longitudinal part -----
# Coefficients
b0 <- 40; b1 <- -10; b22 <- 5; b23 <- 15; b3 <- 0.1
Bl <- matrix(c(b0, b1, b22, b23, b3), nrow = 1) # Cast to matrix
# Baseline covariates
x1_l <- rep(x1, each = n_i)
x2_l <- rep(x2, each = n_i)
x3_l <- rep(x3s, each = n_i)
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
XsBs <- cbind(x1,x3s) %*% c(-0.1, 0.05)

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
  
  candidate.times <- temp[lambda.tdt > temp[,7], 3]
  candidate.times <- c(candidate.times, (n_i-1)) # Adding truncation time in for safety.
  surv.times[i] <- min(candidate.times)
  message(i, " done----\n")
}
hist(surv.times)

# Censoring
cens.times <- rexp(m, 0.1) # Generate censoring times
f <- cens.times < surv.times # Flagging variable for when X_i < C_i
status <- rep(1,m) # Status vector (1 = dead, 0 = censored)
status[f] <- 0 # Replace those censored


surv.times[f] <- cens.times[f] # Change times
status[surv.times == (n_i-1)] <- 0 # If time = truncation time then censored

cbind(surv.times, cens.times, status, f) # Make sure this looks right


surv.data <- data.frame(surv.times, x1, x3, x3s, id, status)
phs <- coxph(Surv(surv.times, status) ~ x1 + x3s, data = surv.data)
summary(phs)


# Single-run of jointdata() and joint()
library(joineR)
jd <- jointdata(
  longitudinal = long_data,
  survival = surv.data,
  baseline = surv.data[,c("id", "x1", "x3s")],
  id.col = "id", time.col = "time"
)

jfit <- joint(jd, model = "quad", max.it = 500,
              long.formula = Y ~ x1_l + x2_l + x3_l + time,
              surv.formula = Surv(surv.times, status) ~ x1 + x3s,
              sepassoc = T)

summary(jfit)



