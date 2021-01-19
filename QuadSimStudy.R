#' ##########################
#' QuadSimStudy.R   Author: James Murray
#' --------------------------
#' Simulation study into model estimates provided by
#' joineR::joint() 
#' --------------------------
#' Another simulation study will be done into
#' effects of parameter choice in separate script
#' ##########################

# Prerequisites ------------------------------------------------------------
dev.off()
rm(list = ls())
library(tidyverse)
theme_set(theme_light())
library(lme4)
library(survival)

source("~/Documents/PhD/SMMR_Simulations/simdata.R")

# Set out covariance matrix -----
sigma.0 <- 1   # Intercept RE
sigma.1 <- 0.5   # Slope RE
sigma.2 <- 0.15   # Quadratic RE 

Sigma <- matrix(
  c(sigma.0 ^ 2,       sigma.0 * sigma.1, sigma.0 * sigma.2,
    sigma.1 * sigma.0, sigma.1 ^ 2,       sigma.1 * sigma.2,
    sigma.2 * sigma.0, sigma.2 * sigma.1, sigma.2 ^ 2), nrow = 3, byrow = T
)

# Test run nsims = 5
dat <- suppressMessages(replicate(5, sim.data(Sigma = Sigma), simplify = F))

# Do this in stages -
# Step 1. Cast simulated data to class "jointdata" -----

cast.joint <- function(x){
  long <- left_join(x$long.data, x$surv.data, "id") %>% 
    filter(time <= surv.time)
  
  jd <- jointdata(
    longitudinal = long,
    survival = x$surv.data,
    baseline = x$surv.data[, c("id", "x1", "x3")],
    id.col = "id", time.col = "time"
  )
  
  
  return(jd)
}

dat.joint <- lapply(dat, cast.joint)


# Step 2. Fitting the joint model -----------------------------------------
joint.fit <- function(x){
  joint.fit <- joint(x, model = "quad", max.it = 300,
        long.formula = Yl ~ x1l + x2l + x3l + time,
        surv.formula = Surv(surv.time, status) ~ x1 + x3,
        sepassoc = T)
  return(joint.fit)
}

fits <- lapply(dat.joint, joint.fit)


# Step 3. Extracting Parameter Estimates ----------------------------------
extract.coefs <- function(fit){
  convergence <- fit$convergence
  sigma.e <- sqrt(fit$sigma.z) # sigma.e
  sigma.u <- sqrt(diag(fit$sigma.u)) # sigma.U
  gamma <- fit$coefficients$latent # Latent association
  long.coefs <- t(fit$coefficients$fixed$longitudinal) # Longitudinal coefs
  surv.coefs <- fit$coefficients$fixed$survival # Survival coefs
  
  # All beta (model coefficient) estimates
  beta.ests <- cbind(long.coefs,t(surv.coefs))
  
  # Combine all parameter estimates into one data frame
  param.ests <- data.frame(beta.ests, sigma.e, t(sigma.u), t(gamma), convergence)
  # And return
  return(param.ests)
}

params <- lapply(fits, extract.coefs)
params <- bind_rows(params)

# Step 4. Plot! -----------------------------------------------------------
# Transform parameter set into more workable (long) format
tparams <- params %>% select(-time) %>% 
  gather(-convergence, key = "Parameter", value = "Estimate") 

# Get true values
true.values <- tparams %>% distinct(Parameter)
true.values$xint <- c(40, -10, 5, 15, 0.1, -0.1, 0.01, 1.5, 1, 0.5, 0.15, 1, 0.5, 0.1)
true.values

# Now go about plotting the rest from tparams
ex <- expression
tparams2 <- tparams %>% 
  mutate(
    Param = factor(Parameter,
                   levels = true.values$Parameter,
                   labels = c(ex(beta[0]), ex(beta[1*"L"]), ex(beta[22]), ex(beta[23]),
                              ex(beta[3*"L"]), ex(beta[1*"S"]), ex(beta[3*"S"]),
                              ex(sigma[epsilon]), ex(sigma[U[0]]), ex(sigma[U[1]]),
                              ex(sigma[U[2]]), ex(gamma[0]), ex(gamma[1]), ex(gamma[2]))
    )
  )

true.values$Param <- tparams2 %>% distinct(Param) %>% pull(Param)
true.values <- true.values %>% select(-Parameter)

tparams2 %>% 
  filter(convergence) %>% 
  ggplot(aes(x = Estimate)) + 
  geom_density(fill = "grey20", alpha = .2) + 
  geom_vline(data = true.values, aes(xintercept = xint), lty = 2, colour = "magenta") + 
  facet_wrap(~Param, scales = "free", labeller = label_parsed,
             nrow = 5, ncol = 3) + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black")
  )
  
