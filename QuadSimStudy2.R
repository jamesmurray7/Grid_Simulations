#' ##########################
#' QuadSimStudy2.R   Author: James Murray
#' --------------------------
#' Simulation study into effect of changing parameter values
#' on estimations provided by joineR::joint()
#' --------------------------
#' ##########################

# Prerequisites ------------------------------------------------------------
dev.off()
rm(list = ls())
library(tidyverse)
theme_set(theme_light())
library(lme4)
library(survival)

source("~/Documents/PhD/SMMR_Simulations/jointsim.R")

# Define function to generate Covariance matrix -----
SigmaGen <- function(sigma.0, sigma.1, sigma.2){
  Sigma <- matrix(
    c(sigma.0 ^ 2,       sigma.0 * sigma.1, sigma.0 * sigma.2,
      sigma.1 * sigma.0, sigma.1 ^ 2,       sigma.1 * sigma.2,
      sigma.2 * sigma.0, sigma.2 * sigma.1, sigma.2 ^ 2), nrow = 3, byrow = T
  )
  return(Sigma)
}

# Test run nsims = 5
nsims <- 5
dat <- suppressMessages(replicate(nsims, sim.data(Sigma = Sigma), simplify = F))

# Do this in stages -
# Step 1. Cast simulated data to class "jointdata" -----

dat.joint <- lapply(dat, cast.joint)


# Step 2. Fitting the joint model -----------------------------------------

fits <- lapply(dat.joint, joint.fit)


# Step 3. Extracting Parameter Estimates ----------------------------------

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
  ) + 
  labs(caption = nrow(tparams[tparams$convergence,])/nrow(tparams) * 100)
ggsave("~/Documents/PhD/SMMR_Simulations/QuadParameterEstimates.png")  
