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


