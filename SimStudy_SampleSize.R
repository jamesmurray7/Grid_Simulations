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
library(joineR)

source("~/Documents/SMMR Simulations _ LINUX/jointsim.R")

# Generate parameter permutations on N ----
simstudy <- crossing(
  m = c(250, 500),
  n = c(5, 10, 15),
) %>% 
  mutate(id = glue::glue("m = {m}, n = {n}")) %>% 
  group_by(id) %>% 
  nest(data = c(m, n)) %>%  ungroup

Sigma <- SigmaGen(1, 0.5, 0.15)

# Simulate data for each permutation of num_subj and num_times

sim200 <- function(x){
  replicate(200, sim.data(num_subj = x$m, num_times = x$n,
                          Sigma = Sigma), simplify = F)
}

simstudy.data <- lapply(simstudy$data, sim200)

# Do this in stages -
# Step 1. Cast simulated data to class "jointdata" -----

dat.joint <- lapply(simstudy.data %>% flatten, cast.joint)


# Step 2. Fitting the joint model -----------------------------------------
library(parallel)
fits <- mclapply(dat.joint, joint.fit,
                 mc.preschedule = T, mc.cores = 4)


# Step 3. Extracting Parameter Estimates ----------------------------------

params <- lapply(fits, extract.params)
params <- bind_rows(params)

id <- simstudy %>% pull(id) %>% as.character() %>% rep(., each =200)
params$id <- id

# Step 4. Plot! -----------------------------------------------------------
# Transform parameter set into more workable (long) format
tparams <- params %>% select(-time) %>% 
  gather(-id, -convergence, key = "Parameter", value = "Estimate") 

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
  ggplot(aes(x = Estimate, fill = id)) + 
  geom_density(alpha = .2) + 
  geom_vline(data = true.values, aes(xintercept = xint), lty = 2, colour = "magenta") + 
  facet_wrap(~Param, scales = "free", labeller = label_parsed,
             nrow = 5, ncol = 3) + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black")
  ) + 
  scale_fill_manual(values = c("grey20", "cyan", "orange", "blue",
                               "red", "darkgreen")) + 
  labs(caption = nrow(tparams[tparams$convergence,])/nrow(tparams) * 100)

ggsave("~/Documents/SMMR Simulations _ LINUX/SampleSizeRedux.png")

params %>% 
  filter(convergence) %>% 
  group_by(id) %>% 
  summarise(across(`X.Intercept.`:gamma_2, median))
save(params, file="~/Documents/SMMR Simulations _ LINUX/SampleSizeRedux.RData")
