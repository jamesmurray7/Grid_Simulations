#' #############
#' simdata.R    Author: Jame Murray
#' function to simulate data from the two
#' sub-models prior to analysis in a joint model
#' #############

sim.data <- function(num_subj = 250, num_times = 6,
                     sigma.epsilon = 1.5, gamma = c(1, 0.5, 0.1),
                     Sigma, beta.longit = c(40, -10, 5, 15, 0.1),
                     beta.surv = c(-0.1,0.01), censoring = T, cens.rate = 0.01,
                     dt = 0.001, print.data = F){
  
  # Checks on provided covariance matrix ----
  if(!isSymmetric(Sigma)){
    stop("Covariance matrix Sigma not symmetric")
  }
  if(any(eigen(Sigma)$values < 0) || det(Sigma) <= 0){
    message("Provided covariance matrix is not positive definite, transforming...")
    Sigma <- as.matrix(Matrix::nearPD(Sigma)$mat)
  }
  
  N <- num_subj * num_times
  tau <- num_times - 1
  
  # Generate random effects ----
  REs <- MASS::mvrnorm(num_subj, mu = c(0, 0, 0), Sigma)
  U0 <- REs[,1]; U1 <- REs[,2]; U2 <- REs[,3]
  
  # Baseline Covariates ----
  id <- 1:num_subj
  x1 <- rbinom(num_subj, 1, 0.5) # Treatment received
  x2 <- gl(3, 1, num_subj) # Factor
  x3 <- floor(rnorm(num_subj, 65, 10)); x3 <- scale(x3, scale = F) # Continuous
  
  # Longitudinal part ----
  Bl <- matrix(beta.longit, nrow = 1)
  x1l <- rep(x1, each = num_times)
  x2l <- rep(x2, each = num_times)
  x3l <- rep(x3, each = num_times)
  Xl <- model.matrix(~x1l + x2l + x3l)
  time <- rep(0:tau, num_subj)
  Yl <- Xl %*% t(Bl) + rep(U0, each = num_times) + rep(U1, each = num_times) * time + rnorm(N, 0, sigma.epsilon)
  long.data <- data.frame(id = rep(id, each = num_times), time, x1l, x2l, x3l, Yl)
  
  # Survival Part ----
  theta.0 <- log(0.01) # Gompertz baseline hazard: Scale parameter = exp(theta.0)
  theta.1 <- 0.5 # 0.1 # Gompertz baseline hazard: Shape parameter
  # Fixed effects
  Bs <- matrix(beta.surv, nrow = 1)
  Xs <- model.matrix(~ x1 + x3 - 1)
  XsBs <- Xs %*% t(Bs)
  # Define gridsteps
  gridt <- seq(0, tau, dt)[-1] 
  # Piece together Hazard at each time-point t for each subject
  bl.haz <- exp(XsBs) %*% exp(theta.0 + theta.1 * gridt) 
  gamma.U <- exp(gamma[1] * U0 + gamma[2] * U1 * gridt + gamma[3] * U2 * gridt^2)
  lambda <- gamma.U * bl.haz
  lambda.dt <- lambda * dt # These are the hazards (scaled)
  
  # Creating matrices of candidate times and Uniform RVs to generate survival times
  dims <- dim(lambda.dt)
  candidate.times <- matrix(gridt, nrow = dims[1], ncol = dims[2], byrow = T)
  unifs <- matrix(runif(prod(dims)), nrow = dims[1], ncol = dims[2])
  
  candidate.times[lambda.dt < unifs] <- tau
  surv.time <- apply(candidate.times, 1, min)
  
  # Survival - Censoring ----
  
  if(censoring){
    status <- rep(1, num_subj)
    cens.time <- rexp(num_subj, cens.rate)
    f <- cens.time < surv.time # If censoring time predates survival time then
    status[f] <- 0 # Set status = 0
    surv.time[f] <- cens.time[f] # And set survival time = censoring time.
  }else{
    status <- 1
  }
  status[surv.time == tau] <- 0 # If survival time = truncation time then they survive.
  
  surv.data <- data.frame(id, x1, x3, surv.time, status)
  
  if(print.data){
    message("Longitudinal data:")
    print(long.data[1:25,])
    message("Survival data:")
    print(surv.data[1:25,])
  }
  message("Failure rate of ", round(length(which(status == 1))/num_subj * 100), "%")
  
  return(list(long.data = long.data, surv.data = surv.data))  
}