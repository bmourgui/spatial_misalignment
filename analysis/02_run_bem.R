#######################################%
#
# Script that fits bem
# to simulated data.
#
# Depends on 01_simulate_data.R 
# which creates the data
#
######################################%

# Load simulate data
#load(here::here("results", "simu", "simulated_data.RData"))

# Create a function to fit bem using nimble
run_bem <- function(data, xname, seed = 409) {
  library(nimble)
  
  bem.nimb <- nimbleCode({
    
    #likelihood
    for (i in 1:n){
      x[i] ~ dnorm(w[i], sd = sd.x)
      y[i] ~ dbern(psi[i])
      logit(psi[i]) <- b0 + b1 * x[i] + b2 * pow(x[i], 2)
    }
    
    #priors
    sd.x ~ T(dt(mu = 0, sigma = 5, df = 1), 0, Inf) # Half-Cauchy with large scale values (sigma) following Gelman et. (2006)
    b0 <- logit(mean.psi)
    mean.psi ~ dunif(0, 1)
    b1 ~ dnorm(0, tau = 0.1)
    b2 ~ dnorm(0, tau = 0.1)
    
  })
  
  myconst <- list(n = nrow(data))
  parameters.to.save <- c("b0", "b1", "b2", "sd.x")
  set.seed(seed)
  initial.values <- list(x = data$x_c1,
                         mean.psi = runif(1, 0, 1), 
                         b1 = rnorm(1, 0, 0.1),
                         b2 = rnorm(1, -1, 0.1),  
                         sd.x = rnorm(1, .5, 0.05))
  
  myd <- list(y = data$y_fine, w = data[, xname])
  
  mymod <- nimbleModel(code = bem.nimb,
                       data = myd,
                       constants = myconst,
                       inits = initial.values)
  
  Cmymod <- compileNimble(mymod)
  
  confMCMC <- configureMCMC(Cmymod, monitors = parameters.to.save)
  myMCMC <- buildMCMC(confMCMC)
  CmyMCMC <- compileNimble(myMCMC, project = mymod)
  
  runMCMC(CmyMCMC,
          nchains = 3,
          nburnin = 10000,
          niter = 210000,
          thin = 50,
          setSeed = seed)
}

# Execute the function in parallel for the three coarsening scenarios
l.bem.c1 <- l.bem.c2 <- l.bem.c3 <- NULL

for (i in seq(5, 30, 5)){
  l.bem.c1[(i-4):i] <- parallel::mclapply(train[(i-4):i], 
                                          function(x)run_bem(data = x, xname = "x_c1"),
                                          mc.cores = 5)
  l.bem.c2[(i-4):i] <- parallel::mclapply(train[(i-4):i], 
                                          function(x)run_bem(data = x, xname = "x_c2"),
                                          mc.cores = 5)
  l.bem.c3[(i-4):i] <- parallel::mclapply(train[(i-4):i], 
                                          function(x)run_bem(data = x, xname = "x_c3"),
                                          mc.cores = 5)
}


# Save outputs
save(l.bem.c1, l.bem.c2, l.bem.c3, 
     file = here::here("results", "simu", "output_bem.RData"))