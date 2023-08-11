#######################################%
#
# Script that fits GLM
# to simulated data.
#
# Depends on 01_simulate_data.R 
# which creates the data
#
######################################%

# Load simulate data
#load(here::here("results", "simu", "simulated_data.RData"))

# Create a function to fit glm using nimble
run_glm <- function(data, xname, seed = 409) {
  library(nimble)
  
  glm.nimb <- nimbleCode({
    
    #Likelihood
    for (i in 1:n){
      y[i] ~ dbern(psi[i])
      logit(psi[i]) <- b0 +  b1 * w[i] + b2 * pow(w[i], 2)
    }
    
    #Priors
    b0 <- logit(mean.psi)
    mean.psi ~ dunif(0, 1)
    b1 ~ dnorm(0, 0.1)
    b2 ~ dnorm(0, 0.1)
    
  })
  
  myconst <- list(n = nrow(data))
  parameters.to.save <- c("b0", "b1", "b2")
  set.seed(seed)
  initial.values <- list(mean.psi = rnorm(1, 0.5, 0.1), 
                         b1 = rnorm(1, 0, 0.1),
                         b2 = rnorm(1, -1, 0.1))
  
  myd <- list(y = data$y_fine, w = data[, xname])
  
  mymod <- nimbleModel(code = glm.nimb,
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
l.glm.c1 <- l.glm.c2 <- l.glm.c3 <- NULL

for (i in seq(5, 30, 5)){
  l.glm.c1[(i-4):i] <- parallel::mclapply(train[(i-4):i], 
                                          function(x)run_glm(data = x, xname = "x_c1"),
                                          mc.cores = 5)
  l.glm.c2[(i-4):i] <- parallel::mclapply(train[(i-4):i], 
                                          function(x)run_glm(data = x, xname = "x_c2"),
                                          mc.cores = 5)
  l.glm.c3[(i-4):i] <- parallel::mclapply(train[(i-4):i], 
                                          function(x)run_glm(data = x, xname = "x_c3"),
                                          mc.cores = 5)
}
# Save outputs
save(l.glm.c1, l.glm.c2, l.glm.c3, 
     file = here::here("results", "simu", "output_glm.RData"))
