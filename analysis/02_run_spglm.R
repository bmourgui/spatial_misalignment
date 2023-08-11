#######################################%
#
# Script that fits spGLM
# to simulated data using INLA.
#
# Depends on 01_simulate_data.R 
# which creates the data
#
######################################%

## Load simulated data
#load(here::here("results", "simu", "simulated_data.RData"))

library(INLA)
run_spglm <- function(df = train[[2]],
                      df.inter = test[[2]],
                      df.extra1 = pred.new1[[2]],
                      df.extra2 = pred.new2[[2]],
                      y.name = "y_fine",
                      x.name = "x_c2",
                      lon.name = "lon",
                      lat.name = "lat",
                      max.edge = c(10, 50), cutoff = 5, # chose max.edge < spatial range following recommendation in http://rstudio-pubs-static.s3.amazonaws.com/302639_3ccd091c277d4d6c9ad9c2cf524250b6.html
                      prior.sigma = c(5, 0.05), prior.range = c(1, 0.05), 
                      alpha = 2) {
  
  ## Construct mesh, spde, projector and covariate matrices, and stack for model fitting
  # extract fitting variables
  y.fit <- df[, y.name]
  x.fit <- df[, x.name]
  loc.fit <- as.matrix(df[, c(lon.name, lat.name)])
  
  # make a mesh from coordinates of all sites (fit, test, pred)
  mesh <- INLA::inla.mesh.2d(loc = loc.fit,
                             max.edge = max.edge,
                             cutoff = cutoff
  )
  # spde and u.index
  spde <- INLA::inla.spde2.pcmatern(mesh,
                                  alpha = alpha, # alpha = nu + d/2, where d is the dimension of the spatial domain (in 2D, d = 2)
                                  prior.range = prior.range,
                                  prior.sigma = prior.sigma) 
  u.index <- INLA::inla.spde.make.index(name = 'u',
                                        n.spde = spde$n.spde,
                                        n.group = 1, # refers to spatio-temporal groups
                                        n.repl = 1)
  
  # projector matrix
  A.fit <- INLA::inla.spde.make.A(mesh, loc = loc.fit)
  
  # covariate matrix
  Xm.fit <- model.matrix(~ -1 + x + I(x ^ 2),
                         data = data.frame("x" = x.fit))
  X.fit <- data.frame(x = Xm.fit[, 1],
                      x2 = Xm.fit[, 2])
  
  # stack
  stk.fitAG <- INLA::inla.stack(tag = "fitAG", 
                              data = list(y = y.fit), 
                              A = list(1, 1, A.fit), 
                              effects = list( 
                                Intercept = rep(1, length(x.fit)),
                                X = X.fit, 
                                u = u.index 
                              ))
  
  ## Construct projector and covariate matrices, and stack 
  ## for model prediction at fitting points
  ## with covariates at ecological grain
  # extract fitting variables
  x.fitEG <- df[, "x_fine"]
  
  # projector matrix same as fit
  # covariate matrix
  Xm.fitEG <- model.matrix(~ -1 + x + I(x ^ 2),
                             data = data.frame("x" = x.fitEG))
  X.fitEG <- data.frame(x = Xm.fitEG[, 1],
                          x2 = Xm.fitEG[, 2])
  
  # stack
  stk.fitEG <- INLA::inla.stack(tag = "fitEG", 
                                  data = list(y = NA), 
                                  A = list(1, 1, A.fit), 
                                  effects = list( 
                                    Intercept = rep(1, length(x.fitEG)),
                                    X = X.fitEG, 
                                    u = u.index 
                                  ))
  
  ## Construct projector and covariate matrices, and stack 
  ## for model interpolation
  ## with covariates at analysis grain
  # extract fitting variables
  x.interAG <- df.inter[, x.name]
  loc.inter <- as.matrix(df.inter[, c(lon.name, lat.name)])
  
  # projector matrix
  A.interAG <- INLA::inla.spde.make.A(mesh, 
                                      loc = loc.inter)
  
  # covariate matrix
  Xm.interAG <- model.matrix(~ -1 + x + I(x ^ 2),
                         data = data.frame("x" = x.interAG))
  X.interAG <- data.frame(x = Xm.interAG[, 1],
                      x2 = Xm.interAG[, 2])
  
  # stack
  stk.interAG <- INLA::inla.stack(tag = "interAG", 
                              data = list(y = NA), 
                              A = list(1, 1, A.interAG), 
                              effects = list( 
                                Intercept = rep(1, length(x.interAG)),
                                X = X.interAG, 
                                u = u.index 
                              ))
  
  ## with covariates at ecological grain
  # extract fitting variables
  x.interEG <- df.inter[, "x_fine"]
  
  # projector matrix
  A.interEG <- INLA::inla.spde.make.A(mesh, 
                                      loc = loc.inter)
  
  # covariate matrix
  Xm.interEG <- model.matrix(~ -1 + x + I(x ^ 2),
                             data = data.frame("x" = x.interEG))
  X.interEG <- data.frame(x = Xm.interEG[, 1],
                          x2 = Xm.interEG[, 2])
  
  # stack
  stk.interEG <- INLA::inla.stack(tag = "interEG", 
                                  data = list(y = NA), 
                                  A = list(1, 1, A.interEG), 
                                  effects = list( 
                                    Intercept = rep(1, length(x.interEG)),
                                    X = X.interEG, 
                                    u = u.index 
                                  ))
  
  
  ## Construct projector and covariate matrices, and stack 
  ## for model extrapolation in smoother environment
  ## with covariates at analysis grain
  
  # extract fitting variables
  x.extra1AG <- df.extra1[, x.name]
  loc.extra1AG <- as.matrix(df.extra1[, c(lon.name, lat.name)])
  
  # projector matrix
  A.extra1AG <- INLA::inla.spde.make.A(mesh, 
                                      loc = loc.extra1AG)
  
  # covariate matrix
  Xm.extra1AG <- model.matrix(~ -1 + x + I(x ^ 2),
                             data = data.frame("x" = x.extra1AG))
  X.extra1AG <- data.frame(x = Xm.extra1AG[, 1],
                          x2 = Xm.extra1AG[, 2])
  
  # stack
  stk.extra1AG <- INLA::inla.stack(tag = "extra1AG", 
                                  data = list(y = NA), 
                                  A = list(1, 1, A.extra1AG), 
                                  effects = list( 
                                    Intercept = rep(1, length(x.extra1AG)),
                                    X = X.extra1AG, 
                                    u = u.index 
                                  ))
  
  ## with covariates at ecological grain
  # extract fitting variables
  x.extra1EG <- df.extra1[, "x_fine"]
  
  # projector matrix
  A.extra1EG <- INLA::inla.spde.make.A(mesh, 
                                      loc = loc.extra1AG)
  
  # covariate matrix
  Xm.extra1EG <- model.matrix(~ -1 + x + I(x ^ 2),
                             data = data.frame("x" = x.extra1EG))
  X.extra1EG <- data.frame(x = Xm.extra1EG[, 1],
                          x2 = Xm.extra1EG[, 2])
  
  # stack
  stk.extra1EG <- INLA::inla.stack(tag = "extra1EG", 
                                  data = list(y = NA), 
                                  A = list(1, 1, A.extra1EG), 
                                  effects = list( 
                                    Intercept = rep(1, length(x.extra1EG)),
                                    X = X.extra1EG, 
                                    u = u.index 
                                  ))
  
  
  ## Construct projector and covariate matrices, and stack 
  ## for model extrapolation in more heterogeneous environment
  ## with covariates at analysis grain
  
  # extract fitting variables
  x.extra2AG <- df.extra2[, x.name]
  loc.extra2AG <- as.matrix(df.extra2[, c(lon.name, lat.name)])
  
  # projector matrix
  A.extra2AG <- INLA::inla.spde.make.A(mesh, 
                                       loc = loc.extra2AG)
  
  # covariate matrix
  Xm.extra2AG <- model.matrix(~ -1 + x + I(x ^ 2),
                              data = data.frame("x" = x.extra2AG))
  X.extra2AG <- data.frame(x = Xm.extra2AG[, 1],
                           x2 = Xm.extra2AG[, 2])
  
  # stack
  stk.extra2AG <- INLA::inla.stack(tag = "extra2AG", 
                                   data = list(y = NA), 
                                   A = list(1, 1, A.extra2AG), 
                                   effects = list( 
                                     Intercept = rep(1, length(x.extra2AG)),
                                     X = X.extra2AG, 
                                     u = u.index 
                                   ))
  
  ## with covariates at ecological grain
  # extract fitting variables
  x.extra2EG <- df.extra2[, "x_fine"]
  
  # projector matrix
  A.extra2EG <- INLA::inla.spde.make.A(mesh, 
                                       loc = loc.extra2AG)
  
  # covariate matrix
  Xm.extra2EG <- model.matrix(~ -1 + x + I(x ^ 2),
                              data = data.frame("x" = x.extra2EG))
  X.extra2EG <- data.frame(x = Xm.extra2EG[, 1],
                           x2 = Xm.extra2EG[, 2])
  
  # stack
  stk.extra2EG <- INLA::inla.stack(tag = "extra2EG", 
                                   data = list(y = NA), 
                                   A = list(1, 1, A.extra2EG), 
                                   effects = list( 
                                     Intercept = rep(1, length(x.extra2EG)),
                                     X = X.extra2EG, 
                                     u = u.index 
                                   ))
  
  stk.all <- INLA::inla.stack(stk.fitAG,
                              stk.fitEG,
                              stk.interAG,
                              stk.interEG,
                              stk.extra1AG,
                              stk.extra1EG,
                              stk.extra2AG,
                              stk.extra2EG)
  f <- y ~ -1 + Intercept + x + x2 + f(u, model = spde)
  #cmp <- y_fine ~ Intercept(1) + x(get(x.name)) + x2(get(x.name) ^ 2) + u(coordinates, model = spde)
  #out <- inlabru::bru(cmp, data, family = "binomial", Ntrials = 1)
  
  out <- INLA::inla(f,
                    family = "binomial", Ntrials = 1, # Bernoulli
                    data = INLA::inla.stack.data(stk.all),
                    #control.compute = list(
                    #  config = TRUE, # save the internal GMRF approximations in the created object
                    #  dic = TRUE,
                    #  waic = TRUE
                    #),
                    control.predictor = list(
                      A = INLA::inla.stack.A(stk.all),
                      compute = TRUE, # compute marginals for the linear predictor
                      link = 1
                    )
  )
  spfi <- INLA::inla.spde2.result(inla = out,
                                  name = "u",
                                  spde = spde,
                                  do.transform = TRUE)
  out[[(length(out) + 1)]] <- spfi
  names(out)[length(out)] <- "spatial.random.field"
  return(out)
  
}


# Execute the function in parallel for the three coarsening scenarios
l.spglm.c1 <- l.spglm.c2 <- l.spglm.c3 <- NULL
for (i in 1:30){
  l.spglm.c1[[i]] <- run_spglm(df = train[[i]],
                              df.inter = test[[i]],
                              df.extra1 = pred.new1[[i]],
                              df.extra2 = pred.new2[[i]],
                              x.name = "x_c1")
  
  l.spglm.c2[[i]] <- run_spglm(df = train[[i]],
                               df.inter = test[[i]],
                               df.extra1 = pred.new1[[i]],
                               df.extra2 = pred.new2[[i]],
                               x.name = "x_c2") 
  
  l.spglm.c3[[i]] <- run_spglm(df = train[[i]],
                               df.inter = test[[i]],
                               df.extra1 = pred.new1[[i]],
                               df.extra2 = pred.new2[[i]],
                               x.name = "x_c3") 
}
save(l.spglm.c1, l.spglm.c2, l.spglm.c3,
     file = here::here("results", "simu", "output_spglm+pred.RData"))
