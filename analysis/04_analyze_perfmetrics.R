#################################################################%
##
## Code to evaluate explanatory and predictive power of models
##
## Depends on 03_compute_perfmetrics.R and its dependencies
## 
#################################################################%



source("analysis/01_simulate_data.R") # /!\ re-run 01_simulate_data.R doesn't work with saved simulated_data
load(file = here::here("results", "simu", "output_bem.RData"))
load(file = here::here("results", "simu", "output_glm.RData"))
load(file = here::here("results", "simu", "output_spglm+pred.RData"))
load(here::here("results", "simu", "summary_perf.RData"))

library(magrittr)
library(ggplot2)

# Main analysis ----
## Estimation of species-environment relationships ----
### Bias in SERs parameters ----
# Compute bias and relative bias in SER estimates (optimum, maximum, width)
d.bias.ser <- summary.src |> 
  dplyr::select(model, grain, repli, type, opti:width) |> 
  tidyr::pivot_longer(opti:width, values_to = "value", names_to = "param") |> 
  tidyr::pivot_wider(names_from = "type", values_from = "value") |> 
  dplyr::mutate("bias" = mean - true,
                "relative.bias" = (mean - true)/true) |> 
  dplyr::filter(param != "tol")

# True values of SER parameters
d.true.ser <- d.bias.ser |> 
  dplyr::group_by(grain, param) |> 
  dplyr::summarise_at("true", mean)

# Compute summary statistics (mean, sd) for SER estimates
d.summary.estim.ser <- d.bias.ser |> 
  dplyr::group_by(grain, param, model) |> 
  dplyr::summarise_at("mean", c("avg" = mean, "sd" = sd))

# Compute summary statistics (mean, sd) for bias in SER estimates
d.summary.bias.ser <- d.bias.ser |> 
  dplyr::group_by(grain, param, model) |> 
  dplyr::summarise_at("bias", c("avg" = mean, "sd" = sd))

# Compute summary statistics (mean, sd) for relative bias in SER estimates
d.summary.relbias.ser <- d.bias.ser |> 
  dplyr::group_by(grain, param, model) |> 
  dplyr::summarise_at("relative.bias", c("avg" = mean, "sd" = sd))


### Explanatory performance metrics ----
# Averaging results for combinations: model x metric x covariate grain x prediction grain (grain from which predictions are made)
d.mean.perf.expla <- summary.psi |>
  dplyr::filter(data.pred == "train") |> 
  tidyr::pivot_longer(rmse:maxtss, names_to = "metric", values_to = "value") |>  
  dplyr::group_by(model, metric, fit.grain, pred.grain, data.pred) |> 
  dplyr::summarise_at("value", c("avg" = mean, "sd" = sd))

## Predictive performance metrics ----
# Averaging results for combinations: model x metric x covariate grain x prediction grain (grain from which predictions are made)
d.mean.perf.pred <- summary.psi |>
  dplyr::filter(data.pred != "train") |> 
  tidyr::pivot_longer(rmse:maxtss, names_to = "metric", values_to = "value") |>
  dplyr::group_by(model, metric, fit.grain, pred.grain, data.pred) |> 
  dplyr::summarise_at("value", c("avg" = mean, "sd" = sd))


# Supporting information ----
## Compute SRCs from coefficient estimates ----
x.pred <- seq(min(data$x_fine) - .5, max(data$x_fine) + .5, length.out = 200) # gradient of environmental values
model <- c("glm", "bem", "spglm")
CG <- c("c1", "c2", "c3")

n.mod <- length(model)
n.CG <- length(CG)
n.repli <- 30
n.x <- length(x.pred)

d.src <- data.frame("model" = rep(model, each = n.CG*n.repli*n.x),
                   "fit.grain" = rep(rep(CG, each = n.repli*n.x), n.mod),
                   "repli" = rep(1:n.repli, each = n.x),
                   "x.pred" = x.pred,
                   "psi.true" = plogis(b0 + b1*x.pred + b2*(x.pred^2)),
                   "psi.pred" = NA)

for (i in model){
  for (j in CG){
    for (k in 1:n.repli){
      xc <- summary.src |> 
        dplyr::filter(model == i & grain == j & type == "mean" & repli == k)
      ppred <- plogis(xc$b0 + xc$b1*x.pred + xc$b2*(x.pred^2))
      d.src <- d.src |> 
        dplyr::mutate("psi.pred" = replace(psi.pred, 
                                           model == i & fit.grain == j & repli == k, 
                                           ppred))
    }
  }
}

## Loss of variability in environmental conditions during coarsening ----

# function calculating fine-grain variance within or between coarse-grain cells
var_loss <- function(xfine, # fine-grain raster
                     grid.crop = grid.sple, 
                     d = seq(10, 500, 10)){ # vector of coarsening factors
  
  xcoarse <- lapply(d,
                    function(l)
                      terra::aggregate(xfine, fact = l, fun = mean)
  )
  
  var_inter.cell <- lapply(xcoarse, 
                           function(x){
                             z <- terra::crop(x, grid.crop)
                             var(terra::values(z), na.rm = F)
                           }
  )
  
  var_intra.cell <- lapply(d,
                           function(l){
                             z <- terra::aggregate(xfine, fact = l, fun = var) |> 
                               terra::crop(grid.crop)
                             mean(terra::values(z), na.rm = F)
                           }
  )
  list("var_inter.cell" = var_inter.cell, "var_intra.cell" = var_intra.cell)
}  

grid.crop <- grid.sple # avoid "edge effect" when coarsening data
d <- seq(2, 200, 1) # distance gradient on which calculate variability loss    

# Calculate variances for the three types of environmental heterogeneity
var_loss.medium <- var_loss(x_fine, grid.crop, d)
var_loss.smooth <- var_loss(x.new1_fine, grid.crop, d)
var_loss.hete <- var_loss(x.new2_fine, grid.crop, d)


# Not used ####
# Supplementary analysis not included in the manuscript or supporting information

## Convergence check for BEM ----
did.converge <- function(out.nimble){
  Rhat <- MCMCvis::MCMCsummary(object = out.nimble)[, "Rhat"]
  sum(Rhat > 1.1) == 0
}

lapply(l.bem.c1, did.converge)
lapply(l.bem.c2, did.converge)
lapply(l.bem.c3, did.converge)

## Bias in berkson error estimate ----

# compare error estimate by BEM with "true" error approximated by
# sd(X^EG - X^CG) because a BE is defined as x.true = x.obs + N(O, sd.error)

# approximate the true error
true.error <- data |> 
  dplyr::mutate("diff.c1" = x_fine - x_c1,
                "diff.c2" = x_fine - x_c2,
                "diff.c3" = x_fine - x_c3) |> 
  dplyr::summarise_at(c("diff.c1", "diff.c2", "diff.c3"), sd) |> 
  as.numeric()
d.true_sd.x <- data.frame("grain" = c("c1", "c2", "c3"),
                          "sd.x" = true.error)

# extract error estimates
mean.sd.bem <- function(out.nimble){
  (MCMCvis::MCMCsummary(object = out.nimble)["sd.x", "mean"])
}

d_sd.x <- data.frame("grain" = rep(c("c1", "c2", "c3"), each = length(l.bem.c1)),
                     "sd.x" = c(mapply(mean.sd.bem, l.bem.c1),
                                mapply(mean.sd.bem, l.bem.c2),
                                mapply(mean.sd.bem, l.bem.c3)))

d.mean_sd.x <- d_sd.x |> 
  dplyr::group_by(grain) |> 
  dplyr::summarise_at("sd.x", .funs = c("avg" = mean, "sd" = sd))


## Estimates of spatial parameters by spGLM ----
extract_param_sp <- function(
    inla = l.spglm.c1[[1]]
){
  
  
  spfi <- inla$spatial.random.field
  r <- INLA::inla.emarginal(function(x){x},
                            spfi$marginals.range.nominal[[1]])
  v <- INLA::inla.emarginal(function(x){x},
                            spfi$marginals.variance.nominal[[1]])
  c("range" = r,
    "variance" = v)
  
}

sp_c1 <- t(mapply(function(x, y)extract_param_sp(inla = y),
                  train[1:15], l.spglm.c1))
sp_c2 <- t(mapply(function(x, y)extract_param_sp(inla = y),
                  train[1:15], l.spglm.c2))
sp_c3 <- t(mapply(function(x, y)extract_param_sp(inla = y),
                  train[1:15], l.spglm.c3))

d_sp.x <- data.frame("grain" = rep(c("c1", "c2", "c3"), each = length(l.spglm.c1)),
                     "range.x" = c(sp_c1[, 1],
                                   sp_c2[, 1],
                                   sp_c3[, 1]),
                     "variance" = c(sp_c1[, 2],
                                    sp_c2[, 2],
                                    sp_c3[, 2])) |> 
  tidyr::pivot_longer(range.x:variance,
                      names_to = "param",
                      values_to = "estimate")

d.mean_sp.x <- d_sp.x |> 
  dplyr::group_by(grain, param) |> 
  dplyr::summarise_at("estimate", .funs = c("avg" = mean, "sd" = sd))

## Bias in models parameters (the betas) ----
d.bias.beta <- summary.src |> 
  dplyr::select(model, grain, repli, type, b0:b2) |> 
  tidyr::pivot_longer(b0:b2, values_to = "value", names_to = "param") |> 
  tidyr::pivot_wider(names_from = "type", values_from = "value")

d.true.beta <- d.bias.beta |> 
  dplyr::group_by(grain, param) |> 
  dplyr::summarise_at("true", mean)

d.summary.bias.beta <- d.bias.beta |> 
  dplyr::group_by(grain, param, model) |> 
  dplyr::summarise_at("mean", c("avg" = mean, "sd" = sd))




