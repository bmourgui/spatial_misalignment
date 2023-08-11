######
# Script of modelling workflow used for:
# "Dealing with grain mismatch between response and explanatory variables in species distribution models."
######

dir.create(here::here("results"))
dir.create(here::here("results/img"))
dir.create(here::here("results/table"))
dir.create(here::here("results/simu"))

# 1. Simulate data
simulate <- TRUE
if (simulate == TRUE){
  source(here::here("analysis", "01_simulate_data.R"))
}else{
  load(here::here("results", "simu", "simulated_data.RData"))
}

# 2. Fit data
## 2.1. with bem
run.bem <- TRUE
if (run.bem == TRUE){
  source(here::here("analysis", "02_run_bem.R")) # can be long
}else{
  load(here::here("results", "simu", "output_bem.RData"))
}

## 2.2. with glm
run.glm <- TRUE
if (run.glm == TRUE){
  source(here::here("analysis", "02_run_glm.R")) # can be long
}else{
  load(here::here("results", "simu", "output_glm.RData"))
}

## 2.1. with spglm
run.spglm <- TRUE
if (run.spglm == TRUE){
  source(here::here("analysis", "02_run_spglm.R")) # can be long
}else{
  load(here::here("results", "simu", "output_spglm+pred.RData"))
}

# 3. Compute performance metrics
run.perf <- TRUE
if (run.perf == TRUE){
  source(here::here("analysis", "03_compute_perfmetrics.R")) # can be long
}else{
  load(here::here("results", "simu", "summary_perf.RData"))
}

# 4. Produce summary results
source(here::here("analysis", "simu", "04_analyze_perfmetrics.R"))

# 5. Plot results
source(here::here("analysis", "simu", "05_make_outputs.R"))