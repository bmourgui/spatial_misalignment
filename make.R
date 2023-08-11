######
# Script of modelling workflow used for:
# 'Time for addressing coarse-to-fine spatial misalignment in SDMs'
# B. Mourguiart, M. Chevalier, M. Marzloff, N. Caill-Milly, K. Mengersen, B. Liquet
#
#




# 1. Simulate data
simulate <- TRUE
if (simulate == TRUE){
  source(here::here("analysis", "simu", "simu01_create_rasterize_occu_data.R"))
}else{
  load(here::here("results", "simu", "simulated_data.RData"))
}

# 2. Fit data
## 2.1. with bem
run.bem <- TRUE
if (run.bem == TRUE){
  source(here::here("analysis", "simu", "simu03_nimble_bem.R"))
}else{
  #load(here::here("results", "simu", "output_bem_rep1a15.RData"))
}

## 2.2. with glm
run.glm <- TRUE
if (run.glm == TRUE){
  source(here::here("analysis", "simu", "simu03_nimble_glm.R"))
}else{
  #load(here::here("results", "simu", "output_glm_rep1a15.RData"))
}

## 2.1. with spglm
run.spglm <- TRUE
if (run.spglm == TRUE){
  source(here::here("analysis", "simu", "simu03_run_spglm.R"))
}else{
  #load(here::here("results", "simu", "output_spglm+pred.RData"))
}

# 3. Compute performance metrics
run.perf <- TRUE
if (run.perf == TRUE){
  source(here::here("analysis", "simu", "simu04_compute_perfmetrics.R"))
}else{
  load(here::here("results", "simu", "summary_perf.RData"))
}

# 4. Produce results
source(here::here("analysis", "simu", "simu05_analyze_perfmetrics.R"))