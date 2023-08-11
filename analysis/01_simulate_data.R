#################################################################%
##
## Code to simulate the data used in misalignment paper
##
#################################################################%

# 0. Install the R-package RandomFields ----

# Steps to install R-package RandomFields latest version.
## 0.1. install gfortran 
## 0.2. localise gfortran in R
  ## dir.create('~/.R')
  ## file.create('~/.R/Makevars')
  ## in this file copy-paste the followings (without #): 
  #FC      = /usr/local/gfortran/bin/gfortran # this path might be changed depending on the localisation of gfortran on your computer
  #F77     = /usr/local/gfortran/bin/gfortran # this path might be changed depending on the localisation of gfortran on your computer
  #FLIBS   = -L/usr/local/opt/gcc/lib
## 0.3. restart R 

## 0.4. install packages
  ## install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/RandomFieldsUtils_1.2.5.tar.gz", type="source")
  ## install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFields/RandomFields_3.3.14.tar.gz", type = "source")

# 1. Simulate environment and species data at ecological grain ----
## 1.1. Simulate a Gaussian random field for the environment ----
# Create the raster
# Simulate a larger grid to the one surveyed to avoid scenarios where coarsening the data induce NAs at boarders
x_fine <- terra::rast(nrow = 1500, ncols = 1500, xmin = -250, xmax = 1250, ymin = -250, ymax = 1250)
grid <- terra::xyFromCell(x_fine, 1:terra::ncell(x_fine)) # simulated grid
grid.sple <- cbind(rep(seq(0.5, 1000, 1), 1000),
                   rep(seq(0.5, 1000, 1), each = 1000)) # sampled grid

# Define parameters of the Matern covariance function
variance <- 1
nu <- 1
range <- 50
kappa <- sqrt(8*nu) / range
  # relation between the range and kappa is from Chapt 6, p.194 Blangiardo & Cameletti (2015)
  # range = sqrt(8*nu)/kappa
  # the range is the distance for which the correlation function has fallen to approximately 0.13.
scale <- sqrt(2 * nu) / kappa
  # this relation denotes the difference in parametrization between the package RandomFields and INLA
  # in RandomFields::RMmatern() covariance function is defined by : Cov.f(dist) = var * BesselK(sqrt(2 * nu) * dist / scale ) * (sqrt(2 * nu) * dist / scale )^nu * 2^(1 - nu) / Gamma(nu)
  # while in INLA by: Cov.f(dist) = var * BesselK(dist * kappa ) * (dist * kappa )^nu * 2^(1 - nu) / Gamma(nu)
  # thus: scale = sqrt(2 * nu) / kappa

# Simulate the covariance function and spatial random field
RandomFields::RFoptions(seed = 321,
                        install = "no")
matern <- RandomFields::RMmatern(nu = nu, var = variance, scale = scale)
# Attribute values to the raster according to the spatial random field
terra::values(x_fine) <- RandomFields::RFsimulate(model = matern, 
                                                  x = grid)@data$variable1

#####%
# Alternative simulation with inlabru following https://inlabru-org.github.io/inlabru/articles/web/random_fields_2d.html
# (takes long time to run)
#library(inlabru)
#bnd <- spoly(data.frame(easting = c(-250, 1250, 1250, -250), northing = c(-250, -250, 1250, 1250)))
#mesh_fine <- INLA::inla.mesh.2d(boundary = bnd, max.edge = 2)
#matern_fine <- INLA::inla.spde2.pcmatern(mesh_fine,
#                                   # prior not used to simulate data but needed in the function anyway
#                                   prior.sigma = c(10, 0.01), 
#                                   prior.range = c(1, 0.01) 
#)

#Q <- inla.spde.precision(matern_fine, theta = log(c(range, sqrt(variance))))
#x.field <- inla.qsample(1, Q.u)[, 1]

#terra::values(x_fine) <- fm_evaluate(
#  mesh_fine,
#  loc = grid,
#  field = x.field)

#####%


## 1.2. Simulate species response with ecological fine grain environment ----
## Define coefficients of the species-environment relationship
logit <- function(x)log(x / (1 - x))

opti <- 0.1
psi.max <- 0.95
width <- 2.5

a <- logit(psi.max)
tau <- 0.5 * width / sqrt(2 * (a - logit(.05))) # see appendix 1 in Mourguiart et al. (2023)

b0 <- 3  #a - (opti^2)/(2*tau^2)
b1 <- 1.5 #2*opti/(tau ^ 2)
b2 <- -3 #-1/(2*tau^2)

data <- data.frame("lon" = grid[, 1],
                   "lat" = grid[, 2],
                   "x_fine" = as.numeric(terra::values(x_fine))) |> 
  dplyr::filter(lon %in% grid.sple[, 1] & lat %in% grid.sple[, 2]) |>  
  dplyr::rowwise() |> 
  dplyr::mutate("psi_fine" = plogis(b0 + b1 * x_fine + b2 * (x_fine^2)),
                "y_fine" = rbinom(1, 1, psi_fine))

# 2. Simulate sampling process ----
# Define values of grain- and sample-size of the different scenarios
grain <- c(range/10, range/2, range)
resolution <- grain
sple.size <- 300

# Coarsening environmental covariate for each grain-size scenario
x_coarse <- lapply(resolution,
                   function(l)
                     terra::aggregate(x_fine, fact = l, fun = mean)
                   )

# add coarse covariates in the data
x_c <- lapply(x_coarse,
              function(v){
                terra::extract(v, 
                               data[, c("lon", "lat")], 
                               ID = FALSE)[, 1]
              }
              )
names(x_c) <- paste0("x_c", 1:length(x_c))
data <- cbind(data, do.call(cbind, x_c))

# Define number of replications
n.repli <- 30

# Simulate sampling
n.samp <- 300
n.train <- 300
n.test <- 300

set.seed(321) 
sample <- lapply(rep((n.test + n.train), n.repli), 
                 function(n)data |> dplyr::slice_sample(n = n)
                 )
train <- lapply(sample,
                function(s)s[1:n.train, ])
test <- lapply(sample,
               function(s)s[-c(1:n.train), ])

# 3. Simulate new environments for prediction ----
## 3.1. Extrapolation data 1, smoother environment ----
x.new1_fine <- terra::rast(nrow = 1500, ncols = 1500, xmin = -250, xmax = 1250, ymin = -250, ymax = 1250)

range_x.new1 <- range*5
kappa_x.new1 <- sqrt(8*nu) / range_x.new1
scale_x.new1 <- sqrt(2 * nu) / kappa_x.new1

RandomFields::RFoptions(seed = 321,
                        install = "no")
matern_x.new1 <- RandomFields::RMmatern(nu = nu, var = variance, scale = scale_x.new1)
# Attribute values to the raster according to the spatial random field
terra::values(x.new1_fine) <- RandomFields::RFsimulate(model = matern_x.new1, 
                                                    x = grid)@data$variable1

#####%
# Simulate smooth field with inlabru

#Q.new1 <- inla.spde.precision(matern_fine, theta = log(c(range_x.new1, sqrt(variance))))
#x.new1 <- inla.qsample(1, Q.new1)[, 1]

#terra::values(x.new1_fine) <- fm_evaluate(
#  mesh_fine,
#  loc = grid,
#  field = x.new1)

#####%

d.new1 <- data.frame("lon" = grid[, 1],
                          "lat" = grid[, 2],
                          "x_fine" = as.numeric(terra::values(x.new1_fine))) |> 
  dplyr::filter(lon %in% grid.sple[, 1] & lat %in% grid.sple[, 2]) |>  
  dplyr::rowwise() |> 
  dplyr::mutate("psi_fine" = plogis(b0 + b1 * x_fine + b2 * (x_fine^2)),
                "y_fine" = rbinom(1, 1, psi_fine))


x.new1_coarse <- lapply(resolution,
                   function(l)
                     terra::aggregate(x.new1_fine, fact = l, fun = mean)
)

x.new1_c <- lapply(x.new1_coarse,
              function(v){
                terra::extract(v, 
                               d.new1[, c("lon", "lat")], 
                               ID = FALSE)[, 1]
              }
)
names(x.new1_c) <- paste0("x_c", 1:length(x.new1_c))
d.new1 <- cbind(d.new1, do.call(cbind, x.new1_c))

# Random samples 
pred.new1 <- lapply(rep(n.samp, n.repli), 
               function(n)d.new1 |> dplyr::slice_sample(n = n)
)


## 3.2. Extrapolation data 2, more heterogeneous environment ----
x.new2_fine <- terra::rast(nrow = 1500, ncols = 1500, xmin = -250, xmax = 1250, ymin = -250, ymax = 1250)

range_x.new2 <- range/5
kappa_x.new2 <- sqrt(8*nu) / range_x.new2
scale_x.new2 <- sqrt(2 * nu) / kappa_x.new2

RandomFields::RFoptions(seed = 321,
                        install = "no")
matern_x.new2 <- RandomFields::RMmatern(nu = nu, var = variance, scale = scale_x.new2)
# Attribute values to the raster according to the spatial random field
terra::values(x.new2_fine) <- RandomFields::RFsimulate(model = matern_x.new2, 
                                                    x = grid)@data$variable1

#####%
# Simulate heterogeneous field with inlabru

#Q.new2 <- inla.spde.precision(matern_fine, theta = log(c(range_x.new2, sqrt(variance))))
#x.new2 <- inla.qsample(1, Q.new2)[, 1]

#terra::values(x.new2_fine) <- fm_evaluate(
#  mesh_fine,
#  loc = grid,
#  field = x.new2)

#####%

d.new2 <- data.frame("lon" = grid[, 1],
                       "lat" = grid[, 2],
                       "x_fine" = as.numeric(terra::values(x.new2_fine))) |> 
  dplyr::filter(lon %in% grid.sple[, 1] & lat %in% grid.sple[, 2]) |>  
  dplyr::rowwise() |> 
  dplyr::mutate("psi_fine" = plogis(b0 + b1 * x_fine + b2 * (x_fine^2)),
                "y_fine" = rbinom(1, 1, psi_fine))


x.new2_coarse <- lapply(resolution,
                          function(l)
                            terra::aggregate(x.new2_fine, fact = l, fun = mean)
)

x.new2_c <- lapply(x.new2_coarse,
                     function(v){
                       terra::extract(v, 
                                      d.new2[, c("lon", "lat")], 
                                      ID = FALSE)[, 1]
                     }
)
names(x.new2_c) <- paste0("x_c", 1:length(x.new2_c))
d.new2 <- cbind(d.new2, do.call(cbind, x.new2_c))

# Random samples
pred.new2 <- lapply(rep(n.samp, n.repli), 
                      function(n)d.new2 |> dplyr::slice_sample(n = n)
)

# 4. Save simulated data ----
save(d.new2, d.new1, data, grid, grid.sple, pred.new2, pred.new1, sample, test, train, x_c, x_coarse, x_fine, x.new2_fine, x.new1_fine, x.new2_c, x.new2_coarse, x.new1_c, x.new1_coarse, a, b0, b1, b2, grain, n.repli, n.samp, n.train, nu, range, range_x.new2, range_x.new1, resolution, sple.size, variance, logit,
           file = here::here("results/simu/simulated_data.RData"))
