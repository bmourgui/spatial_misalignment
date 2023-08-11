#######################################%
#
# Script computing performance metrics
# regarding explanatory and predictive abilities
# of BEM, GLM, spGLM fitted on misaligned simulated
# data.
#
# Depends on 01_simulate_data.R 
# and fitting scripts (e.g. 02_run_glm.R)
# 
######################################%


library(magrittr)

load(here::here("results", "simu", "simulated_data.RData"))
load(file = here::here("results", "simu", "output_bem.RData"))
load(file = here::here("results", "simu", "output_glm.RData"))
load(file = here::here("results", "simu", "output_spglm+pred.RData"))

##### NEW ######
R <- min(length(l.bem.c1), 
         length(l.glm.c1), 
         length(l.spglm.c1)) # number of replications
# 1. Explanatory perf ----
## 1.1. Bias in ecological descriptors of SRCs ----

### 1.1.1. Extract model coef estimates ----
extract.beta.nimb <- function(out){
  summ <- MCMCvis::MCMCsummary(out)
  b0 <- c("param" = "b0", summ["b0", c(1, 4, 3, 5)])
  b1 <- c("param" = "b1", summ["b1", c(1, 4, 3, 5)])
  b2 <- c("param" = "b2", summ["b2", c(1, 4, 3, 5)])
  d <- as.data.frame(rbind(b0, b1, b2))
  colnames(d) <- c("param", "mean", "med", "lwr", "upr")
  d
}

compile.dbeta.nimb <- function(list.out = l.bem.c1, 
                               grain = "c1", 
                               true = c("b0" = b0, "b1" = b1, "b2" = b2)) {
  do.call(rbind,
          lapply(list.out,
                 function(out)extract.beta.nimb(out))
          ) |>
    dplyr::mutate("grain" = grain,
                  "repli" = rep(1:R, each = length(true)),
                  "true" = rep(true, R))
}

extract.beta.inla <- function(out){
  b0 <- c("param" = "b0", out$summary.fixed[1, c(1, 4, 3, 5)])
  b1 <- c("param" = "b1", out$summary.fixed[2, c(1, 4, 3, 5)])
  b2 <- c("param" = "b2", out$summary.fixed[3, c(1, 4, 3, 5)])
  d <- as.data.frame(rbind(b0, b1, b2))
  colnames(d) <- c("param", "mean", "med", "lwr", "upr")
  d
}

compile.dbeta.inla <- function(list.out = l.spglm.c1, 
                               grain = "c1", 
                               true = c("b0" = b0, "b1" = b1, "b2" = b2)) {
  do.call(rbind,
          lapply(list.out,
                 function(out)extract.beta.inla(out))
  ) |>
    dplyr::mutate("grain" = grain,
                  "repli" = rep(1:R, each = length(true)),
                  "true" = rep(true, R))
}


model <- c("glm", "bem", "spglm")
AG <- c("c1", "c2", "c3")
true.coef <- c("b0" = b0, "b1" = b1, "b2" = b2)

for (i in model){
  if (i == "spglm"){
    for (j in AG){
      lname <- paste("l", i, j, sep = ".")
      s <- compile.dbeta.inla(get(lname)[1:R], j, true.coef)
      assign(paste0("s.", j), s)
    }
  }else{
    for (j in AG){
      lname <- paste("l", i, j, sep = ".")
      s <- compile.dbeta.nimb(get(lname)[1:R], j, true.coef)
      assign(paste0("s.", j), s)
    }
  }
  s.m <- rbind(s.c1, s.c2, s.c3) |> 
    dplyr::mutate("model" = i)
  outname <- paste0("summary.beta.", i)
  assign(outname, s.m)
}
summary.beta <- rbind(summary.beta.bem, 
                      summary.beta.glm, 
                      summary.beta.spglm)

### 1.1.2. Calculate ecological descriptors ----
# i.e. optimum, ecological breadth, maximum response
summary.beta |> 
  dplyr::select(model, param, mean, grain, repli, true) |> 
  dplyr::mutate_at(c("mean", "true"), as.numeric) |> 
  dplyr::mutate_at(c("param"), as.character) |> 
  tidyr::pivot_longer(cols = c(mean, true), names_to = "type", values_to = "value") |> 
  tidyr::pivot_wider(names_from = param, values_from = value) |> 
  dplyr::rowwise() |> 
  dplyr::mutate("opti" = - b1 / (b2 * 2),
                "max" = plogis(b0 - (b1 ^ 2) / (4 * b2)),
                "tol" = 1/sqrt(-2*b2), 
                # equivalent way of calculating ecological width to that presented in appendix 1
                "width" = 2*tol*sqrt(2*(logit(max) - logit(.05))) # 0.05 represents the response value at which we calculate the SRC width
  ) -> summary.src


## 1.2. Errors in fitted values ----
predict.nimb <- function(out.nimb, x.pred){
  beta <- do.call(rbind, out.nimb)[, -4]
  X <- cbind(rep(1, length(x.pred)), x.pred, x.pred^2)
  posterior <- apply(beta, 1, function(b)as.numeric(plogis(X %*% b)))
  m <- apply(posterior, 1, mean)
  l <- apply(posterior, 1, function(x)quantile(x, 0.025))
  u <- apply(posterior, 1, function(x)quantile(x, 0.975))
  list("mean" = m, "lwr" = l, "upr" = u)
}

predict.inla <- function(out.inla, # model output
                         df,
                         xname) { 
  
  mesh <- out.inla[[length(out.inla)]]
  loc.pred <- df[, c("lon", "lat")]
  x.pred <- df[, xname]
  
  # create matrix A for prediction
  A.pred <- INLA::inla.spde.make.A(mesh, as.matrix(loc.pred))
  # create the spatial structure
  m.nodes <- out.inla$summary.random$u["mean"]
  m.field <- as.numeric(A.pred %*% as.data.frame(m.nodes)[, 1])
  l.nodes <- out.inla$summary.random$u["0.025quant"]
  l.field <- as.numeric(A.pred %*% as.data.frame(l.nodes)[, 1])
  u.nodes <- out.inla$summary.random$u["0.975quant"]
  u.field <- as.numeric(A.pred %*% as.data.frame(u.nodes)[, 1])
  # compute predictions as a function of environment and spatial field
  summ <- out.inla$summary.fixed
  beta <- (summ[1:3, c(1, 3, 5)])
  X <- cbind(rep(1, length(x.pred)), x.pred, x.pred^2)
  m <- as.numeric(plogis(X %*% beta[, 1] + m.field))
  l <- as.numeric(plogis(X %*% beta[, 2] + l.field))
  u <- as.numeric(plogis(X %*% beta[, 3] + u.field))
  list("mean" = m, "lwr" = l, "upr" = u)
}


ISfunction <- function(alpha = 0.05, lwr, upr, true){
  (upr - lwr) + 
    (2/alpha)*(lwr - true) * as.numeric(true < lwr)  + 
    (2/alpha)*(true - upr) * as.numeric(true > upr)
}

TSSmax <- function(y.true, psi.pred, threshold = seq(0.1, 0.95, 0.01)){
  threshold <- seq(min(psi.pred) + .01, max(psi.pred) - .01, .01)
  mat.psi <- matrix(rep(psi.pred, length(threshold)), 
                    ncol = length(threshold))
  mat.thresh <- matrix(rep(threshold, length(psi.pred)),
                       ncol = length(threshold),
                       byrow = T)
  mat.ypred <- apply(mat.psi > mat.thresh, 2, as.numeric)
  tss <- apply(mat.ypred, 2, function(x)spm::pred.acc(as.factor(y.true), factor(x, levels = c(0, 1)))[["tss"]])
  res <- data.frame("threshold" = threshold, "tss" = tss)
  unique(res[res$tss == max(tss), "tss"])
}

perf.metrics <- function(mean.pred, 
                        lwr.pred,
                        upr.pred,
                        psi.true,
                        y.true){
  sq.e <- (psi.true - mean.pred) ^ 2 # squared error
  qu <- quantile(psi.true, probs = c(.25, .75))
  rmse <- sqrt(mean(sq.e))
  rrmse <- sqrt(mean(sq.e)) / diff(qu) # 'relative' root mean squared error
  corr <- cor(psi.true, mean.pred)
  pis <- ISfunction(alpha = .05, lwr = lwr.pred, upr = upr.pred, true = psi.true)
  rpis <- ISfunction(alpha = .05, lwr = lwr.pred, upr = upr.pred, true = psi.true) / diff(qu)
  iw <- (upr.pred - lwr.pred) / diff(qu)
  roc <- pROC::roc(response = y.true, predictor = mean.pred, quiet = TRUE)
  maxtss <- TSSmax(y.true = y.true, psi.pred = mean.pred)
  brier <- DescTools::BrierScore(resp = y.true, pred = mean.pred)
  cbind("rmse" = rmse, 
        "rrmse" = rrmse,
        "corr" = corr, 
        "pis" = mean(pis), 
        "rpis" = mean(rpis), 
        "iw" = mean(iw),
        "brier" = brier,
        "auc" = as.numeric(pROC::auc(y.true, mean.pred)),
        "maxTSS" = maxtss)
}

model <- c("glm", "bem", "spglm")
AG <- c("c1", "c2", "c3")
PG <- c("coarse", "fine")
preddata <- c("train", "test", "pred.new1", "pred.new2")

n.mod <- length(model)
n.AG <- length(AG)
n.PG <- length(PG)
n.data <- length(preddata)
n.repli <- R

# store index values of inla predictions ##
id <- seq(0, 2400, 300)
id.fitAG <- (id[1]+1):id[2]
id.fitEG <- (id[2]+1):id[3]
id.interAG <- (id[3]+1):id[4]
id.interEG <- (id[4]+1):id[5]
id.extra1AG <- (id[5]+1):id[6]
id.extra1EG <- (id[6]+1):id[7]
id.extra2AG <- (id[7]+1):id[8]
id.extra2EG <- (id[8]+1):id[9]
##

summary.psi <- data.frame("model" = rep(model, each = n.AG*n.PG*n.data*n.repli),
                           "fit.grain" = rep(rep(AG, each = n.PG*n.data*n.repli), n.mod),
                           "pred.grain" = rep(rep(PG, each = n.data*n.repli), n.mod*n.AG),
                           "data.pred" = rep(rep(preddata, each = n.repli), n.mod*n.AG*n.PG),
                           "repli" = 1:n.repli,
                           "rmse" = NA,
                          "rrmse" = NA,
                           "corr" = NA,
                          "pis" = NA,
                          "rpis" = NA,
                          "iw" = NA,
                          "brier" = NA,
                          "auc" = NA,
                          "maxtss" = NA)

for (i in model){
  for (j in AG){
    for(k in preddata){
      for (l in PG){
        xname <- paste0("x_", ifelse(l == "coarse", j, l))
        PGname <- ifelse(l == "coarse", j, l)
        lname <- paste("l", i, j, sep = ".")
        fullname <- paste("rmse", k, l, i, j, sep = ".")
        
        if (i == "spglm"){
          metrics <- mapply(
            function(x, y){
              eqPG <- data.frame("d.name" = c("coarse", "fine"),
                                 "inla.name" = c("AG", "EG"))
              eqPreddata <- data.frame("d.name" = c("train", "test", "pred.new1", "pred.new2"),
                                       "inla.name" = c("fit", "inter", "extra1", "extra2"))
              id.name <- paste0("id.", eqPreddata[eqPreddata$d.name == k, "inla.name"], eqPG[eqPG$d.name == l, "inla.name"])
              p <- x$summary.fitted.values[get(id.name), ]
              perf.metrics(mean.pred = (p$mean),
                           lwr.pred = (p$`0.025quant`),
                           upr.pred = (p$`0.975quant`),
                           psi.true = y$psi_fine,
                           y.true = y$y_fine)
            },
            get(lname), get(k)[1:n.repli]
          )
          RMSE <- metrics[1, ]
          RRMSE <- metrics[2, ]
          CORR <- metrics[3, ]
          PIS <- metrics[4, ]
          RPIS <- metrics[5, ]
          IW <- metrics[6, ]
          BRIER <- metrics[7, ]
          AUC <- metrics[8, ]
          MAXTSS <- metrics[9, ]
        }else{
          metrics <- mapply(
            function(x, y){
              p <- predict.nimb(x, y[, xname])
              perf.metrics(mean.pred = p$mean,
                           lwr.pred = p$lwr,
                           upr.pred = p$upr,
                           psi.true = y$psi_fine,
                           y.true = y$y_fine)
            },
            get(lname), get(k)[1:n.repli]
          )
          RMSE <- metrics[1, ]
          RRMSE <- metrics[2, ]
          CORR <- metrics[3, ]
          PIS <- metrics[4, ]
          RPIS <- metrics[5, ]
          IW <- metrics[6, ]
          BRIER <- metrics[7, ]
          AUC <- metrics[8, ]
          MAXTSS <- metrics[9, ]
        }
        
        
        summary.psi <- summary.psi |> 
          dplyr::mutate("rmse" = replace(rmse, model == i & fit.grain == j & pred.grain == l & data.pred == k, RMSE),
                        "rrmse" = replace(rrmse, model == i & fit.grain == j & pred.grain == l & data.pred == k, RRMSE),
                        "corr" = replace(corr, model == i & fit.grain == j & pred.grain == l & data.pred == k, CORR),
                        "pis" = replace(pis, model == i & fit.grain == j & pred.grain == l & data.pred == k, PIS),
                        "rpis" = replace(rpis, model == i & fit.grain == j & pred.grain == l & data.pred == k, RPIS),
                        "iw" = replace(iw, model == i & fit.grain == j & pred.grain == l & data.pred == k, IW),
                        "brier" = replace(brier, model == i & fit.grain == j & pred.grain == l & data.pred == k, BRIER),
                        "auc" = replace(auc, model == i & fit.grain == j & pred.grain == l & data.pred == k, AUC),
                        "maxtss" = replace(maxtss, model == i & fit.grain == j & pred.grain == l & data.pred == k, MAXTSS))
      }
    }
  }
}


save(summary.beta,
     summary.src,
     summary.psi, 
     file = here::here("results", "simu", "summary_perf.RData"))
