#################################################################%
##
## Code to create outputs (e.g. figures) presented in manuscript
## and supporting information.
##
## Depends on 04_analyze_perfmetrics.R and its dependencies
## 
#################################################################%

library(ggplot2)

source(here::here("analysis", "04_analyze_perfmetrics.R"))

col.mod <- c("#DA8459", "#AB1488", "cadetblue")

# Main ----
## Figure 1: Illustration of area-to-point misalignment effect ----
# Produce the components for figure 1 which are then assembled on ppt

# Simulate true SER
xp <- seq(-3, 3, 0.1)
p <- plogis(3 + 1.5*xp + -3*(xp^2))
ser <- data.frame("grain" = rep(c("x1", "x2"), each = length(xp)),
                  "x" = xp,
                  "psi" = p)

df <- train[[2]]

# Plot raster of X at fine grain
data |> 
  dplyr::select(lon, lat, x_fine) |> 
  ggplot(aes(x = lon, y = lat, fill = x_fine)) +
  geom_raster() +
  theme_void() +
  annotate("text",
           x = 500, y = 500, 
           label = expression(X^{EG}),
           size = 6) +
  theme(legend.position = "none") +
  scale_fill_gradientn(colours = rev(grDevices::terrain.colors(50))) -> rast.x_eg

# Plot raster of probability of presence at EG
data |> 
  dplyr::select(lon, lat, psi_fine) |> 
  ggplot(aes(x = lon, y = lat, fill = psi_fine)) +
  geom_raster() +
  theme_void() +
  annotate("text",
           x = 500, y = 500, 
           label = expression(psi^{EG}),
           size = 6,
           color = "grey10") +
  theme(legend.position = "none") +
  scale_fill_gradientn(colours = rev(grDevices::terrain.colors(50))) -> rast.psi_eg

# Plot raster of X at coarser grain than EG
data |> 
  dplyr::select(lon, lat, paste0("x_c", 3)) |> 
  ggplot(aes(x = lon, y = lat, fill = get(paste0("x_c", 3)))) +
  geom_raster() +
  annotate("text",
           x = 500, y = 500, 
           label = bquote("X"^{"CG"}),
           size = 6,
           color = "black") +
  theme_void() +
  scale_fill_gradientn(colours = rev(grDevices::terrain.colors(50))) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5, vjust = - 2)) -> rast.x_cg

# Plot simulated probabilities at EG in function of X at fine grain
df |>
  ggplot() +
  geom_line(data = ser[ser$grain == "x1",], aes(x = x, y = psi)) +
  geom_point(aes(x = x_c1, y = psi_fine, color = as.factor(y_fine)), alpha = .4) +
  xlim(range(xp)) +
  scale_color_manual(values = c("grey20", "forestgreen"), name = bquote(Y^{EG})) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylab(bquote(psi^{EG})) +
  xlab(bquote(X^{CG})) -> ser.cg1

# Plot simulated probabilities at EG in function of X at coarse grain
df |>
  ggplot() +
  geom_line(data = ser[ser$grain == "x2",], aes(x = x, y = psi)) +
  geom_point(aes(x = x_c3, y = psi_fine, color = as.factor(y_fine)), alpha = .4) +
  xlim(range(xp)) +
  scale_color_manual(values = c("grey20", "forestgreen"), name = bquote(Y^{EG})) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none") +
  ylab(bquote(psi^{EG})) +
  xlab(bquote(X^{CG})) -> ser.cg2

leg <- cowplot::get_legend(
  ser.cg2 +
    theme(legend.position = "right")
)

# Combine the different plots in one final plot
comb.rast.cg1 <- cowplot::plot_grid(NULL, rast.x_eg, rast.psi_eg, 
                                    ncol = 3)

full.ser.cg1 <- cowplot::plot_grid(comb.rast.cg1, ser.cg1, 
                                   nrow = 2, 
                                   rel_heights = c(0.2, 1))
comb.rast.cg2 <- cowplot::plot_grid(NULL, rast.x_cg, rast.psi_eg, 
                                    ncol = 3)

full.ser.cg2 <- cowplot::plot_grid(comb.rast.cg2, ser.cg2, 
                                   nrow = 2, 
                                   rel_heights = c(0.2, 1))
comb.sers <- cowplot::plot_grid(full.ser.cg1, full.ser.cg2, leg, ncol = 3, rel_widths = c(1, 1, .45))

ggsave(here::here("results", "img", "figure1.pdf"),
       plot = comb.sers,
       width = 1476,
       height = 1476,
       units = "px",
       dpi = 300)

## Figure 2: Simulation and analysis framework ----
# Create the different components of the figure and then assembled them on ppt

### Fig. 2.1.: Spatial autocorrelation curve with rho = 50 ----
dis <- seq(0.01, range*3, length.out = 200)
corr <- 2 ^ (1 - nu) / gamma(nu) * (kappa * dis) ^ nu * besselK(x = dis * kappa, nu = nu)

data.frame("dis" = dis, "corr" = corr) |> 
  ggplot(aes(x = dis, y = corr)) +
  geom_line(linewidth = 1.2, colour = "blue") +
  geom_segment(aes(x = 0, xend = range, 
                   y = 0.13, yend = 0.13),
               linetype = 2) +
  geom_segment(aes(x = range, xend = range, 
                   y = 0, yend = 0.13),
               linetype = 2) +
  #annotate("text", x = range/2, y = 0.15,
  #         label = paste0("Spatial range = ", range, " units"),
  #         size = 4) +
  scale_x_continuous(expand = c(0, 0.2),
                     limits = c(0, range*3+1)) +
  scale_y_continuous(expand = c(0.005, 0),
                     limits = c(0, 1.05)) +
  ylab("Spatial autocorrelation") +
  xlab("Distance between two points") +
  theme_classic(base_size = 1.54) +
  theme(text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), type = "closed")))

ggsave(here::here("results", "img", "fig2-1_autocorr_curve.png"), width = 50, height = 50, units = "mm")

### Fig. 2.2.: Raster X^EG ----
data |> 
  dplyr::select(lon, lat, x_fine) |> 
  ggplot(aes(x = lon, y = lat, fill = x_fine)) +
  geom_raster() +
  theme_void() +
  #annotate("text",
  #         x = 500, y = 500, 
  #         label = expression(X^{EG}),
  #         size = 1.54) +
  theme(legend.position = "none") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")

ggsave(here::here("results", "img", "fig2-2_Xeg_field.png"), width = 100, height = 100, units = "mm")

### Fig 2.3.: Raster Psi^EG ----
data |> 
  dplyr::select(lon, lat, psi_fine) |> 
  ggplot(aes(x = lon, y = lat, fill = psi_fine)) +
  geom_raster() +
  theme_void() +
  #annotate("text",
  #         x = 500, y = 500, 
  #         label = expression(psi^{EG}),
  #         size = 1.54,
  #         color = "grey10") +
  theme(legend.position = "none") +
  scale_fill_gradient(low = "white", high = "#00A6BB")
ggsave(here::here("results", "img", "fig2-3_PSIeg_field.png"), width = 80, height = 80, units = "mm")

### Fig 2.4.: Raster Y^EG ----
data |> 
  dplyr::select(lon, lat, y_fine) |> 
  dplyr::mutate("y_fine" = as.factor(y_fine)) |> 
  ggplot(aes(x = lon, y = lat, fill = y_fine)) +
  geom_tile() +
  theme_void() +
  #annotate("text",
  #         x = 500, y = 500, 
  #         label = expression(Y^{EG}),
  #         size = 1.54,
  #         color = "grey10") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("white", "forestgreen"))
ggsave(here::here("results", "img", "fig2-4_Yeg_field.png"), width = 100, height = 100, units = "mm")


### Fig. 2.5-7.: Three Rasters of X^CG ----
CG <- paste0("CG = spatial range / ", c(10, 2, 1))
for (i in 1:3){
  set.seed(423)
  dt <- data |> 
    dplyr::slice_sample(n = 500)
  
  dt_coarse <- terra::extract(x_coarse[[i]],
                              dt[, c("lon", "lat")],
                              ID = FALSE,
                              xy = TRUE) |> 
    as.data.frame()
  colnames(dt_coarse) <- c(paste0("x_c", i), "lon", "lat")
  
  data |> 
    dplyr::select(lon, lat, paste0("x_c", i)) |> 
    ggplot(aes(x = lon, y = lat, fill = get(paste0("x_c", i)))) +
    #ggtitle(label = AG[i]) +
    geom_raster(alpha = .8) +
    #geom_tile(data = dt_coarse, colour="grey20", linewidth = 0.5) +
    #geom_point(data = dt, colour="grey40", size = .5) +
    #annotate("text",
    #         x = 500, y = 500, 
    #         label = bquote("X"^{"AG"[.(i)]}),
    #         size = 1.54,
    #         color = "black") +
    theme_void() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5, vjust = - 2)) -> p
  print(p)
  ggsave(here::here("results", "img", paste0("fig2-", (i+4), "_Xag_field.png")), width = 100, height = 100, units = "mm")
}

### Fig. 2.8.: Raster X^CG with sampling cells ----
dt <- data |> 
  dplyr::slice_sample(n = 50)

dt_coarse <- terra::extract(x_coarse[[3]],
                            dt[, c("lon", "lat")],
                            ID = FALSE,
                            xy = TRUE) |> 
  as.data.frame()
colnames(dt_coarse) <- c(paste0("x_c", 3), "lon", "lat")

data |> 
  dplyr::select(lon, lat, paste0("x_c", 3)) |> 
  ggplot(aes(x = lon, y = lat, fill = get(paste0("x_c", 3)))) +
  geom_raster(alpha = 1) +
  geom_tile(data = dt_coarse, colour="grey20", linewidth = 0.5) +
  theme_void() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5, vjust = - 2))

ggsave(here::here("results", "img", "fig2-8_Xag3_sampled.png"), width = 100, height = 100, units = "mm")

### Fig. 2.9.: Raster Y^EG with sampling points ----
data |> 
  dplyr::select(lon, lat, y_fine) |> 
  dplyr::mutate("y_fine" = as.factor(y_fine)) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(alpha = 1, aes(fill = y_fine)) +
  geom_point(data = dt, colour = "darkred", size = 4) +
  theme_void() +
  scale_fill_manual(values = c("white", "forestgreen")) +
  theme(legend.position = "none")

ggsave(here::here("results", "img", "fig2-9_Yeg_sampled.png"), width = 100, height = 100, units = "mm")


### Fig. 2.10-12: 3 "Estimated" SER ----
x <- seq(min(data$x_fine), max(data$x_fine), length.out = 200)

for (i in 1:3){
  opti.e <- (opti + rnorm(1, 0, 0.3))
  a.e <- a + rnorm(1, 0, 0.5)
  psi.max.e <- plogis(a.e)
  psi2 <- plogis(a.e - (x - opti.e)^2 / (2 * (tau + rnorm(1, 0, 0.2)) ^ 2))
  y2 <- rbinom(length(psi2), 1, psi2)
  data.frame("X" = x, "psi" = psi2) |> 
    ggplot(aes(x = X)) +
    geom_line(aes(y = psi), linewidth = 1.2, color = col.mod[i]) +
    #geom_jitter(data = dt, aes(x = x_fine, y = y_fine), 
    #            color = alpha("darkred", .6), 
    #            size = 1.5.4,
    #            width = 0,
    #            height = .01) +
    theme_classic() +
    theme(legend.position = "none",
          axis.ticks = element_blank(),
          text = element_blank(),
          axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), type = "closed"))) +
    scale_x_continuous(expand = c(0, 0.2)) +
    scale_y_continuous(expand = c(0.005, 0),
                       limits = c(0, 1.05)) -> p
  print(p)
  ggsave(here::here("results", "img", paste0("fig2-", (i+9), "_SER_estim.png")), width = 60, height = 60, units = "mm")
}


### Fig. 2.13: Explanatory Power ----
# Comparison between true and estimates SERs
opti.e <- (opti-.5)
a.e <- a-1
psi.max.e <- plogis(a.e)
psi <- plogis(a - (x - opti)^2 / (2 * (tau) ^ 2))
psi2 <- plogis(a.e - (x - opti.e)^2 / (2 * (tau + 0.3) ^ 2))
y2 <- rbinom(length(psi2), 1, psi2)
data.frame("X" = c(x, x), "psi" = c(psi, psi2), "type" = rep(c("t", "e"), each = length(x))) |> 
  #dplyr::filter(type == "e") |> 
  ggplot(aes(x = X)) +
  geom_line(aes(color = type, y = psi), linewidth = 1.2) +
  #geom_point(aes(y = as.numeric(rep(y2, 2))), color = "darkred") +
  scale_color_manual(values = c(col.mod[2], "forestgreen")) +
  #geom_segment(aes(x = opti.e, xend = opti.e, y = psi.max.e, yend = 0),
  #             linewidth = .2,
  #             linetype = "longdash",
  #             colour = "grey50",
  #             arrow = grid::arrow(length = unit(0.25, "cm"), type = "open")) +
  #geom_segment(aes(x = opti.e, xend = -5, y = psi.max.e, yend = psi.max.e),
  #             linewidth = .2,
  #             linetype = "longdash",
  #             colour = "grey50",
  #             arrow = grid::arrow(length = unit(0.25, "cm"), type = "open")) +
  #geom_segment(aes(x = (opti.e-(width/2)), xend = (opti.e+(width/2)), y = 0.05, yend = 0.05),
#             linewidth = .2,
#             linetype = "longdash",
#             colour = "grey50",
#             arrow = grid::arrow(length = unit(0.25, "cm"), ends = "both", type = "open")) +
#annotate("text", x = opti.e+.5, y = .03, label = expression(hat(theta)), size = 7) +
#annotate("text", x = opti.e, y = .08, label = expression(hat(omega)), size = 7) +
#annotate("text", x = -4, y = (psi.max.e + .05), label = expression(hat(psi)["max"]), size = 7) +
theme_classic() +
  theme(legend.position = "none",
        #axis.title = element_text(size = 1.58),
        axis.ticks = element_blank(),
        text = element_blank(),
        axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), type = "closed"))) +
  scale_x_continuous(expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0.005, 0),
                     limits = c(0, 1.05)) 

ggsave(here::here("results", "img", "fig2-13_SRC_comparison.png"), width = 100, height = 100, units = "mm")


### Fig. 2.14-15: curves autocorrelation  of the two new environments ----
R <- c("smooth" = range_x.new1, "hete" = range_x.new2)
dis <- seq(0.01, R[1]*3, length.out = 400)
for (i in 1:2){
  kappa <- sqrt(8*nu) / R[i]
  corr <- 2 ^ (1 - nu) / gamma(nu) * (kappa * dis) ^ nu * besselK(x = dis * kappa, nu = nu)
  
  data.frame("dis" = dis, "corr" = corr) |> 
    ggplot(aes(x = dis, y = corr)) +
    geom_line(linewidth = 1.2, colour = "blue") +
    geom_segment(aes(x = 0, xend = R[i], 
                     y = 0.13, yend = 0.13),
                 linetype = 2) +
    geom_segment(aes(x = R[i], xend = R[i], 
                     y = 0, yend = 0.13),
                 linetype = 2) +
    scale_x_continuous(expand = c(0, 0.2),
                       limits = c(0, max(R[i]*2, range*3)+1)) +
    scale_y_continuous(expand = c(0.005, 0),
                       limits = c(0, 1.05)) +
    theme_classic() +
    theme(text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), type = "closed"))) -> p
  print(p)
  ggsave(here::here("results", "img", paste0("fig2-", (i+13), "autocorr_curve", names(R)[i], ".png")), width = 50, height = 50, units = "mm")
}

### Fig. 2.16-18: Rasters X^EG of the three environments ----
data <- data |> 
  dplyr::mutate("heterogeneity" = "Intermediate")
d.new1 <- d.new1 |> 
  dplyr::mutate("heterogeneity" = "Low")
d.new2 <- d.new2 |> 
  dplyr::mutate("heterogeneity" = "High")

df <- rbind(data, d.new1, d.new2) |> 
  dplyr::select(-c(psi_fine:x_c2)) |> 
  tidyr::pivot_longer(x_fine:x_c3, names_to = "grain", values_to = "x")

for (h in unique(df$heterogeneity)){
  df |> 
    dplyr::filter(grain == "x_fine", heterogeneity == h) |> 
    ggplot(aes(x = lon, y = lat)) +
    geom_raster(aes(fill = x)) +
    #geom_point(data = dt, colour = "darkred", size = 1.5) +
    theme_void() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) -> p
  print(p)
  ggsave(here::here("results", "img", paste0("Xpred_", h, ".png")), width = 100, height = 100, units = "mm")
}

### Fig. 2.19.: Raster X_new2^CG with sampling cells ----
set.seed(23)
sp <- sample(1:nrow(data), 50)
dt <- df |> 
  dplyr::group_by(heterogeneity, grain) |> 
  dplyr::slice(sp)

xcoarse <- x.new2_coarse
dt_coarse <- terra::extract(xcoarse[[3]],
                            dt[, c("lon", "lat")],
                            ID = FALSE,
                            xy = TRUE) |> 
  as.data.frame()
colnames(dt_coarse) <- c("x", "lon", "lat")

df |> 
  dplyr::filter(grain == "x_c3" & heterogeneity == "High") |> 
  ggplot(aes(x = lon, y = lat, fill = x)) +
  geom_raster() +
  geom_tile(data = dt_coarse, colour="grey20", linewidth = 0.5) +
  theme_void() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave(here::here("results", "img", "fig2-19_Xhete_AG.png"),  width = 100, height = 100, units = "mm")

### Fig. 2.20.: Schema performance metrics ----
# Plot Fake Predicted vs Observed probabilities
set.seed(234)
x <- rnorm(300)
y <- x + rnorm(300, 0, 0.5)
data.frame(x, y) |> 
  ggplot(aes(x = x, y = y)) +
  geom_point(colour = "darkred", size = 1.5) +
  ylab(expression(hat(psi)^{EG})) +
  xlab(expression(psi^{EG})) +
  theme_classic() +
  theme(legend.position = "none",
        #axis.title = element_text(size = 1.58),
        text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), type = "closed"))) 
ggsave(here::here("results", "img", "fig2-20_EGpred_vs_true.png"), width = 100, height = 100, units = "mm")


## Estimation of species-environment relationships ----
### Figure 3: Relative bias in SER estimates ----
d.bias.ser |> 
  ggplot() +
  geom_point(aes(x = grain, y = relative.bias, color = model),
             position = position_jitterdodge(dodge.width = 0.3,
                                             jitter.width = .1),
             alpha = 0.2) +
  facet_wrap(.~param, scales = "free",
             labeller = labeller(param = c("max" = "a) Maximum",
                                           "opti" = "b) Optimum",
                                           "width" = "c) Ecological width"))) +
  ylab("Relative bias") +
  xlab("Grain size") +
  scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
  geom_hline(aes(yintercept = 0)) +
  theme_bw() +
  geom_point(data = d.summary.relbias.ser, 
             aes(x = grain, y = avg, color = model),
             size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(data = d.summary.relbias.ser, 
                aes(x = grain, ymin = avg - sd, ymax = avg + sd, color = model),
                width = 0.2,
                position = position_dodge(width = 0.3)) +
  theme(legend.position = "top",
        text = element_text(size = 8)) +
  scale_color_manual(values = col.mod)

ggsave(here::here("results/img/figure3.pdf"),
       width = 1476,
       height = 1204,
       units = "px",
       dpi = 300)

# Create a table with summary results for relative bias
# (used for presenting results in the text)
d.summary.relbias.ser |> 
  dplyr::mutate("avg.sd" = paste0(round(avg, 2), " (", round(sd, 2), ")")) |> 
  dplyr::select(-c(avg, sd)) |> 
  tidyr::pivot_wider(values_from = avg.sd, names_from = model) |> 
  write.csv(file = here::here("results", "table", "relative_bias_eco_param.csv"))

### Figure 4: AUC and Brier for train data ----
for (m in c("brier", "auc")) {
    
  dMean <- d.mean.perf.expla |> 
    dplyr::filter(pred.grain == "coarse", # predictions from covariate grain environmental values
                  metric == m)
  
    summary.psi |>
      dplyr::filter() |> 
      tidyr::pivot_longer(rmse:maxtss, names_to = "metric", values_to = "value") |> 
      dplyr::filter(data.pred == "train" & 
                      pred.grain == "coarse" & # predictions from covariate grain environmental values
                      metric == m) |> 
      ggplot() +
      geom_point(aes(x = fit.grain, y = value, color = model),
                 position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.1),
                 alpha = 0.2)  +
      scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
      ylab(m) +
      xlab("Grain size") +
      theme_bw() +
      geom_point(data = dMean,
                 aes(x = fit.grain, y = avg, color = model),
                 size = 2,
                 position = position_dodge(width = 0.4)) +
      geom_errorbar(data = dMean,
                    aes(x = fit.grain, ymin = avg - sd, ymax = avg + sd, color = model),
                    width = 0.2,
                    position = position_dodge(width = 0.4)) +
      theme(legend.position = "top",
            text = element_text(size = 8)) +
      scale_color_manual(values = col.mod) -> gp
    assign(paste0("fig4.", m), gp)
  }

fig4.auc +
  ylab("AUC") +
  xlab("") -> fig4.A

fig4.brier +
  ylab("Brier score") +
  theme(legend.position = "none") +
  ggtitle("") -> fig4.B

cowplot::plot_grid(fig4.A, fig4.B,
                   labels = c("(a)", "(b)"),
                   ncol = 1)
ggsave(here::here("results", "img", "figure4.pdf"),
       width = 1476,
       height = 1476,
       units = "px",
       dpi = 300)

### Figure 5: AUC, Brier and PIS for test data ----
for (m in c("pis", "brier", "auc")) {
    dMean.metric.pred_c <- d.mean.perf.pred |>
      dplyr::filter(data.pred != "train" & 
                      pred.grain == "coarse" & # predictions from CG
                      metric == m)
    
    summary.psi |>
      dplyr::filter() |> 
      tidyr::pivot_longer(rmse:maxtss, names_to = "metric", values_to = "value") |> 
      dplyr::filter(data.pred != "train" & 
                      pred.grain == "coarse" & 
                      metric == m) |> 
      ggplot() +
      geom_point(aes(x = fit.grain, y = value, color = model),
                 position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1),
                 alpha = 0.2) +
      facet_grid(. ~ factor(data.pred, 
                            levels = c("pred.new1", "test", "pred.new2"),
                            labels = c("pred.new1" = "Low ~~ (X[new1])",
                                       "test" = "Intermediate ~~ (X)",
                                       "pred.new2" = "High ~~ (X[new2])")),
                 scales = "free",
                 labeller = labeller(.default = label_parsed))  +
      scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
      ylab(m) +
      xlab("Grain size") +
      theme_bw() +
      geom_point(data = dMean.metric.pred_c,
                 aes(x = fit.grain, y = avg, color = model),
                 size = 2,
                 position = position_dodge(width = 0.65)) +
      geom_errorbar(data = dMean.metric.pred_c,
                    aes(x = fit.grain, ymin = avg - sd, ymax = avg + sd, color = model),
                    width = 0.2,
                    position = position_dodge(width = 0.65)) +
      theme(legend.position = "top",
            text = element_text(size = 8)) +
      scale_color_manual(values = col.mod) -> gp
    assign(paste0("fig5.", m), gp)
  }

fig5.auc +
  ggtitle("Environmental spatial heterogeneity") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.65),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank()) +
  ylab("AUC") -> fig5.A

fig5.brier +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("Brier score") -> fig5.B

fig5.pis +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("PIS") -> fig5.C

# combine the three metrics
fig5 <- cowplot::plot_grid(fig5.A, fig5.B, fig5.C,
                           ncol = 1,
                           rel_heights = c(1, 0.75, 0.9))

# add legend
legend_b <- cowplot::get_legend(
  fig5.C + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
cowplot::plot_grid(fig5, legend_b, ncol = 1, rel_heights = c(1, .05))

# save
ggsave(here::here("results", "img", "figure5.pdf"),
       width = 1476,
       height = 1676,
       units = "px",
       dpi = 300)


# Supporting information ----

## Illustration of SER parameters (opt, max, width) ----
x <- seq(min(data$x_fine), max(data$x_fine), length.out = 200)
psi <- plogis(a - (x - opti)^2 / (2 * (tau) ^ 2))
data.frame("X" = x, "psi" = psi) |> 
  ggplot(aes(x = X, y = psi)) +
  geom_line(color = "forestgreen", linewidth = 1.2) +
  geom_segment(aes(x = opti, xend = opti, y = psi.max, yend = 0),
               linewidth = .2,
               linetype = "longdash",
               colour = "grey50",
               arrow = grid::arrow(length = unit(0.25, "cm"), type = "open")) +
  geom_segment(aes(x = opti, xend = -5, y = psi.max, yend = psi.max),
               linewidth = .2,
               linetype = "longdash",
               colour = "grey50",
               arrow = grid::arrow(length = unit(0.25, "cm"), type = "open")) +
  geom_segment(aes(x = (opti-(width/2)), xend = (opti+(width/2)), y = 0.05, yend = 0.05),
               linewidth = .2,
               linetype = "longdash",
               colour = "grey50",
               arrow = grid::arrow(length = unit(0.25, "cm"), ends = "both", type = "open")) +
  annotate("text", x = opti+.3, y = .02, label = expression(theta), size = 7) +
  annotate("text", x = opti, y = .08, label = expression(omega), size = 7) +
  annotate("text", x = -4, y = (psi.max + .04), label = expression(psi["max"]), size = 7) +
  theme_classic() +
  ylab(label = expression(psi)) +
  xlab(label = expression(X)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 22),
        axis.ticks = element_blank(),
        axis.line = element_line(arrow = grid::arrow(length = unit(0.1, "cm"), type = "closed"))) +
  scale_x_continuous(expand = c(0, 0.2)) +
  scale_y_continuous(expand = c(0.005, 0),
                     limits = c(0, 1.05)) 
ggsave(here::here("results", "img", "figS1.png"), width = 1476, height = 1600, units = "px")

## Extended results - RMSE and CrI width ----
### Explanatory performances ----
for (m in c("rmse", "iw")) {
  
  dMean <- d.mean.perf.expla |> 
    dplyr::filter(pred.grain == "coarse", # predictions from EG environmental values
                  metric == m)
  
  summary.psi |>
    dplyr::filter() |> 
    tidyr::pivot_longer(rmse:maxtss, names_to = "metric", values_to = "value") |> 
    dplyr::filter(data.pred == "train" & 
                    pred.grain == "coarse" & # predictions from EG environmental values
                    metric == m) |> 
    ggplot() +
    geom_point(aes(x = fit.grain, y = value, color = model),
               position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.1),
               alpha = 0.2)  +
    scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
    ylab(m) +
    xlab("Grain size") +
    theme_bw() +
    geom_point(data = dMean,
               aes(x = fit.grain, y = avg, color = model),
               size = 2,
               position = position_dodge(width = 0.4)) +
    geom_errorbar(data = dMean,
                  aes(x = fit.grain, ymin = avg - sd, ymax = avg + sd, color = model),
                  width = 0.2,
                  position = position_dodge(width = 0.4)) +
    theme(legend.position = "top",
          text = element_text(size = 8)) +
    scale_color_manual(values = col.mod) -> gp
  assign(paste0("figS4.", m), gp)
}

figS4.rmse +
  ylab("RMSE") +
  xlab("") -> figS4.A

figS4.iw +
  ylab("Credible interval width") +
  theme(legend.position = "none") +
  ggtitle("") -> figS4.B

cowplot::plot_grid(figS4.A, figS4.B,
                   labels = c("(a)", "(b)"),
                   ncol = 1)
ggsave(here::here("results", "img", "figureS2.pdf"),
       width = 1476,
       height = 1476,
       units = "px",
       dpi = 300)

### Predictive performances ----
for (m in c("rmse", "iw")) {
  dMean.metric.pred_c <- d.mean.perf.pred |>
    dplyr::filter(data.pred != "train" & 
                    pred.grain == "coarse" & # predictions from CG
                    metric == m)
  
  summary.psi |>
    dplyr::filter() |> 
    tidyr::pivot_longer(rmse:maxtss, names_to = "metric", values_to = "value") |> 
    dplyr::filter(data.pred != "train" & 
                    pred.grain == "coarse" & 
                    metric == m) |> 
    ggplot() +
    geom_point(aes(x = fit.grain, y = value, color = model),
               position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1),
               alpha = 0.2) +
    facet_grid(. ~ factor(data.pred, 
                          levels = c("pred.new1", "test", "pred.new2"),
                          labels = c("pred.new1" = "Low ~~ (X[new1])",
                                     "test" = "Intermediate ~~ (X)",
                                     "pred.new2" = "High ~~ (X[new2])")),
               scales = "free",
               labeller = labeller(.default = label_parsed))  +
    scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
    ylab(m) +
    xlab("Grain size") +
    theme_bw() +
    geom_point(data = dMean.metric.pred_c,
               aes(x = fit.grain, y = avg, color = model),
               size = 2,
               position = position_dodge(width = 0.65)) +
    geom_errorbar(data = dMean.metric.pred_c,
                  aes(x = fit.grain, ymin = avg - sd, ymax = avg + sd, color = model),
                  width = 0.2,
                  position = position_dodge(width = 0.65)) +
    theme(legend.position = "top",
          text = element_text(size = 8)) +
    scale_color_manual(values = col.mod) -> gp
  assign(paste0("figS5.", m), gp)
}

figS5.rmse +
  ggtitle("Environmental spatial heterogeneity") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.65),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank()) +
  ylab("RMSE") -> figS5.A

figS5.iw +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("Credible interval width") -> figS5.B

# combine the three metrics
figS5 <- cowplot::plot_grid(figS5.A, figS5.B,
                           ncol = 1,
                           rel_heights = c(1, 0.9))

# add legend
cowplot::plot_grid(figS5, legend_b, ncol = 1, rel_heights = c(1, .05))

# save
ggsave(here::here("results", "img", "figureS3.pdf"),
       width = 1476,
       height = 1676,
       units = "px",
       dpi = 300)

## Predictive performances when predicting from EG covariate ----
### Explanatory performances ----
for (m in c("brier", "auc")) {
  
  dMean <- d.mean.perf.expla |> 
    dplyr::filter(pred.grain == "fine", # predictions from EG environmental values
                  metric == m)
  
  summary.psi |>
    dplyr::filter() |> 
    tidyr::pivot_longer(rmse:maxtss, names_to = "metric", values_to = "value") |> 
    dplyr::filter(data.pred == "train" & 
                    pred.grain == "fine" & # predictions from EG environmental values
                    metric == m) |> 
    ggplot() +
    geom_point(aes(x = fit.grain, y = value, color = model),
               position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.1),
               alpha = 0.2)  +
    scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
    ylab(m) +
    xlab("Grain size") +
    theme_bw() +
    geom_point(data = dMean,
               aes(x = fit.grain, y = avg, color = model),
               size = 2,
               position = position_dodge(width = 0.4)) +
    geom_errorbar(data = dMean,
                  aes(x = fit.grain, ymin = avg - sd, ymax = avg + sd, color = model),
                  width = 0.2,
                  position = position_dodge(width = 0.4)) +
    theme(legend.position = "top",
          text = element_text(size = 8)) +
    scale_color_manual(values = col.mod) -> gp
  assign(paste0("figS4.", m), gp)
}

figS4.auc +
  ylab("AUC") +
  xlab("") -> figS4.A

figS4.brier +
  ylab("Brier score") +
  theme(legend.position = "none") +
  ggtitle("") -> figS4.B

cowplot::plot_grid(figS4.A, figS4.B,
                   labels = c("(a)", "(b)"),
                   ncol = 1)
ggsave(here::here("results", "img", "figureS4.pdf"),
       width = 1476,
       height = 1476,
       units = "px",
       dpi = 300)

### Predictive performances ----
for (m in c("pis", "brier", "auc")) {
  dMean.metric.pred_c <- d.mean.perf.pred |>
    dplyr::filter(data.pred != "train" & 
                    pred.grain == "fine" & # predictions from CG
                    metric == m)
  
  summary.psi |>
    dplyr::filter() |> 
    tidyr::pivot_longer(rmse:maxtss, names_to = "metric", values_to = "value") |> 
    dplyr::filter(data.pred != "train" & 
                    pred.grain == "fine" & 
                    metric == m) |> 
    ggplot() +
    geom_point(aes(x = fit.grain, y = value, color = model),
               position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1),
               alpha = 0.2) +
    facet_grid(. ~ factor(data.pred, 
                          levels = c("pred.new1", "test", "pred.new2"),
                          labels = c("pred.new1" = "Low ~~ (X[new1])",
                                     "test" = "Intermediate ~~ (X)",
                                     "pred.new2" = "High ~~ (X[new2])")),
               scales = "free",
               labeller = labeller(.default = label_parsed))  +
    scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
    ylab(m) +
    xlab("Grain size") +
    theme_bw() +
    geom_point(data = dMean.metric.pred_c,
               aes(x = fit.grain, y = avg, color = model),
               size = 2,
               position = position_dodge(width = 0.65)) +
    geom_errorbar(data = dMean.metric.pred_c,
                  aes(x = fit.grain, ymin = avg - sd, ymax = avg + sd, color = model),
                  width = 0.2,
                  position = position_dodge(width = 0.65)) +
    theme(legend.position = "top",
          text = element_text(size = 8)) +
    scale_color_manual(values = col.mod) -> gp
  assign(paste0("figS5.", m), gp)
}

figS5.auc +
  ggtitle("Environmental spatial heterogeneity") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.65),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank()) +
  ylab("AUC") -> figS5.A

figS5.brier +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("Brier score") -> figS5.B

figS5.pis +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank()) +
  ylab("PIS") -> figS5.C

# combine the three metrics
figS5 <- cowplot::plot_grid(figS5.A, figS5.B, figS5.C,
                           ncol = 1,
                           rel_heights = c(1, 0.75, 0.9))

# add legend
cowplot::plot_grid(figS5, legend_b, ncol = 1, rel_heights = c(1, .05))

# save
ggsave(here::here("results", "img", "figureS5.pdf"),
       width = 1476,
       height = 1676,
       units = "px",
       dpi = 300)


## Estimates of species response curves ----
d.src |> 
  ggplot(aes(x = x.pred, group = interaction(repli, model))) + 
  geom_line(aes(y = psi.pred, col = model)) +
  geom_line(aes(y = psi.true), col = "forestgreen") +
  facet_grid(fit.grain ~ model,
             labeller = labeller(fit.grain = c("c1" = "Fine",
                                               "c2" = "Medium",
                                               "c3" = "Coarse"),
                                 model = c("bem" = "BEM",
                                           "glm" = "GLM",
                                           "spglm" = "spGLM"))) +
  xlab(expression(X^EG)) +
  ylim(c(0, 1)) +
  ylab(expression(psi^EG)) +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size = 8)) +
  scale_color_manual(values = col.mod)

ggsave(here::here("results", "img", "figS6.pdf"),
       width = 1476,
       height = 1204,
       units = "px")

## Loss of fine-grain variability when coarsening the covariate ----
### Figure SX: Loss of intra CG cell variability with coarsening ----
var_intra.medium <- var_loss.medium[["var_intra.cell"]]
var_intra.smooth <- var_loss.smooth[["var_intra.cell"]]
var_intra.hete <- var_loss.hete[["var_intra.cell"]]

pdf(here::here("results", "img", "figS7.pdf"),
    width = 8)
plot(d, 1 - unlist(var_intra.medium), 
     ylim = c(0, 1), 
     xlab = "Grain size", 
     ylab = "1 - var(intracell)", 
     main = "Loss of fine-grain variability", 
     pch = 19, 
     cex = 0.4)
points(d, 1 - unlist(var_intra.smooth), 
       pch = 19, 
       cex = 0.4, 
       col = "forestgreen")
points(d, 1 - unlist(var_intra.hete), 
       pch = 19, 
       cex = 0.4, 
       col = "deeppink3")
abline(v = grain, lty = 2)
#abline(h = (1 - unlist(var_intra.medium))[grain-1], lty = 2)
mtext(text = c("Fine", "Medium", "Coarse"),
      at = grain, 
      side = 1)
legend(x = 150,
       y = 0.95,
       title = "Spatial heterogeneity",
       legend = c("Low", "Intermediate", "High"),
       col = c("forestgreen", "black", "deeppink3"),
       pch = 19)
dev.off()

# Not used ----
## Other version of Figure 1 with estimated SERs ----
# Plot simulated probabilities of presence at EG in function of covariate at CG with SERs estimated
d.src |> 
  dplyr::filter(repli == 2) |> 
  dplyr::mutate("grain" = paste0("x_", fit.grain)) |> 
  dplyr::select(-c(psi.true, fit.grain)) -> src.train1

train[[2]] |> 
  tidyr::pivot_longer(cols = c(x_c1:x_c3),
                      names_to = "grain",
                      values_to = "x") |> 
  ggplot() +
  geom_point(aes(x = x, y = psi_fine),
             alpha = 0.4,
             color = "darkgreen") +
  geom_line(data = src.train1, 
            aes(x = x.pred, y = psi.pred, color = model),
            size = 1) +
  scale_color_manual(values = col.mod) + 
  facet_grid(. ~ factor(grain, 
                        levels = c("x_fine", "x_c1", "x_c2", "x_c3"),
                        labels = c("CG = EG (not fitted)", "Fine CG", "Medium CG", "Coarse CG"))) +
  ylab(expression(psi^EG)) +
  xlab(expression(X^CG)) +
  theme_bw()

## Figure 3 with regression coefficients ----
d.bias.beta |> 
  ggplot() +
  geom_point(aes(x = grain, y = mean, color = model),
             position = position_jitterdodge(dodge.width = 0.3,
                                             jitter.width = .1),
             alpha = 0.2) +
  facet_wrap(.~param, scales = "free",
             labeller = labeller(param = c("b0" = "beta[0]",
                                           "b1" = "beta[1]",
                                           "b2" = "beta[2]"),
                                 .default = label_parsed)) +
  geom_hline(data = d.true.beta, aes(yintercept = true)) +
  theme_bw() +
  ylab("Mean estimate") +
  xlab("Grain size") +
  scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
  geom_point(data = d.summary.bias.beta, 
             aes(x = grain, y = avg, color = model),
             size = 2,
             position = position_dodge(width = 0.3)) +
  geom_errorbar(data = d.summary.bias.beta, 
                aes(x = grain, ymin = avg - sd, ymax = avg + sd, color = model),
                width = 0.2,
                position = position_dodge(width = 0.3)) +
  theme(legend.position = "top",
        text = element_text(size = 8)) +
  scale_color_manual(values = col.mod)


## BEM error estimates and true error approximation ----
d_sd.x |> 
  ggplot() +
  geom_point(aes(x = grain, y = sd.x), 
             color = col.mod[1], 
             alpha = 0.2,position = position_jitter(width = .1,
                                                    height = 0)) +
  ylab(expression(hat(sigma[X]))) +
  xlab("Grain size") +
  theme_bw() +
  geom_point(data = d.mean_sd.x, 
             aes(x = grain, y = avg),
             color = col.mod[1]) +
  geom_errorbar(data = d.mean_sd.x, 
                aes(x = grain, ymin = avg - sd, ymax = avg + sd),
                width = 0.1,
                color = col.mod[1]) +
  geom_point(data = d.true_sd.x,
             aes(x = grain, y = sd.x),
             color = "forestgreen",
             shape = 18,
             size = 4) +
  scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
  theme(legend.position = "none",
        text = element_text(size = 8)) 

## Estimates of spatial parameters (range, var) by spGLM ----
d_sp.x |> 
  ggplot() +
  geom_point(aes(x = grain, y = estimate), 
             color = col.mod[3], 
             alpha = 0.2,position = position_jitter(width = .1,
                                                    height = 0)) +
  facet_wrap(.~ param, scales = "free") +
  ylab("Estimate") +
  xlab("Grain size") +
  theme_bw() +
  geom_point(data = d.mean_sp.x, 
             aes(x = grain, y = avg),
             color = col.mod[3]) +
  geom_errorbar(data = d.mean_sp.x, 
                aes(x = grain, ymin = avg - sd, ymax = avg + sd),
                width = 0.1,
                color = col.mod[3]) +
  scale_x_discrete(labels = c("c1" = "Fine", "c2" = "Medium", "c3" = "Coarse")) +
  theme(legend.position = "none",
        text = element_text(size = 8)) 



#### Schema types of misaligned data ----
col.x <- "deeppink3"
col.y <- "black"

## Aligned case ##
sple <- sample(1:nrow(data), size = 200)
dt <- data |> 
  dplyr::slice(sple)

data |> 
  dplyr::select(lon, lat, x_fine) |> 
  ggplot(aes(x = lon, y = lat, fill = x_fine)) +
  geom_raster() +
  theme_void() +
  geom_point(data = dt, color = col.x, size = .8) +
  theme(legend.position = "none") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") -> gg.x1

data |> 
  dplyr::select(lon, lat, y_fine) |> 
  dplyr::mutate("y_fine" = as.factor(y_fine)) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_tile(aes(fill = y_fine)) +
  theme_void() +
  geom_point(data = dt, color = col.y, size = .8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("white", "forestgreen")) -> gg.y1

predd1 <- data.frame("x_fine" = seq(-3, 3, length.out = 200))
mod1 <- glm(y_fine ~ x_fine + I(x_fine^2), data = dt, family = binomial)
pred1 <- predict(mod1, predd1, type = "response")

dt |> 
  ggplot() +
  geom_jitter(aes(x = x_fine, y = y_fine), 
              alpha = .4,
             width = 0,
             height = .01) +
  theme_classic() +
  geom_line(aes(x = predd1$x_fine, y = pred1)) +
  theme(legend.position = "none",         
        axis.title = element_blank()) -> gg.src1

gg.src1 +
  geom_point(aes(x = x_fine, y = psi_fine),
             alpha = .4,
             col = "forestgreen")
  
## Point-to-point ##
sple.x <- sample(1:nrow(data), size = 200)
dt.x <- data |> 
  dplyr::slice(sple.x)
dt$x.err <- terra::extract(x_fine, 
                           dt[, c("lon", "lat")],
                           ID = FALSE)[, 1]

data |> 
  dplyr::select(lon, lat, x_fine) |> 
  ggplot(aes(x = lon, y = lat, fill = x_fine)) +
  geom_raster() +
  theme_void() +
  geom_point(data = dt, aes(x = lon, y = lat), size = .8, alpha = .4, color = col.y) +
  geom_point(data = dt.x, aes(x = lon, y = lat), size = .8, color = col.x) +
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,.5,"cm")) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") -> gg.x2

data |> 
  dplyr::select(lon, lat, y_fine) |> 
  dplyr::mutate("y_fine" = as.factor(y_fine)) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_tile(aes(fill = y_fine)) +
  theme_void() +
  geom_point(data = dt, aes(x = lon, y = lat), size = .8, color = col.y) +
  geom_point(data = dt.x, aes(x = lon, y = lat), size = .8, color = col.x, alpha = .4) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("white", "forestgreen")) -> gg.y2

predd <- data.frame("x.err" = seq(-3, 3, length.out = 200))
mod2 <- glm(y_fine ~ x.err + I(x.err^2), data = dt, family = binomial)
pred2 <- predict(mod2, predd, type = "response")
dt.p <- data.frame("yp" = c(pred1, pred2),
                   "x" = seq(-3, 3, length.out = 200),
                   "type" = rep(c("x_fine", "x.err"), each = 200))

dt |> 
  tidyr::pivot_longer(c(x_fine, x.err), names_to = "type", values_to = "x") |> 
  ggplot(aes( color = type)) +
  geom_jitter(aes(x = x, y = y_fine), 
             #color = "deeppink3",
              alpha = .4,
             position = ggstance::position_jitterdodgev(jitter.height = .02, 
                                                        dodge.height = .05)) +
  theme_classic() +
  scale_color_manual(values = c("black", "deeppink3")) +
  geom_line(data = dt.p, aes(x = x, y = yp)) +
  theme(legend.position = "none",         
        axis.title = element_blank()) -> gg.src2

## Point-to-area ##
psi_fine <- terra::rast(nrow = 1500, ncols = 1500, xmin = -250, xmax = 1250, ymin = -250, ymax = 1250)
terra::values(psi_fine) <- plogis(b0 + b1*terra::values(x_fine) + b2*(terra::values(x_fine)^2))
y_fine <- psi_fine
terra::values(y_fine) <- rbinom(length(terra::values(psi_fine)), 1, terra::values(psi_fine))
y_coarse <- terra::aggregate(y_fine, fact = 50, fun = function(x)max(sample(x, 1)))
y_c <- terra::extract(y_coarse, 
                      data[, c("lon", "lat")],
                      ID = FALSE)[, 1]
names(y_c) <- paste0("y_c")
data$y_c <- y_c

dt <- data |> 
  dplyr::slice(sple)

dt_coarse <- terra::extract(y_coarse,
                            dt[, c("lon", "lat")],
                            ID = FALSE,
                            xy = TRUE) |> 
  as.data.frame() |> 
  dplyr::mutate("lyr.1" = as.factor(lyr.1))
colnames(dt_coarse) <- c("y_c", "lon", "lat")

data |> 
  dplyr::select(lon, lat, x_fine) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = x_fine)) +
  theme_void() +
  geom_point(data = dt, color = col.x, size = .8) +
  #annotate("text",
  #         x = 500, y = 500, 
  #         label = expression(X^{EG}),
  #         size = .84) +
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,.5,"cm")) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") -> gg.x3

data |> 
  dplyr::select(lon, lat, y_c) |> 
  dplyr::mutate("y_c" = as.factor(y_c)) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_raster(aes(fill = y_c)) +
  geom_tile(data = dt_coarse, aes(fill = y_c), colour= col.y, linewidth = 0.5) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("white", "forestgreen")) -> gg.y3


predd <- data.frame("x_fine" = seq(-3, 3, length.out = 200))
mod3 <- glm(y_c ~ x_fine + I(x_fine^2), data = dt, family = binomial)
pred3 <- predict(mod3, predd, type = "response")
dt.p <- data.frame("yp" = c(pred1, pred3),
                   "x" = seq(-3, 3, length.out = 200),
                   "type" = rep(c("y_fine", "y_c"), each = 200))

dt |> 
  tidyr::pivot_longer(c(y_fine, y_c), names_to = "type", values_to = "y") |> 
  ggplot(aes( color = type)) +
  geom_jitter(aes(x = x_fine, y = y), 
              alpha = .4,
              position = ggstance::position_jitterdodgev(jitter.height = .02, 
                                                         dodge.height = .05)) +
  theme_classic() +
  scale_color_manual(values = c("deeppink3", "black")) +
  geom_line(data = dt.p, aes(x = x, y = yp)) +
  theme(legend.position = "none",         
        axis.title = element_blank()) -> gg.src3


## Area-to-point ##
dt_coarse <- terra::extract(x_coarse[[3]],
                            dt[, c("lon", "lat")],
                            ID = FALSE,
                            xy = TRUE) |> 
  as.data.frame() 
colnames(dt_coarse) <- c("x_c3", "lon", "lat")

data |> 
  dplyr::select(lon, lat, x_c3) |> 
  ggplot(aes(x = lon, y = lat, fill = x_c3)) +
  geom_raster(aes(fill = x_c3)) +
  geom_tile(data = dt_coarse, aes(fill = x_c3), colour=col.x, linewidth = 0.5) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,.5,"cm"))  +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") -> gg.x4

data |> 
  dplyr::select(lon, lat, y_fine) |> 
  dplyr::mutate("y_fine" = as.factor(y_fine)) |> 
  ggplot(aes(x = lon, y = lat)) +
  geom_tile(aes(fill = y_fine)) +
  theme_void() +
  geom_point(data = dt, color = col.y, size = .8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("white", "forestgreen")) -> gg.y4

predd <- data.frame("x_c2" = seq(-3, 3, length.out = 200))
mod4 <- glm(y_fine ~ x_c2 + I(x_c2^2), data = dt, family = binomial)
pred4 <- predict(mod4, predd, type = "response")
dt.p <- data.frame("yp" = c(pred1, pred4),
                   "x" = seq(-3, 3, length.out = 200),
                   "type" = rep(c("x_fine", "x_c2"), each = 200))

dt |> 
  tidyr::pivot_longer(c(x_fine, x_c2), names_to = "type", values_to = "x") |> 
  ggplot(aes( color = type)) +
  geom_jitter(aes(x = x, y = y_fine), 
              alpha = .4,
              position = ggstance::position_jitterdodgev(jitter.height = .02, 
                                                         dodge.height = .05)) +
  theme_classic() +
  scale_color_manual(values = c("deeppink3", "black")) +
  geom_line(data = dt.p, aes(x = x, y = yp)) +
  theme(legend.position = "none",
        axis.title = element_blank()) -> gg.src4

gg.misaligned1 <- cowplot::plot_grid(gg.x2, gg.y2, ncol = 2)
gg.misaligned2 <- cowplot::plot_grid(gg.x3, gg.y3, ncol = 2)
gg.misaligned3 <- cowplot::plot_grid(gg.x4, gg.y4, ncol = 2)

gg.misaligned <- cowplot::plot_grid(gg.misaligned1, gg.misaligned2, gg.misaligned3,
                                    nrow = 3,
                                    labels = c("a)", "b)", "c)"),
                                    label_size = 10)

#ggsave(here::here("results", "img", "spatial_misalignment_scheme.pdf"),
#       plot = gg.misaligned,
#       width = 1476,
#       height = 1204,
#       units = "px")
cowplot::plot_grid(#gg.x1, gg.y1, gg.src1, 
                   gg.x2, gg.y2, #gg.src2, 
                   gg.x3, gg.y3, #gg.src3, 
                   gg.x4, gg.y4, #gg.src4,
                   ncol = 2)
