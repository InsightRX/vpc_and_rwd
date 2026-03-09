# Script to create Figure 1

library(vpc)
library(ggplot2)
library(dplyr)
library(patchwork)

# read in data sets ------------------------------------------------------------
# case 0:
obs_0 <- readRDS("data/case0_obs.rds")
sim_0 <- readRDS("data/case0_sim.rds")

# case 1
obs_1 <- readRDS("data/case1_obs.rds")
sim_1 <- readRDS("data/case1_sim.rds")

# case 2
obs_2 <- readRDS("data/case2_obs.rds")
sim_2 <- readRDS("data/case2_sim.rds")

# case 3
obs_3 <- readRDS("data/case3_obs.rds")
sim_3 <- readRDS("data/case3_sim.rds")


# case 2b
obs_2b <- readRDS("data/case2b_obs.rds") |>
  rename(id = ID, pred = preds)
sim_2b <- readRDS("data/case2b_sim.rds") |>
  rename(id = ID, pred = preds)
sim_2bc <- readRDS("data/case2b_sim_censor.rds") |>
  rename(id = ID, pred = preds)
obs_2bf <- readRDS("data/case2b_obs_full.rds") |>
  rename(id = ID, pred = preds)
sim_2bf <- readRDS("data/case2b_sim_full.rds") |>
  rename(id = ID, pred = preds)

# plotting convenience function ------------------------------------------------
plot_vpc <- function(obs, sim, pc, idv) {
  pi <- vpc::vpc(
    obs = obs, 
    sim = sim, 
    obs_cols = list(
      dv = "y", idv = idv, id = "id", pred = "pred"
    ), 
    sim_cols = list(
      dv = "y", idv = idv, id = "id", pred = "pred", sim = "replicate"
    ),
    pred_corr = pc,
    show = list(obs_dv = TRUE),
    vpcdb = FALSE
  ) 
  pi$layers[[7]]$aes_params$shape <- 16 # change obs points to filled circles
  pi$layers[[7]]$aes_params$alpha <- 0.3
  pi + 
    labs(
      x = ifelse(idv == "t", "Time (h)", "Time after dose (h)"),
      y = "U/L",
      title = ifelse(pc, "pcVPC", "VPC")
    ) +
    theme_bw() +
    theme(
      plot.tag.position = c(0.01, 0.99)
    ) +
    scale_y_log10() 
}

# figure 1 ---------------------------------------------------------------------
p0a <- plot_vpc(obs_0, sim_0[sim_0$t < 120,], FALSE, "t") +
  labs(tag = "(a)", title = "Case 0, VPC")
p0b <- plot_vpc(obs_0, sim_0[sim_0$t < 120,], TRUE, "t") +
  labs(tag = "(b)", title = "Case 0, pcVPC")
p1a <- plot_vpc(obs_1, sim_1[sim_1$t < 216,], FALSE, "t") +
  labs(tag = "(c)", title = "Case 1, VPC")
p1b <- plot_vpc(obs_1, sim_1[sim_1$t < 216,], TRUE, "t") +
  labs(tag = "(d)", title = "Case 1, pcVPC")
p2a <- plot_vpc(obs_2, sim_2[sim_2$t < 120,], FALSE, "t") +
  labs(tag = "(e)", title = "Case 2, VPC")
p2b <- plot_vpc(obs_2, sim_2[sim_2$t < 120,], TRUE, "t") +
  labs(tag = "(f)", title = "Case 2, pcVPC")
p3a <- plot_vpc(obs_3, sim_3, FALSE, "t") +
  labs(tag = "(g)", title = "Case 3, VPC")
p3b <- plot_vpc(obs_3, sim_3, TRUE, "t") +
  labs(tag = "(h)", title = "Case 3, pcVPC")

p0a + p0b + p1a + p1b + p2a + p2b + p3a + p3b +
  plot_layout(axes = "collect_x", ncol = 2)
  
ggsave("figures/all_plots.png", scale = 1.5, width = 6, height = 7)

# supplemental figures ---------------------------------------------------------
# methotrexate case study
p2ba <- plot_vpc(obs_2bf, sim_2bf, FALSE, "TAD") +
  labs(
    y = "MTX (µmol/L.hr)", tag = "(a)", title = "Case 2b, VPC, no censoring"
  )
p2bb <- plot_vpc(obs_2bf, sim_2bf, TRUE, "TAD") +
  labs(
    y = "MTX (µmol/L.hr)", tag = "(b)", title = "Case 2b, pcVPC, no censoring"
  )
p2bc <- plot_vpc(obs_2b, sim_2b, FALSE, "TAD") +
  labs(y = "MTX (µmol/L.hr)", tag = "(c)", title = "Case 2b, VPC")
p2bd <- plot_vpc(obs_2b, sim_2b, TRUE, "TAD") +
  labs(y = "MTX (µmol/L.hr)", tag = "(d)", title = "Case 2b, pcVPC")
p2be <- plot_vpc(obs_2b, sim_2bc, FALSE, "TAD") +
  labs(
    y = "MTX (µmol/L.hr)", tag = "(e)", title = "Case 2b, VPC, dropout model"
  )
p2bf <- plot_vpc(obs_2b, sim_2bc, TRUE, "TAD") +
  labs(
    y = "MTX (µmol/L.hr)", tag = "(f)", title = "Case 2b, pcVPC, dropout model"
  )

p2ba + p2bb + p2bc + p2bd + p2be + p2bf +
  plot_layout(axes = "collect_x", ncol = 2)
ggsave(
  "figures/case2b.png", 
  scale = 1.5, width = 6, height = 6
)

# Case 0X: VPCs for Case 0 with dose rounding applied --------------------------
obs_0X <- readRDS("data/case0X_obs.rds")
sim_0X <- readRDS("data/case0X_sim.rds")
p0Xa <- plot_vpc(obs_0X, sim_0X[sim_0X$t < 120,], FALSE, "t") +
  labs(tag = "(a)", title = "Case 0, VPC")
p0Xb <- plot_vpc(obs_0X, sim_0X[sim_0X$t < 120,], TRUE, "t") +
  labs(tag = "(b)", title = "Case 0, pcVPC")
p0Xa + p0Xb
ggsave("figures/case0X.png", scale = 1.5, width = 6, height = 2.5)

# Case 1X: VPCs for Case 1 simulated with different seeds ----------------------
folds <- readRDS("data/case1X.rds")
plot_fold_n <- function(fold, pc) {
  title_string <- paste0(
    "Case 1, ", ifelse(pc, "pc", ""), "VPC, Seed = ", fold
  )
  plot_vpc(
    folds[[fold]]$obs, 
    folds[[fold]]$sim[folds[[fold]]$sim$t < 216,], 
    pc, "t"
  ) +
    labs(title = title_string)
}
case1_plots <- c(
  lapply(1:10, plot_fold_n, FALSE), # VPC
  lapply(1:10, plot_fold_n, TRUE)   # pcVPC
)
p_case1X  <- wrap_plots(case1_plots, ncol = 4, byrow = FALSE)

ggsave("figures/case1X.png", scale = 2, width = 6, height = 7)
