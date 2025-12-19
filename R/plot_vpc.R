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
obs_3 <- readRDS("data/case3_obs.rds") |>
  rename(id = ID, pred = preds)
sim_3 <- readRDS("data/case3_sim.rds") |>
  rename(id = ID, pred = preds)
sim_3c <- readRDS("data/case3_sim_censor.rds") |>
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


p0a <- plot_vpc(obs_0, sim_0, FALSE, "t") +
  labs(tag = "(a)", title = "Case 0, VPC")
p0b <- plot_vpc(obs_0, sim_0, TRUE, "t") +
  labs(tag = "(b)", title = "Case 0, pcVPC")
p1a <- plot_vpc(obs_1, sim_1, FALSE, "t") +
  labs(tag = "(c)", title = "Case 1, VPC")
p1b <- plot_vpc(obs_1, sim_1, TRUE, "t") +
  labs(tag = "(d)", title = "Case 1, pcVPC")
p2a <- plot_vpc(obs_2, sim_2, FALSE, "t") +
  labs(tag = "(e)", title = "Case 2, VPC")
p2b <- plot_vpc(obs_2, sim_2, TRUE, "t") +
  labs(tag = "(f)", title = "Case 2, pcVPC")
p3a <- plot_vpc(obs_3, sim_3, FALSE, "TAD") +
  labs(y = "MTX (µmol/L.hr)", tag = "(g)", title = "Case 3, VPC")
p3b <- plot_vpc(obs_3, sim_3, TRUE, "TAD") +
  labs(y = "MTX (µmol/L.hr)", tag = "(h)", title = "Case 3, pcVPC")
p3c <- plot_vpc(obs_3, sim_3c, FALSE, "TAD") +
  labs(
    y = "MTX (µmol/L.hr)", tag = "(i)", title = "Case 3, VPC, dropout model"
  )


r02 <- p0a + p0b + p1a + p1b + p2a + p2b + 
  plot_layout(axes = "collect_x", ncol = 2)
r3 <- p3a + p3b + p3c + plot_layout(axes = "collect_x")

r02 / r3 + plot_layout(heights = c(3, 1)) &
  theme(plot.title = element_text(size = 10))
  
ggsave("figures/all_plots.png", scale = 1.5, width = 6, height = 7)
