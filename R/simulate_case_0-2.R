# This script simulates data for adaptive dosing scenarios.

library(mipdtrial)
library(PKPDsim)
library(dplyr)

# - Case 0: Adaptive dosing with variation in dosing ---------------------------
# (Bergstrand et al. The AAPS Journal (2010)). Briefly, 200 simulated patients 
# are given once-daily bolus doses, beginning with 10,000 U and subsequently 
# adjusted each day based on measured trough concentrations. Troughs are sampled
# 23.5 hours after each dose, and individual PK parameters are estimated using
# MAP Bayesian estimation. Doses are adapted to achieve a target trough of 
# 12.5 U/L. The popPK model is a one-compartment intravenous model with a 
# clearance of 8.66 L/h (IIV of 50%), a volume of distribution of 100L (IIV of 
# 30%), and a proportional error of 20%.


# define the model
mod <- PKPDsim::new_ode_model(
  code = "dAdt[0] = -(CLi/Vi)*A[0]; ",
  pk_code = "CLi = CL; Vi = V",
  parameters = c("CL", "V"),
  declare_variables = c("CLi", "Vi"),
  obs = list(cmt = 1, scale = "Vi"),
  dose = list(cmt = 1, bioav = 1),
  cpp = FALSE,
  iiv = list(CL = 0.5, V = 0.3),
  omega_matrix = c(0.5, 0, 0.3) ** 2,
  ruv = list(prop = 0.2)
)

model_design <- create_model_design(
  model = mod,
  parameters = list(CL = 8.66, V = 100),
  omega_matrix = c(0.5, 0.0225, 0.3) ** 2,
  ruv = list(prop = 0.20, add = 0)
)

# Set up dose adaptation trial. Initial dose: 10,000 U, trough concentration
# plus dose adaptation once per day to target a trough of 12.5 U/L. 
t_obs <- seq(23.5, by = 24, length.out = 6)

tdm_design <- create_sampling_design(
  time = t_obs
)
target_design <- create_target_design(
  targettype = "trough", targetvalue = 12.5,
  at = 14, anchor = "day"
)
dose_update_design <- create_regimen_update_design(
  at = 2:6, anchor = "dose", update_type = "dose",
  dose_optimization_method = map_adjust_dose
)
flat_dose <- function(...) {
  new_regimen(n = 24, amt = 10000, interval = 24, t_inf = 0, type = "bolus")
}
reg_initial_design <- create_initial_regimen_design(
  method = flat_dose 
)

# set up data: 200 adult patients
dat <- data.frame(ID = 1:200, WT = 70)

# simulate the trial:
design <- create_trial_design(
  sampling_design = tdm_design,
  target_design = target_design,
  regimen_update_design = dose_update_design,
  initial_regimen_design = reg_initial_design,
  sim_design = model_design, 
  est_design = model_design
)
study <- run_trial(
  data = dat,
  design = design,
  cov_mapping = list(),
  progress = FALSE,
  seed = 9
)

# VPC: observed data
obs <- study$tdms[, c("id", "t", "y")]

# Get population predictions given the final regimen
final_reg <- study$final_reg
pop_pred <- data.frame(id = numeric(0), t = numeric(0), pred = numeric(0))
for (i in dat$ID) {
  doses <- final_reg[final_reg$id == i,]
  reg <- new_regimen(time = doses$dose_times, amt = doses$dose_amts, t_inf = 0)
  res <- sim(
    mod,
    parameters = list(CL = 8.66, V = 100),
    regimen = reg,
    only_obs = TRUE,
    t_obs = t_obs
  )
  res$id = i
  pop_pred <- bind_rows(
    pop_pred,
    select(res, id, t, pred = y)
  )
}

# VPC: simulations
set.seed(7)
sim <- data.frame(
  id = numeric(0), t = numeric(0), y = numeric(0), replicate = numeric(0)
)
n_vpc_sim <- 1000
for (i in dat$ID) {
  doses <- final_reg[final_reg$id == i,]
  reg <- new_regimen(time = doses$dose_times, amt = doses$dose_amts, t_inf = 0)
  for (j in seq_len(n_vpc_sim)) {
    res <- sim(
      mod,
      parameters = list(CL = 8.66, V = 100),
      regimen = reg,
      only_obs = TRUE,
      t_obs = t_obs[1:5],
      omega = c(0.5, 0, 0.3) ** 2,
      res_var = list(prop = 0.20, add = 0)
    )
    res$id <- i
    res$replicate <- j
    sim <- bind_rows(
      sim,
      select(res, id, t, y, replicate)
    )
  } 
}

# add population predictions
obs_pc <- left_join(obs, pop_pred)
sim_pc <- left_join(sim, pop_pred)

saveRDS(obs_pc, "data/case0_obs.rds")
saveRDS(sim_pc, "data/case0_sim.rds")


# - Case 1: Adaptive dosing with variation in dosing intervals -----------------
# Case 0 but now the dosing interval is adapted while keeping dose quantity 
# fixed. Patients are started on a regimen of 4000 U every 12 hours. Drug levels
# are collected at 23.5 hours and then every subsequent 48 hours, and the 
# interval is adjusted to achieve a trough concentration of 12.5 U/L. Allowed 
# intervals are 4, 6, 8, 12, 24, and 48 hours. This design ensures uniform 
# sampling time across all patients and equal numbers of samples per patient, 
# isolating the dosing interval as the only source of variation.

t_obs <- seq(23.5, by = 48, length.out = 6)
flat_dose <- function(...) { # 4000 allows more variation in dosing intervals
  new_regimen(n = 300, amt = 4000, interval = 12, t_inf = 0, type = "bolus")
}
reg_initial_design <- create_initial_regimen_design(
  method = flat_dose 
)
dose_update_design <- create_regimen_update_design(
  at = c(2, 4, 6, 8, 10),
  anchor = "day",
  update_type = "interval",
  dose_optimization_method = map_adjust_interval, # allow interval to vary
  grid = c(4, 6, 8, 12, 24, 48) # allowable intervals
)
tdm_design <- create_sampling_design(
  time = t_obs
)
eval_design <- create_eval_design(
  evaltype = "trough",
  time = c(seq(23.99, by = 48, length.out = 6))
)

# simulate the trial:
design <- create_trial_design(
  sampling_design = tdm_design,
  target_design = target_design,
  regimen_update_design = dose_update_design,
  initial_regimen_design = reg_initial_design,
  sim_design = model_design, 
  est_design = model_design,
  eval_design = eval_design
)
study <- run_trial(
  data = dat,
  design = design,
  cov_mapping = list(),
  progress = FALSE,
  seed = 9
)

# VPC: observed data
obs <- study$tdms[, c("id", "t", "y")]

# VPC: population predictions
final_reg <- study$final_reg
pop_pred <- data.frame(id = numeric(0), t = numeric(0), pred = numeric(0))
for (i in dat$ID) {
  doses <- final_reg[final_reg$id == i,]
  reg <- new_regimen(time = doses$dose_times, amt = doses$dose_amts, t_inf = 0)
  res <- sim(
    mod,
    parameters = list(CL = 8.66, V = 100),
    regimen = reg,
    only_obs = TRUE,
    t_obs = t_obs # avoid bolus
  )
  res$id = i
  pop_pred <- bind_rows(
    pop_pred,
    select(res, id, t, pred = y)
  )
}

# VPC: simulations
set.seed(7)
sim <- data.frame(
  id = numeric(0), t = numeric(0), y = numeric(0), replicate = numeric(0)
)
n_vpc_sim <- 1000
for (i in dat$ID) {
  doses <- final_reg[final_reg$id == i,]
  reg <- new_regimen(time = doses$dose_times, amt = doses$dose_amts, t_inf = 0)
  for (j in seq_len(n_vpc_sim)) {
    res <- sim(
      mod,
      parameters = list(CL = 8.66, V = 100),
      regimen = reg,
      only_obs = TRUE,
      t_obs = t_obs, 
      omega = c(0.5, 0.0225, 0.3) ** 2,
      res_var = list(prop = 0.20, add = 0)
    )
    res$id <- i
    res$replicate <- j
    sim <- bind_rows(sim, res[, c("id", "t", "y", "replicate")])
  }
  
}

# add population predictions
obs_pc <- left_join(obs, pop_pred)
sim_pc <- left_join(sim, pop_pred)

# add time after dose
tad <- function(id, t, final_reg) {
  dt <- final_reg$dose_times[final_reg$id == id]
  t - dt[max(which(dt < t))]
}

saveRDS(obs_pc, "data/case1_obs.rds")
saveRDS(sim_pc, "data/case1_sim.rds")


# - Case 2: Adaptive dosing (quantity) with drop-out ---------------------------
# Here we can take the Case 0 data, but remove observations once a patient is
# at target, i.e., once an appropriate maintenance dose has been found. We will
# deem a patient as in-target if their trough concentration is within 20% of the
# target of 12.5 U/L (10-15 U/L).

obs0 <- readRDS("data/case0_obs.rds")
sim0 <- readRDS("data/case0_sim.rds")

obs2 <- obs0 |>
  group_by(id) |>
  mutate(
    reached_target = cumsum(between(y, 10, 15)),
    # once on target, stop sampling even if later samples are off target
    first_target = cumsum(reached_target)
  ) |>
  filter(first_target <= 1) |>
  mutate(n_tdm = n())
sim2 <- sim0 |>
  # only simulate levels that were observed
  inner_join(
    select(obs2, id, t)
  )

saveRDS(obs2, "data/case2_obs.rds")
saveRDS(sim2, "data/case2_sim.rds")
