# This script simulates data for adaptive dosing scenarios.

library(mipdtrial)
library(PKPDsim)
library(dplyr)


# Shared setup ------------------------------------------------------------
pop_params <- list(CL = 8.66, V = 100)
omega <- c(0.5, 0, 0.3)^2

mod <- PKPDsim::new_ode_model(
  code = "dAdt[0] = -(CLi/Vi)*A[0]; ",
  pk_code = "CLi = CL; Vi = V",
  parameters = c("CL", "V"),
  declare_variables = c("CLi", "Vi"),
  obs = list(cmt = 1, scale = "Vi"),
  dose = list(cmt = 1, bioav = 1),
  cpp = FALSE,
  iiv = list(CL = 0.5, V = 0.3),
  omega_matrix = omega,
  ruv = list(prop = 0.2)
)

model_design <- create_model_design(
  model = mod,
  parameters = pop_params,
  omega_matrix = omega,
  ruv = list(prop = 0.20, add = 0)
)

dat <- data.frame(ID = 1:200, WT = 70)

target_design <- create_target_design(
  targettype = "trough", targetvalue = 12.5,
  at = 14, anchor = "day"
)


# VPC simulation helpers -------------------------------------------------------

compute_pop_pred <- function(final_reg, t_obs) {
  pop_pred <- data.frame(id = numeric(0), t = numeric(0), pred = numeric(0))
  for (i in dat$ID) {
    doses <- final_reg[final_reg$id == i, ]
    reg <- new_regimen(
      time = doses$dose_times, 
      amt = doses$dose_amts, 
      t_inf = 0
    )
    res <- PKPDsim::sim(
      mod, 
      parameters = pop_params, 
      regimen = reg, 
      only_obs = TRUE, t_obs = t_obs
    )
    res$id <- i
    pop_pred <- bind_rows(pop_pred, select(res, id, t, pred = y))
  }
  pop_pred
}

run_vpc_sims <- function(final_reg, t_obs, n_vpc_sim = 1000, seed = 7) {
  set.seed(seed)
  n_obs <- length(t_obs)
  n_pts <- length(dat$ID)
  sim_df <- data.frame(
    id        = integer(n_pts * n_vpc_sim * n_obs),
    t         = numeric(n_pts * n_vpc_sim * n_obs),
    y         = numeric(n_pts * n_vpc_sim * n_obs),
    replicate = integer(n_pts * n_vpc_sim * n_obs)
  )
  k <- 1L
  for (i in dat$ID) {
    dos <- final_reg[final_reg$id == i, ]
    reg <- new_regimen(time = dos$dose_times, amt = dos$dose_amts, t_inf = 0)
    for (j in seq_len(n_vpc_sim)) {
      res <- PKPDsim::sim(
        mod, parameters = pop_params, regimen = reg, only_obs = TRUE,
        t_obs = t_obs, omega = omega, res_var = list(prop = 0.20, add = 0)
      )
      idx <- k:(k + n_obs - 1L)
      sim_df[idx, "id"]        <- i
      sim_df[idx, "t"]         <- res$t
      sim_df[idx, "y"]         <- res$y
      sim_df[idx, "replicate"] <- j
      k <- k + n_obs
    }
  }
  sim_df
}


# Case 0: Adaptive dosing with variation in dosing quantities --------------
# (Bergstrand et al. The AAPS Journal (2010)). Briefly, 200 simulated patients
# are given once-daily bolus doses, beginning with 10,000 U and subsequently
# adjusted each day based on measured trough concentrations. Troughs are sampled
# 23.5 hours after each dose, and individual PK parameters are estimated using
# MAP Bayesian estimation. Doses are adapted to achieve a target trough of
# 12.5 U/L. The popPK model is a one-compartment intravenous model with a
# clearance of 8.66 L/h (IIV of 50%), a volume of distribution of 100L (IIV of
# 30%), and a proportional error of 20%.

t_obs_case0 <- seq(23.5, by = 24, length.out = 6)

design_case0 <- create_trial_design(
  sampling_design = create_sampling_design(time = t_obs_case0),
  target_design = target_design,
  regimen_update_design = create_regimen_update_design(
    at = 2:6, anchor = "dose", update_type = "dose",
    dose_optimization_method = map_adjust_dose
  ),
  initial_regimen_design = create_initial_regimen_design(
    method = function(...) {
      new_regimen(n = 24, amt = 10000, interval = 24, t_inf = 0, type = "bolus")
    }
  ),
  sim_design = model_design,
  est_design = model_design
)

study0 <- run_trial(
  data = dat, 
  design = design_case0, 
  cov_mapping = list(), progress = FALSE, seed = 9
)

final_reg0 <- study0$final_reg
pop_pred0  <- compute_pop_pred(final_reg0, t_obs_case0)
vpc_sim0   <- run_vpc_sims(final_reg0, t_obs_case0)

saveRDS(
  left_join(study0$tdms[, c("id", "t", "y")], pop_pred0),
  "data/case0_obs.rds"
)
saveRDS(left_join(vpc_sim0, pop_pred0), "data/case0_sim.rds")


# Case 1: Adaptive dosing with variation in dosing intervals ---------------
# Case 0 but now the dosing interval is adapted while keeping dose quantity
# fixed. Patients are started on a regimen of 4000 U every 12 hours. Drug levels
# are collected at 23.5 hours and then every subsequent 48 hours, and the
# interval is adjusted to achieve a trough concentration of 12.5 U/L. Allowed
# intervals are 4, 6, 8, 12, 24, and 48 hours. This design ensures uniform
# sampling time across all patients and equal numbers of samples per patient,
# isolating the dosing interval as the only source of variation.

t_obs_case1 <- seq(23.5, by = 48, length.out = 6)

run_case1 <- function(update_design, seed = 9, vpc_seed = 7, n_vpc_sim = 1000) {
  design <- create_trial_design(
    sampling_design = create_sampling_design(time = t_obs_case1),
    target_design = target_design,
    regimen_update_design = update_design,
    initial_regimen_design = create_initial_regimen_design(
      method = function(...) {
        new_regimen(
          n = 300, amt = 4000, interval = 12,
          t_inf = 0, type = "bolus"
        )
      }
    ),
    sim_design = model_design,
    est_design = model_design,
    eval_design = create_eval_design(
      evaltype = "trough",
      time = seq(23.99, by = 48, length.out = 6)
    )
  )
  study <- run_trial(
    data = dat, design = design, cov_mapping = list(), 
    progress = FALSE, seed = seed
  )
  final_reg <- study$final_reg
  pop_pred  <- compute_pop_pred(final_reg, t_obs_case1)
  vpc_sims  <- run_vpc_sims(
    final_reg, t_obs_case1, n_vpc_sim = n_vpc_sim, seed = vpc_seed
  )
  list(
    obs = left_join(study$tdms[, c("id", "t", "y")], pop_pred),
    sim = left_join(vpc_sims, pop_pred)
  )
}

dose_update_design_case1 <- create_regimen_update_design(
  at = c(2, 4, 6, 8, 10), anchor = "day", update_type = "interval",
  dose_optimization_method = map_adjust_interval, # allow interval to vary
  grid = c(4, 6, 8, 12, 24, 48)                   # allowable intervals
)

case1 <- run_case1(dose_update_design_case1)
saveRDS(case1$obs, "data/case1_obs.rds")
saveRDS(case1$sim, "data/case1_sim.rds")


# Case 2: Adaptive dosing (quantity) with drop-out -------------------------
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
  inner_join(select(obs2, id, t))

saveRDS(obs2, "data/case2_obs.rds")
saveRDS(sim2, "data/case2_sim.rds")

# Case 3: Flat dosing with sample time variation -------------------------------
# Here we repeat Case 0, except we will not use dose adaptation. Instead, we
# will vary sampling based on CL estimates. Patients with low clearance (<50% 
# of population estimate) will provide both a peak and a trough on the next 
# dosing interval. All other patients will provide only a trough.

reg <- new_regimen(
  n = 5, amt = 10000, interval = 24, t_inf = 0, type = "bolus"
)
CL_thr <- pop_params$CL * 0.50
obs <- list()

set.seed(12)
for (i in dat$ID) {
  # first sample is a trough
  res <- sim(
    ode = mod,
    parameters = pop_params,
    regimen = reg,
    omega = omega,
    res_var = list(prop = 0.2),
    only_obs = TRUE,
    t_obs = t_obs_case0[1],
    output_include = list(variables = TRUE, parameters = TRUE)
  )
  tdms <- res[,c("t", "y")]
  iPK <- as.list(res[,names(pop_params)])
  for (d in seq_along(reg$dose_times)[-1]) {
    # Using available levels, perform a MAP Bayesian estimation
    fit <- simulate_fit(
      est_mod = mod, 
      parameters = pop_params, ruv = list(prop = 0.2), omega = omega,
      covariates = list(), regimen = reg,
      tdms = tdms
    )$parameters
    # if the CLi estimate is below (TVCL-1SD), next interval collect 2 samples
    # and otherwise collect a trough. Use individual PK for simulation, not
    # estimated parameters.
    if (fit$CL < CL_thr) {
      t_obs <- max(tdms$t) + c(1, 24)
    } else {
      t_obs <- max(tdms$t) + 24
    }
    res <- sim(
      ode = mod,
      parameters = iPK,
      regimen = reg,
      res_var = list(prop = 0.2),
      only_obs = TRUE,
      t_obs = t_obs
    )
    tdms <- rbind(tdms, res[,c("t", "y")])
  }
  tdms$CLi <- iPK$CL
  obs[[i]] <- tdms
}
obs <- bind_rows(obs, .id = "id")

# compute population predictions
pop_preds <- list()
for (i in unique(obs$id)) {
  t_obs <- obs$t[obs$id == i]
  pop_preds[[i]] <- sim(
    ode = mod,
    parameters = pop_params,
    regimen = reg,
    only_obs = TRUE,
    t_obs = t_obs
  )[, c("t", "y")]
}
pop_preds <- bind_rows(pop_preds, .id = "id") |>
  rename(pred = y)


# compute simulations for VPC
sims <- list()
n_vpc_sim <- 1000
set.seed(120)
for (i in unique(obs$id)) {
  t_obs <- obs$t[obs$id == i]
  sim_i <- list()
  for (rep in 1:n_vpc_sim) {
    sim_i[[rep]] <- sim(
      ode = mod,
      parameters = pop_params,
      regimen = reg,
      omega = omega,
      res_var = list(prop = 0.2),
      only_obs = TRUE,
      t_obs = t_obs
    )[, c("t", "y")]
  }
  sims[[i]] <- bind_rows(sim_i, .id = "replicate")
}
sims <- bind_rows(sims, .id = "id")
obs3 <- left_join(obs, pop_preds, by = join_by(id, t))
sim3 <- left_join(sims, pop_preds, by = join_by(id, t))

saveRDS(obs3, "data/case3_obs.rds")
saveRDS(sim3, "data/case3_sim.rds")

# Case 1X: Adaptive dosing with variation in dosing intervals, N = 10 -----
# Here we repeat Case 1 exactly, but 10 times, allowing the assessment of
# variability between simulated trials.

out <- lapply(
  1:10,
  function(seed) {
    run_case1(dose_update_design_case1, seed = seed, vpc_seed = seed * 100)
  }
)
saveRDS(out, "data/case1X.rds")

# Case 0X: Adaptive dosing with variation in dosing quantities, dose rounding --
# Here we repeat Case 0, but require doses to be given in multiples of 2000 U 
# instead of 1 U (default behaviour in mipdtrial)

design_case0X <- create_trial_design(
  sampling_design = create_sampling_design(time = t_obs_case0),
  target_design = target_design,
  regimen_update_design = create_regimen_update_design(
    at = 2:6, anchor = "dose", update_type = "dose",
    dose_optimization_method = map_adjust_dose,
    settings = list(dose_resolution = 2000)
  ),
  initial_regimen_design = create_initial_regimen_design(
    method = function(...) {
      new_regimen(n = 24, amt = 10000, interval = 24, t_inf = 0, type = "bolus")
    }
  ),
  sim_design = model_design,
  est_design = model_design
)

study0X <- run_trial(
  data = dat, 
  design = design_case0X, 
  cov_mapping = list(), progress = FALSE, seed = 9
)

final_reg0X <- study0X$final_reg
pop_pred0X  <- compute_pop_pred(final_reg0X, t_obs_case0)
vpc_sim0X   <- run_vpc_sims(final_reg0X, t_obs_case0)

saveRDS(
  left_join(study0X$tdms[, c("id", "t", "y")], pop_pred0X), 
  "data/case0X_obs.rds"
)
saveRDS(
  left_join(vpc_sim0X, pop_pred0X),
  "data/case0X_sim.rds"
)
