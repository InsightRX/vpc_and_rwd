# This script simulates data for patient drop-out, using methotrexate as a
# case study.

# Load required packages (all on CRAN)
library(NHANES)
library(dplyr)
library(tidyr)
library(PKPDsim) 
library(clinPK)  

# - Case 2b: Drug concentration monitoring for patient discharge
# In high-dose methotrexate therapy, patients remain in the hospital until their
# drug concentration drops below a certain threshold, with regular monitoring 
# until that level is reached. Here, we use the model by Pauley et al. for
# simulations (Pauley et al., Cancer Chemother Pharmacol. (2013)). Covariate 
# values for individuals under the age of 18 years are sampled from the NHANES R 
# package (N = 200). Patients are dosed at 1 g/m2 infused over 24 hours every 
# 14 days for four cycles, without dose adjustment. Levels were collected at 
# 6, 23, and 42 hours after the start of infusion for all patients, and then 
# every 6 hours until the patient’s drug concentration drops below 0.1 µmol/L.

# define the model
mod <- PKPDsim::new_ode_model(
  code = "\
    KEi = KE * exp(kappa_KE)  \
    dAdt[0] = -KEi*A[0] - KCPi*A[0] + KPCi*A[1] \
    dAdt[1] = +KCPi*A[0] - KPCi*A[1] \
    dAdt[2] = A[0] / (V * BSA * (454.44/1000.0)) \
  ",
  pk_code = " \
    Vi = V * BSA; \
    KCPi = KCP; \
    KPCi = KPC; \
  ",
  parameters = c("KE", "V", "KCP", "KPC"),
  declare_variables = c("Vi", "KEi", "KCPi", "KPCi"),
  covariates = c("BSA"),
  obs = list(cmt = 1, scale = "Vi * (454.44/1000.0)"),
  dose = list(cmt = 1, bioav = 1),
  cpp = FALSE,
  iiv = list(KE = 0.314, V = 0.52, KCP = 0.625, KPC = 0.0345),
  omega_matrix = c(0.0987755, 0, 0.270907, 0, 0, 0.390625, 0, 0, 0, 0.0011934),
  ruv = list(prop = 0.3, add = 0.01),
  iov = list(
    cv = list(KE = 0.169),
    n_bins = 5,
    bins = c(0, 336, 672, 1008, 1344, 9999)
  )
)

# parametrize VPC / simulate
n_vpc_sim <- 1000

# data: NHANES
set.seed(10)
pts <- NHANES |>
  select(ID, Gender, Age, Height, Weight) |>
  filter(Age < 18) |>
  drop_na() |>
  group_by(ID) |>
  slice(1) |> # first row per ID only
  ungroup() |>
  mutate(BSA = clinPK::calc_bsa(weight = Weight, height = Height)$value) |>
  sample_n(200) # more reasonable data size

# tdms: 6, 23 and 42 h from the infusion start, then every 6 hours until
# concentration < 0.1 micromol/l
t_obs <- c(6, 23, 42, seq(48, 120, 6))
t_obs <- rep(t_obs, 4) + rep(c(0, 1, 2, 3) * 24 * 14, each = length(t_obs))

# regimen: four courses of 1 g/m2 infused over 24h every 14 days.
get_regimen <- function(bsa) {
  new_regimen(
    n = 4, t_inf = 24, type = "infusion", interval = 14 * 24,
    amt = 1000 * bsa
  )
}

pop_par <- list(
  # main parameters
  "KE" = 0.7, "V" = 9.03, "KCP" = 0.08, "KPC" = 0.11,
  # iov terms
  kappa_KE_1 = 0, kappa_KE_2 = 0, kappa_KE_3 = 0, kappa_KE_4 = 0, kappa_KE_5 = 0
)

# function for performing simulation
sim_mtx <- function(
  dat,
  t_obs,
  mod,
  f_reg = get_regimen,
  par = pop_par,
  iiv = c(0.0987755, 0, 0.270907, 0, 0, 0.390625, 0, 0, 0, 0.0011934),
  ruv = list(prop = 0.3, add = 0.01),
  iov_bins = c(0, 336, 672, 1008, 1344, 9999)
) {
  # output initialization:
  doses <- data.frame(
    ID = numeric(0),
    dose_times = numeric(0),
    dose_amts = numeric(0),
    t_inf = numeric(0)
  )
  obs <- data.frame(
    ID = numeric(0),
    t = numeric(0),
    y = numeric(0),
    TAD = numeric(0)
  )
  
  for (id in dat$ID) {
    # set up simulation
    pt <- dat[dat$ID == id,]
    reg <- get_regimen(pt$BSA)
    covs <- list(BSA = pt$BSA)

    # set up iov object
    if (!is.null(iiv)) {
      iov <- mipdtrial:::get_iov_specification (mod, par, iiv)
    } else {
      iov <- list(omega = NULL, omega_type = NULL)
    }
    
    
    # simulation
    res <- sim(
      ode = mod,
      parameters = par,
      regimen = reg,
      covariates = covs,
      omega = iov$omega,
      omega_type = iov$omega_type,
      res_var = ruv,
      only_obs = TRUE,
      t_obs = t_obs,
      iov_bins = iov_bins
    )
    
    # save results
    reg <- as.data.frame(reg)
    reg$ID <- id
    doses <- bind_rows(doses, reg[colnames(doses)])
    res$ID <- id
    res$TAD <- res$t %% 336 # doses are every 14 days = 336 hours
    obs <- bind_rows(obs, res[colnames(obs)])
  }
  list(doses = doses, obs = obs)
}

# when patients have a concentration below `thresh`, they can go home, so
# remove observations "collected" after they met this threshold
send_patients_home <- function(obs, thresh = 0.1) {
  # keep all obs > thresh and fixed samples
  keep1 <- filter(obs, TAD <= 42 | y > thresh)
  
  # also keep the FIRST observation per patient per interval <= 0.1
  keep2 <- obs |>
    group_by(ID) |>
    mutate(dose_nr = cumsum(TAD == 0)) |>
    group_by(ID, dose_nr) |>
    filter(y <= thresh & TAD > 0) |>
    slice(1) |>
    ungroup() |>
    select(-dose_nr) |>
    filter(TAD > 42) # if 42-hr time point already < 0.1, is included above
  
  # combine and return
  bind_rows(keep1, keep2) |>
    arrange(ID, t)
}


# --- generate "real" data -----------------------------------------------------
# we will create two sets of data, one with all collected observations, and one
# with only the first 4 TDM time points, which all patients provide. We will
# append with population predictions.
set.seed(7)
full_data <- sim_mtx(pts, t_obs, mod)
pop_preds <- sim_mtx(pts, t_obs, mod, iiv = NULL, ruv = NULL)$obs |> # pop preds
  select(ID, t, preds = y)
obs_data <- send_patients_home(full_data$obs, thresh = 0.1) |>
  left_join(pop_preds)
obs_data$y[obs_data$y < 0] <- 0 # LOQ
full_data$obs$y[full_data$obs$y < 0] <- 0

# ---- now run VPC simulations -------------------------------------------------
# initialize simulation output
full_sim <- data.frame(
  matrix(nrow = n_vpc_sim * nrow(full_data$obs), ncol = 5)
)
names(full_sim) <- c("ID", "t", "y", "TAD", "replicate")

set.seed(22)
for (i in seq_len(n_vpc_sim)) {
  obs <- sim_mtx(pts, t_obs, mod)$obs
  obs$replicate <- i

  idx <- ((i - 1) * nrow(obs) + 1):(i * nrow(obs))
  full_sim[idx, ] <- obs[, c("ID","t","y","TAD","replicate")]
}
full_sim$y[full_sim$y < 0] <- 0
sim_match_obs <- inner_join(full_sim, select(obs_data, ID, t, preds))
sim_censor <- send_patients_home(full_sim, thresh = 0.1) |>
  left_join(pop_preds)

# save data
saveRDS(left_join(full_sim, pop_preds), "data/case2b_sim_full.rds")
saveRDS(left_join(full_data$obs, pop_preds), "data/case2b_obs_full.rds")
saveRDS(obs_data, "data/case2b_obs.rds")
saveRDS(sim_match_obs, "data/case2b_sim.rds")
saveRDS(sim_censor, "data/case2b_sim_censor.rds")
