if(!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(npmv, MANOVA.RM, tidyverse, MASS, pbmcapply, rankFD, dplyr)

set.seed(2025)

# ===========================
# Parameter & Optionen
# ===========================
R <- 1000                  # Anzahl Simulationen pro Setting
alpha <- 0.05
n_cores <- parallel::detectCores() - 1
# Szenarienauswahl
run_scenarios <- c("type1")  

# Stichprobensettings
sample_settings <- list(
  extrasmall_bal   = list(N = 12,  treat = 6,   placeb = 6),
  extrasmall_unbal = list(N = 12,  treat = 9,  placeb = 3),
  small_bal   = list(N = 18,  treat = 9,   placeb = 9),
  small_unbal = list(N = 18,  treat = 12,  placeb = 6),
  mod1_bal     = list(N = 36,  treat = 18,  placeb = 18),
  mod1_unbal   = list(N = 36,  treat = 24,  placeb = 12),
  mod2_bal     = list(N = 72,  treat = 36,  placeb = 36),
  mod2_unbal   = list(N = 72,  treat = 48,  placeb = 24),
  mod3_bal   = list(N = 150, treat = 75,  placeb = 75),
  mod3_unbal = list(N = 150, treat = 100, placeb = 50),
  large_bal   = list(N = 360, treat = 180, placeb = 180),
  large_unbal = list(N = 360, treat = 240, placeb = 120)
)
run_sample_settings <- names(sample_settings)

# Neue Effekt-Settings
effect_settings <- list(
  eff80_20 = list(treat_resp = 0.8, placeb_resp = 0.2)
  #eff60_40 = list(treat_resp = 0.6, placeb_resp = 0.4)
)

# Shift-Settings nur für power2 relevant
shift_settings <- list(
  shiftA21 = list(mean_shift = 2, sd_shift = 1),
  shiftB151 = list(mean_shift = 1.5, sd_shift = 1)
)

# Debug-Optionen
show_data_example <- F
calc_means <- F

# ===========================
# Datengenerator
# ===========================
simulate_data <- function(N, t_size, p_size, scenario, treat_resp, placeb_resp, mean_shift, sd_shift) {
  sim <- data.frame(
    patientID = 1:N,
    treatment = factor(c(rep("Treatment", t_size), rep("Placebo", p_size))),
    treatment_n = c(rep(1, t_size), rep(0, p_size))
  )
  
  # wound_size_pre
  std <- 60
  sim$wound_size_pre <- round((rchisq(N, 4) / sqrt(8)) * std)
  sim$wound_size_pre[sim$wound_size_pre < 20] <- 20
  #sim$wound_size_pre <- round(rlnorm(N, meanlog=4, sdlog=1.2))
  # wound_size_post
  sim$wound_size_post <- sim$wound_size_pre
  if (scenario %in% c("power1", "power2")) {
    heal <- logical(N)
    heal[sim$treatment_n == 1] <- rbinom(sum(sim$treatment_n == 1), 1, treat_resp) == 1
    heal[sim$treatment_n == 0] <- rbinom(sum(sim$treatment_n == 0), 1, placeb_resp) == 1
    if (any(heal)) {
      sim$wound_size_post[heal] <- round(sim$wound_size_pre[heal] * 0.10 + rnorm(sum(heal), mean = 0, sd = 2))
    }
  } else {
    sim$wound_size_post <- round(sim$wound_size_pre + rnorm(N, mean = 0, sd = 2))
  }
  sim$wound_size_post[sim$wound_size_post < 0] <- 0
  
  # Pain & Itch baseline
  
  if (scenario == "power2") {
    
    Sigma_t <- matrix(c(9, 3.6,
                        3.6, 9), 2)
    
    Sigma_p <- matrix(c(16, 8,
                        8, 25), 2)
    
    pi_mat <- matrix(NA, N, 2)
    
    idx_t <- sim$treatment_n == 1
    idx_p <- sim$treatment_n == 0
    
    pi_mat[idx_t, ] <- MASS::mvrnorm(sum(idx_t),
                                     mu = c(5, 5),
                                     Sigma = Sigma_t)
    
    pi_mat[idx_p, ] <- MASS::mvrnorm(sum(idx_p),
                                     mu = c(5, 5),
                                     Sigma = Sigma_p)
    
  } else {
    
    # Kovarianzmatrizen Treatment vs Placebo
    Sigma_t <- matrix(c(
      9,  3.6,
      3.6, 9
    ), 2)
    
    Sigma_p <- matrix(c(16, 8,
                        8, 25), 2)
    
    pi_mat <- matrix(NA, N, 2)
    
    idx_t <- sim$treatment_n == 1
    idx_p <- sim$treatment_n == 0
    
    # WICHTIG:
    # gleiche Mittelwerte -> Type-I-H0 bleibt wahr
    pi_mat[idx_t, ] <- MASS::mvrnorm(
      n = sum(idx_t),
      mu = c(5, 5),
      Sigma = Sigma_t
    )
    
    pi_mat[idx_p, ] <- MASS::mvrnorm(
      n = sum(idx_p),
      mu = c(5, 5),
      Sigma = Sigma_p
    )
  }
  
  pain <- pmin(pmax(round(pi_mat[,1]), 0), 10)
  itch <- pmin(pmax(round(pi_mat[,2]), 0), 10)
  
  sim$pain <- pain
  sim$itch <- itch
  
  return(sim)
}

# ===========================
# Tests
# ===========================
run_tests <- function(sim, alpha = 0.05, permreps = 10000) {
  
  out <- list()
  # 2) MANOVA.RM
  out$p_manovarm_WTS <- tryCatch({
    resm <- MANOVA.wide(cbind(wound_size_post, pain, itch) ~ treatment_n,
                        data = sim, subject = "patientID")
    as.numeric(resm$WTS[3])
  }, error = function(e) NA)
  
  out$p_manovarm_WTSparamBS <- tryCatch({
    resm <- MANOVA.wide(cbind(wound_size_post, pain, itch) ~ treatment_n,
                        data = sim, resampling= "paramBS", subject = "patientID")
    as.numeric(resm$resampling[,"paramBS (WTS)"])
  }, error = function(e) NA)
  
  out$p_manovarm_WTSwildBS <- tryCatch({
    resm <- MANOVA.wide(cbind(wound_size_post, pain, itch) ~ treatment_n,
                        data = sim, resampling= "WildBS", subject = "patientID")
    as.numeric(resm$resampling[,"WildBS (WTS)"])
  }, error = function(e) NA)
  
  # 2) MANOVA.RM
  out$p_manovarm_MATSparamBS <- tryCatch({
    resm <- MANOVA.wide(cbind(wound_size_post, pain, itch) ~ treatment_n,
                        data = sim, subject = "patientID")
    as.numeric(resm$resampling[,"paramBS (MATS)"])
  }, error = function(e) NA)
  
  out$p_manovarm_MATSwildBS <- tryCatch({
    resm <- MANOVA.wide(cbind(wound_size_post, pain, itch) ~ treatment_n,
                        data = sim, resampling= "WildBS", subject = "patientID")
    as.numeric(resm$resampling[,"WildBS (MATS)"])
  }, error = function(e) NA)
  
  return(out)
  
}

# ===========================
# Job-Plan
# ===========================
jobs <- expand.grid(
  rep = 1:R,
  sample_setting = run_sample_settings,
  scenario = run_scenarios,
  effect = names(effect_settings),
  shift = names(shift_settings),
  stringsAsFactors = FALSE
)

jobs <- jobs %>%
  filter(
    !(scenario == "type1" & shift != names(shift_settings)[1])
  ) %>%
  filter(!(scenario == "power2" ))
cat("\nStarte Simulation mit", n_cores, "Kernen... (R =", R, ")\n")

res_list <- pbmcapply::pbmclapply(seq_len(nrow(jobs)), function(i) {
  row <- jobs[i, ]
  set <- sample_settings[[ row$sample_setting ]]
  eff <- effect_settings[[ row$effect ]]
  shf <- shift_settings[[ row$shift ]]
  
  sim <- simulate_data(
    N = set$N,
    t_size = set$treat,
    p_size = set$placeb,
    scenario = row$scenario,
    treat_resp = eff$treat_resp,
    placeb_resp = eff$placeb_resp,
    mean_shift = ifelse(row$scenario == "power2", shf$mean_shift, 0),
    sd_shift = ifelse(row$scenario == "power2", shf$sd_shift, 0)
  )
  run_tests(sim, alpha = alpha, permreps = permreps)
}, mc.cores = n_cores)

# ===========================
# Ergebnisse umwandeln
# ===========================
results_df <- do.call(rbind, lapply(seq_along(res_list), function(i) {
  pvals <- res_list[[i]]
  row <- jobs[i, ]
  data.frame(
    sample_setting = row$sample_setting,
    scenario = row$scenario,
    effect = row$effect,
    shift = row$shift,
    rep = row$rep,
    
    p_WTS = pvals$p_manovarm_WTS,
    p_WTS_paramBS = pvals$p_manovarm_WTSparamBS,
    p_WTS_wildBS = pvals$p_manovarm_WTSwildBS,
    p_MATS_paramBS = pvals$p_manovarm_MATSparamBS,
    p_MATS_wildBS = pvals$p_manovarm_MATSwildBS
  )
}))


results_df <- results_df %>%
  mutate(
    
    WTS =
      ifelse(!is.na(p_WTS) &
               p_WTS < alpha,1,0),
    
    WTS_paramBS =
      ifelse(!is.na(p_WTS_paramBS) &
               p_WTS_paramBS < alpha,1,0),
    
    WTS_wildBS =
      ifelse(!is.na(p_WTS_wildBS) &
               p_WTS_wildBS < alpha,1,0),
    
    MATS_paramBS =
      ifelse(!is.na(p_MATS_paramBS) &
               p_MATS_paramBS < alpha,1,0),
    
    MATS_wildBS =
      ifelse(!is.na(p_MATS_wildBS) &
               p_MATS_wildBS < alpha,1,0)
  )

df_summary <- results_df %>%
  group_by(sample_setting, scenario, effect, shift) %>%
  summarise(
    across(
      c(
        WTS,
        WTS_paramBS,
        WTS_wildBS,
        MATS_paramBS,
        MATS_wildBS
      ),
      mean
    ),
    
    .groups="drop"
  )

# Ergebnisse als TXT speichern
write.table(df_summary, file = "simulation_results_servEbrevision_allvariantsmanova_covariancehiher_typ1.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)