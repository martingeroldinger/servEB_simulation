# if(!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(npmv, MANOVA.RM, tidyverse, MASS, pbmcapply, rankFD, dplyr)

setwd("/home/wlauth/Documents/Martin/")
set.seed(2025)

# ===========================
# Parameter & Optionen
# ===========================
R <- 1000                   # Anzahl Simulationen pro Setting
alpha <- 0.05
permreps <- 1000             # Permutationswiederholungen für kleine N
n_cores <- parallel::detectCores() - 1

# Szenarienauswahl
run_scenarios <- c("type1", "power1", "power2")  

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
  eff80_20 = list(treat_resp = 0.8, placeb_resp = 0.2),
  eff60_40 = list(treat_resp = 0.6, placeb_resp = 0.4)
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
 # std <- 60
  sim$wound_size_pre <- round(rlnorm(N, meanlog=4, sdlog=1.2))
  
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
  Sigma <- matrix(c(9, 3.6, 3.6, 9), nrow = 2)
  pi_mat <- MASS::mvrnorm(n = N, mu = c(5, 5), Sigma = Sigma)
  pain <- pmin(pmax(round(pi_mat[,1]), 0), 10)
  itch <- pmin(pmax(round(pi_mat[,2]), 0), 10)
  
  # power2: Shift bei Pain & Itch
  if (scenario == "power2") {
    idx_pl <- which(sim$treatment_n == 0)
    if (length(idx_pl) > 0) {
      z_vals <- rnorm(length(idx_pl), mean = mean_shift, sd = sd_shift)
      z_vals[z_vals < 0] <- 0
      pain[idx_pl] <- pmin(pmax(pain[idx_pl] + z_vals, 0), 10)
      itch[idx_pl] <- pmin(pmax(itch[idx_pl] + z_vals, 0), 10)
    }
  }
  
  sim$pain <- pain
  sim$itch <- itch
  
  return(sim)
}

# ===========================
# Tests
# ===========================
run_tests <- function(sim, alpha = 0.05, permreps = 1000) {
  out <- list()
  
  # 1) npmv
  out$p_npmv <- tryCatch({
    res <- nonpartest(wound_size_post | pain | itch ~ treatment, sim,
                      test = c(1,0,0,0),
                      permtest = (nrow(sim) < 30),
                      permreps = ifelse(nrow(sim) < 30, permreps, 0),
                      plots = FALSE)
    res$results$`P-value`[1]
  }, error = function(e) NA)
  
  # 2) MANOVA.RM
  out$p_manovarm <- tryCatch({
    resm <- MANOVA.wide(cbind(wound_size_post, pain, itch) ~ treatment_n,
                        data = sim, subject = "patientID")
    as.numeric(resm$WTS[3])
  }, error = function(e) NA)
  
  # 3) Mann-Whitney + Holm
  out$mw_reject <- tryCatch({
    p1 <- wilcox.test(wound_size_post ~ treatment, data = sim)$p.value
    p2 <- wilcox.test(as.numeric(pain) ~ treatment, data = sim)$p.value
    p3 <- wilcox.test(as.numeric(itch) ~ treatment, data = sim)$p.value
    padj <- p.adjust(c(p1,p2,p3), method = "holm")
    any(padj < alpha)
  }, error = function(e) NA)
  
  # 4) rankFD::rank.two.samples (Mann-Whitney)
  out$mw_rankfd <- tryCatch({
    p1 <- rankFD::rank.two.samples(wound_size_post ~ treatment, data = sim)$ANOVA.test[1, "p.value"]
    p2 <- rankFD::rank.two.samples(pain ~ treatment, data = sim)$ANOVA.test[1, "p.value"]
    p3 <- rankFD::rank.two.samples(itch ~ treatment, data = sim)$ANOVA.test[1, "p.value"]
    padj <- p.adjust(c(p1, p2, p3), method = "holm")
    any(padj < alpha)
  }, error = function(e) NA)
  
  
  # 5) klassische MANOVA
  out$p_manova_param <- tryCatch({
    fit <- manova(cbind(wound_size_post, pain, itch) ~ treatment, data = sim)
    s <- summary(fit, test = "Pillai")
    s$stats[1, "Pr(>F)"]
  }, error = function(e) NA)
  
  attr(out, "data") <- sim
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

# Unnötige Kombinationen filtern: Shift nur für power2
jobs <- jobs %>%
  filter(!(scenario != "power2" & shift != names(shift_settings)[1])) %>%
  filter(!(scenario == "type1" & effect != names(effect_settings)[1])) # z.B. nur eff80_20 bei type1

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
    p_npmv = pvals$p_npmv,
    p_manovarm = pvals$p_manovarm,
    p_manova_param = pvals$p_manova_param,
    mw_reject = ifelse(is.na(pvals$mw_reject), FALSE, pvals$mw_reject),
    mw_rankfd_reject = ifelse(is.na(pvals$mw_rankfd), FALSE, pvals$mw_rankfd),
    stringsAsFactors = FALSE
  )
}))


results_df <- results_df %>%
  mutate(
    npmv = ifelse(!is.na(p_npmv) & p_npmv < alpha, 1, 0),
    manovarm = ifelse(!is.na(p_manovarm) & p_manovarm < alpha, 1, 0),
    manova = ifelse(!is.na(p_manova_param) & p_manova_param < alpha, 1, 0),
    mw = ifelse(mw_reject, 1, 0),
    mw_rankfd = ifelse(mw_rankfd_reject, 1, 0))


df_summary <- results_df %>%
  group_by(sample_setting, scenario, effect, shift) %>%
  summarise(across(c(npmv, manovarm, manova, mw, mw_rankfd), mean), .groups = "drop")

# ===========================
# Mittelwerte pro Gruppe
# ===========================
if (calc_means) {
  means_all <- do.call(rbind, lapply(res_list, function(x) {
    sim_data <- attr(x, "data")
    if (is.null(sim_data)) return(NULL)
    sim_data %>%
      group_by(treatment) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
      pivot_longer(-treatment, names_to = "variable", values_to = "mean") %>%
      pivot_wider(names_from = treatment, values_from = mean, names_prefix = "mean_")
  }))
  mean_mat <- means_all %>%
    group_by(variable) %>%
    summarise(across(starts_with("mean_"), mean, na.rm = TRUE), .groups = "drop")
  cat("\n=== Gruppen-Mittelwerte (über alle Replikationen) ===\n")
  print(mean_mat)
}

# ===========================
# Ausgabe
# ===========================
cat("\n=== Zusammenfassung (Anteil signifikanter Tests) ===\n")
print(df_summary)

if (show_data_example) {
  cat("\n=== Beispiel-Datensatz ===\n")
  print(head(attr(res_list[[1]], "data"), 10))
}

write.table(df_summary, paste0("df_summary_rlnorm",format(Sys.Date(), "%d%m%Y"),".txt"))
