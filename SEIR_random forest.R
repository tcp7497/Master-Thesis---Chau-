# ============================================================
# FULL RANDOM FOREST SENSITIVITY ANALYSIS (CHATGPT)
# ============================================================

rm(list = ls())

# ---------------------------
# 1. Packages
# ---------------------------
packs <- c("randomForest", "dplyr", "ggplot2")
for (p in packs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(randomForest)
library(dplyr)
library(ggplot2)

set.seed(20260401)

# ---------------------------
# 2. Load simulation functions
# ---------------------------
source("C:\\Users\\Chau\\Downloads\\F\\Statistics UHasselt\\Semester 4\\Master Thesis\\Master-Thesis---Chau-\\SEIR_contact-rates.R")  

# ---------------------------
# 3. Fixed settings
# ---------------------------
n <- 28
t.stop <- 360
roles <- c(rep("vol", 24), rep("staff", 4))
eps <- 1e-6

logit_clip <- function(x) {
  x <- pmin(pmax(x, eps), 1 - eps)
  log(x / (1 - x))
}

# ---------------------------
# 4. Scenario setup
# ---------------------------
get_scenario <- function(name) {
  
  if (name == "free") {
    list(group_assign = NULL, M = NULL)
  }
  
  else if (name == "staff_bridge") {
    group_assign <- c(rep(1,12), rep(2,12), rep(3,4))
    M <- matrix(c(
      0.90,0.00,0.10,
      0.00,0.90,0.10,
      0.40,0.40,0.20
    ),3,byrow=TRUE)
    list(group_assign = group_assign, M = M)
  }
  
  else {
    group_assign <- c(rep(1,12), rep(2,12), rep(3,2), rep(4,2))
    M <- matrix(c(
      0.90,0.00,0.10,0.00,
      0.00,0.90,0.00,0.10,
      0.50,0.00,0.50,0.00,
      0.00,0.50,0.00,0.50
    ),4,byrow=TRUE)
    list(group_assign = group_assign, M = M)
  }
}

# ---------------------------
# 5. Extract outcomes
# ---------------------------
get_AR <- function(sim) {
  
  sm <- sim$status.matrix
  infected <- !is.na(sm$time.of.infection)
  
  vol <- sm$roles == "vol"
  staff <- sm$roles == "staff"
  
  AR_total <- mean(infected)
  AR_vol <- mean(infected[vol])
  AR_staff <- mean(infected[staff])
  
  data.frame(
    AR_total = AR_total,
    AR_vol = AR_vol,
    AR_staff = AR_staff,
    logit_AR_total = logit_clip(AR_total),
    logit_AR_vol = logit_clip(AR_vol),
    logit_AR_staff = logit_clip(AR_staff)
  )
}

# ---------------------------
# 6. Parameter design
# ---------------------------
n_design <- 200
n_rep <- 20

design <- data.frame(
  scenario = sample(c("free","staff_bridge","strict"), n_design, TRUE),
  R0 = runif(n_design, 1.5, 4.0),
  ct.rates = runif(n_design, 1.5, 5),
  nSeeds = sample(1, n_design, TRUE),
  k = runif(n_design, 0.1, 1.5),
  k_ct = runif(n_design, 0.5, 2),
  k_sus = runif(n_design, 0.5, 2),
  mean_e = runif(n_design, 4, 7),
  sd_e = runif(n_design, 1.5, 3.5),
  mean_i = runif(n_design, 7, 12),
  sd_i = runif(n_design, 2.5, 5)
)

# constraint
design$sd_e <- pmin(design$sd_e, design$mean_e - 0.2)
design$sd_i <- pmin(design$sd_i, design$mean_i - 0.2)

# ---------------------------
# 7. Run simulations
# ---------------------------
results <- list()

for (i in 1:n_design) {
  
  if (i %% 20 == 0) cat("Running:", i, "/", n_design, "\n")
  
  pars <- design[i, ]
  sc <- get_scenario(pars$scenario)
  
  rep_out <- list()
  
  for (r in 1:n_rep) {
    
    sim <- sim.epi_seir(
      n = n,
      ct.rates = pars$ct.rates,
      nSeeds = pars$nSeeds,
      R0 = pars$R0,
      t.stop = t.stop,
      roles = roles,
      average_mod = average_mod,
      k = pars$k,
      k_ct = pars$k_ct,
      k_sus = pars$k_sus,
      mean_e = pars$mean_e,
      sd_e = pars$sd_e,
      mean_i = pars$mean_i,
      sd_i = pars$sd_i,
      group_assign = sc$group_assign,
      M = sc$M
    )
    
    rep_out[[r]] <- get_AR(sim)
  }
  
  df <- bind_rows(rep_out)
  
  results[[i]] <- cbind(
    pars,
    summarise(df,
              logit_AR_total = mean(logit_AR_total),
              logit_AR_vol = mean(logit_AR_vol),
              logit_AR_staff = mean(logit_AR_staff)
    )
  )
}

rf_data <- bind_rows(results)
rf_data$scenario <- as.factor(rf_data$scenario)

# ---------------------------
# 8. Random Forest
# ---------------------------
rf_total <- randomForest(
  logit_AR_total ~ scenario + R0 + ct.rates + nSeeds + 
    k + k_ct + k_sus + mean_e + sd_e + mean_i + sd_i, data = rf_data,
  ntree = 500, importance = TRUE, proximity = TRUE
)

rf_vol <- randomForest(
  logit_AR_vol ~ scenario + R0 + ct.rates +
    k + k_ct + k_sus + mean_e + sd_e + mean_i + sd_i, data = rf_data,
  ntree = 500, importance = TRUE
)

rf_staff <- randomForest(
  logit_AR_staff ~ scenario + R0 + ct.rates +
    k + k_ct + k_sus + mean_e + sd_e + mean_i + sd_i, data = rf_data,
  ntree = 500, importance = TRUE
)

# ---------------------------
# 9. Variable Importance
# ---------------------------
get_imp <- function(rf, name) {
  imp <- importance(rf)[, "%IncMSE"]
  data.frame(
    parameter = names(imp),
    importance = imp,
    outcome = name
  )
}

imp_df <- bind_rows(
  get_imp(rf_total, "total"),
  get_imp(rf_vol, "volunteers"),
  get_imp(rf_staff, "staff")
)

# Plot
ggplot(imp_df, aes(x=reorder(parameter, importance), y=importance)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~outcome, scales="free") +
  theme_minimal() +
  ggtitle("Variable Importance")

# ---------------------------
# 10. Proximity plot
# ---------------------------
prox <- rf_total$proximity
mds <- cmdscale(1 - prox)

plot_df <- data.frame(
  x = mds[,1],
  y = mds[,2],
  scenario = rf_data$scenario
)

ggplot(plot_df, aes(x=x, y=y, color=scenario)) +
  geom_point(size=2) +
  theme_minimal() +
  ggtitle("Proximity Plot")