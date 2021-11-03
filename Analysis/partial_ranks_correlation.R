
library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(RColorBrewer)
library(ppcor)
library(dplyr)

source(here::here("Model", "model.R"))

abx_params = read.csv(here::here("Parameters", "abx_params.csv"))
pha_params = read.csv(here::here("Parameters", "pha_params.csv"))
bac_params = read.csv(here::here("Parameters", "bac_params.csv"))

all_results = data.frame(mu_e = rnorm(1000, bac_params$mu_e[1], bac_params$mu_e[1]/100),
                         mu_t = rnorm(1000, bac_params$mu_t[1], bac_params$mu_t[1]/100),
                         mu_et = rnorm(1000, bac_params$mu_et[1], bac_params$mu_et[1]/100),
                         Nmax = rnorm(1000, bac_params$Nmax[1], bac_params$Nmax[1]/100),
                         beta = rnorm(1000, pha_params$beta, pha_params$beta/100),
                         L = rnorm(1000, pha_params$L, pha_params$L/100),
                         tau = rnorm(1000, pha_params$tau, pha_params$tau/100),
                         alpha = rnorm(1000, pha_params$alpha, pha_params$alpha/100),
                         gamma = runif(1000, 0, 0.1),
                         ery_kill_max_BE = abx_params$kmax[1],
                         ery_kill_max_BT = abx_params$kmax[3],
                         tet_kill_max_BE = abx_params$kmax[2],
                         tet_kill_max_BT = abx_params$kmax[4],
                         EC_ery_BE = abx_params$EC50[1],
                         EC_ery_BT = abx_params$EC50[3],
                         EC_tet_BE = abx_params$EC50[2],
                         EC_tet_BT = abx_params$EC50[4],
                         pow_ery_BE = abx_params$pow[1],
                         pow_ery_BT = abx_params$pow[3],
                         pow_tet_BE = abx_params$pow[2],
                         pow_tet_BT = abx_params$pow[4],
                         gamma_ery = runif(1000, 0, 0.1),
                         gamma_tet = runif(1000, 0, 0.1),
                         max_bet = 0,
                         end_bacteria = 0)


times = seq(0, 48, 0.1)

yinit = c(Be = 1e9,
          Bt = 1e9,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)


for(i in 1:nrow(all_results)){
  
  if(i %% round(nrow(all_results)/10) == 0) cat(i/10, "% done\n")
  
  parameters = as.vector(all_results[i,c(1:23)])
  
  event_dat = data.frame(var = c("ery", "tet", "Pl"),
                         time = c(1, 1, 1),
                         value = c(1, 1, 1e9),
                         method = c("add", "add", "add"))
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  all_results$max_bet[i] = max(max(results$Bet, na.rm = T), 0.01)
  all_results$end_bacteria[i] = max(tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1), 0.01)
  
}

all_results = all_results %>%
  select("mu_e", "mu_t", "mu_et", "Nmax", "beta", "L", "tau",
         "alpha", "gamma", "gamma_ery", "gamma_tet", "max_bet", "end_bacteria")

tt = pcor(all_results)

estimates = tt$estimate
colnames(estimates) = colnames(all_results)
rownames(estimates) = colnames(all_results)

pvals = tt$p.value
colnames(pvals) = colnames(all_results)
rownames(pvals) = colnames(all_results)
