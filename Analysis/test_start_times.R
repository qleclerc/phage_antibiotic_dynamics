
library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)

source(here::here("Model", "model.R"))

abx_params = read.csv(here::here("Parameters", "abx_params.csv"))
pha_params = read.csv(here::here("Parameters", "pha_params.csv"))
bac_params = read.csv(here::here("Parameters", "bac_params.csv"))

parameters = c(mu_e = bac_params$mu_e[1],
               mu_t = bac_params$mu_t[1],
               mu_et = bac_params$mu_et[1],
               Nmax = bac_params$Nmax[1],
               beta = pha_params$beta,
               L = pha_params$L,
               tau = pha_params$tau,
               alpha = pha_params$alpha,
               gamma = 0,
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
               gamma_ery = 0,
               gamma_tet = 0)

times = seq(0, 24, 1)

yinit = c(Be = 1e9,
          Bt = 1e9,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

all_results = data.frame()


for(abx_con in c(0.5,1,1.5,2,2.5,3)){
  
  results_con = data.frame(abx_start = c(rep(1,25), 2:25),
                           phage_start = c(1:25, rep(1,24)),
                           max_bet = 0,
                           duration_bet = 0,
                           end_bacteria = 0)
  
  for(i in 1:nrow(results_con)){
    
    event_dat = data.frame(var = c("ery", "tet", "Pl"),
                           time = c(results_con$abx_start[i], results_con$abx_start[i], results_con$phage_start[i]),
                           value = c(abx_con, abx_con, 1e9),
                           method = c("add", "add", "add"))
    
    times = seq(0, max(24+results_con$abx_start[i], 24+results_con$phage_start[i]), 1)
    
    results = phage_tr_model(parameters, yinit, times, event_dat)
    
    results_con$max_bet[i] = max(results$Bet, na.rm = T)
    
    results_con$duration_bet[i] = length(which(results$Bet>1))
    
    results_con$end_bacteria[i] = tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1)
    
  }
  
  results_con$diff = results_con$abx_start - results_con$phage_start
  results_con$abx_con = abx_con
  
  all_results = rbind(all_results, results_con)
  
}


p_duration = ggplot(all_results) +
  geom_line(aes(diff, duration_bet, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotic addition time, relative to phage addition time",
       y = "Double-resistant bacteria presence time",
       colour = "Antibiotic\nconcentration\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_max = ggplot(all_results) +
  geom_line(aes(diff, max_bet, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotic addition time, relative to phage addition time",
       y = "Maximum double-resistant bacteria",
       colour = "Antibiotic\nconcentration\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 5e6)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_remain = ggplot(all_results) +
  geom_line(aes(diff, end_bacteria, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotic addition time, relative to phage addition time",
       y = "Total remaining bacteria after treatment",
       colour = "Antibiotic\nconcentration\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 5e6)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

# ggplot(all_results) +
#   geom_tile(aes(phage_start, abx_start, fill = duration_bet)) +
#   theme_bw()
# 
# ggplot(all_results) +
#   geom_tile(aes(phage_start, abx_start, fill = log10(max_bet))) +
#   theme_bw()

plot_grid(p_duration, p_max, p_remain)

ggsave(here::here("Figures", "timings.png"), height = 8, width = 10)
