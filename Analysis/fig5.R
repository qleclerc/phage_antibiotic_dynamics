
library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(RColorBrewer)

source(here::here("Model", "model.R"))

palette = brewer.pal(n = 9, name = "BuGn")[4:9]

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

times = seq(0, 48, 0.1)

yinit = c(Be = 1e9,
          Bt = 1e9,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

all_results = data.frame()


for(abx_con in c(0,0.5,1,2,4,8)){
  
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
    
    #times = seq(0, max(24+results_con$abx_start[i], 24+results_con$phage_start[i]), 1)
    
    results = phage_tr_model(parameters, yinit, times, event_dat)
    
    results_con$max_bet[i] = max(max(results$Bet, na.rm = T), 0.01)
    
    results_con$duration_bet[i] = sum(results$Bet>1)/10
    
    results_con$end_bacteria[i] = max(tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1), 0.01)
    
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
  scale_color_manual(values = palette) +
  # scale_color_manual(values = c("grey", "green", "grey", "grey", "grey", "grey")) +
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
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette) +
  # scale_color_manual(values = c("grey", "green", "grey", "grey", "grey", "grey")) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
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
       y = "Total remaining bacteria after 48h",
       colour = "Antibiotic\nconcentration\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette) +
  # scale_color_manual(values = c("grey", "green", "grey", "grey", "grey", "grey")) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
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

pa = plot_grid(p_remain, p_max, p_duration, ncol = 1,
               labels = c("a)", "c)", "e)"))


all_results = data.frame()


for(pha_con in c(10^6, 10^7, 10^8, 10^9, 10^10, 10^11)){
  
  results_con = data.frame(abx_start = c(rep(1,25), 2:25),
                           phage_start = c(1:25, rep(1,24)),
                           max_bet = 0,
                           duration_bet = 0,
                           end_bacteria = 0)
  
  for(i in 1:nrow(results_con)){
    
    event_dat = data.frame(var = c("ery", "tet", "Pl"),
                           time = c(results_con$abx_start[i], results_con$abx_start[i], results_con$phage_start[i]),
                           value = c(1, 1, pha_con),
                           method = c("add", "add", "add"))
    
   #times = seq(0, max(24+results_con$abx_start[i], 24+results_con$phage_start[i]), 1)
    
    results = phage_tr_model(parameters, yinit, times, event_dat)
    
    results_con$max_bet[i] = max(max(results$Bet, na.rm = T), 0.01)
    
    results_con$duration_bet[i] = sum(results$Bet>1)/10
    
    results_con$end_bacteria[i] = max(tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1), 0.01)
    
  }
  
  results_con$diff = results_con$abx_start - results_con$phage_start
  results_con$pha_con = pha_con
  
  all_results = rbind(all_results, results_con)
  
}


p_duration = ggplot(all_results) +
  geom_line(aes(diff, duration_bet, colour = as.factor(pha_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotic addition time, relative to phage addition time",
       y = "Double-resistant bacteria presence time",
       colour = "Phage\nconcentration\n(pfu/mL):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+06", "1e+07", "1e+08", "1e+09", "1e+10", "1e+11"),
                     labels = c(bquote(10^6),
                                bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10),
                                bquote(10^11))) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_max = ggplot(all_results) +
  geom_line(aes(diff, max_bet, colour = as.factor(pha_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotic addition time, relative to phage addition time",
       y = "Maximum double-resistant bacteria",
       colour = "Phage\nconcentration\n(pfu/mL):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                              expression(10^2), expression(10^4), expression(10^6),
                              expression(10^8), "x")) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+06", "1e+07", "1e+08", "1e+09", "1e+10", "1e+11"),
                     labels = c(bquote(10^6),
                                bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10),
                                bquote(10^11))) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_remain = ggplot(all_results) +
  geom_line(aes(diff, end_bacteria, colour = as.factor(pha_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotic addition time, relative to phage addition time",
       y = "Total remaining bacteria after 48h",
       colour = "Phage\nconcentration\n(pfu/mL):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+06", "1e+07", "1e+08", "1e+09", "1e+10", "1e+11"),
                     labels = c(bquote(10^6),
                                bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10),
                                bquote(10^11))) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
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

pb = plot_grid(p_remain, p_max, p_duration, ncol = 1,
               labels = c("b)", "d)", "f)"))

plot_grid(pa, pb, ncol = 2)

ggsave(here::here("Figures", "fig5.png"), height = 12, width = 11)
