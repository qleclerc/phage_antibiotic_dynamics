
library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(RColorBrewer)
library(epiR)

source(here::here("Model", "model.R"))
set.seed(42)

palette = c(brewer.pal(n = 9, name = "BuGn")[c(3,4,6,8)], "#000000")

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
               ery_kill_max_BET = abx_params$kmax[5],
               tet_kill_max_BE = abx_params$kmax[2],
               tet_kill_max_BT = abx_params$kmax[4],
               tet_kill_max_BET = abx_params$kmax[6],
               EC_ery_BE = abx_params$EC50[1],
               EC_ery_BT = abx_params$EC50[3],
               EC_ery_BET = abx_params$EC50[5],
               EC_tet_BE = abx_params$EC50[2],
               EC_tet_BT = abx_params$EC50[4],
               EC_tet_BET = abx_params$EC50[6],
               pow_ery_BE = abx_params$pow[1],
               pow_ery_BT = abx_params$pow[3],
               pow_ery_BET = abx_params$pow[5],
               pow_tet_BE = abx_params$pow[2],
               pow_tet_BT = abx_params$pow[4],
               pow_tet_BET = abx_params$pow[6],
               gamma_ery = 0,
               gamma_tet = 0)

times = seq(0, 24, 0.1)

yinit = c(Be = 1e9,
          Bt = 1e9,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)


event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 0, 0),
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

all_results = data.frame()

for(tr_param in c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10)){
  
  parameters[["alpha"]] = tr_param
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  results  = results %>%
    select(time, Bet, Be, Bt) %>%
    mutate(tr_param = tr_param)
  
  all_results = rbind(all_results, results)
  
}

pa = ggplot(all_results) +
  geom_line(aes(time, Bet, colour = as.factor(tr_param), linetype = "Double-resistant"), size = 0.8) +
  geom_line(aes(time, Be+Bt, colour = as.factor(tr_param), linetype = "Single-resistant"),
            alpha = 0.5,size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("segment", x = 0+4, xend = 0, y = 1e4, yend = 1e4) +
  geom_label(x = 0+4, y = 4, label = "Phage +", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("segment", x = 0+5, xend = 0, y = 1e9, yend = 1e9) +
  geom_label(x = 0+5, y = 9, label = "Antibiotics +", size = 3) +
  theme_bw() +
  labs(x = "Time (hours)",
       y = "Double-resistant bacteria (cfu/mL)",
       colour = "Transducing phage probability:",
       linetype = "Bacteria:") +
  scale_x_continuous(breaks = seq(0,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_linetype_manual(breaks = c("Single-resistant", "Double-resistant"),
                        values = c("dashed", "solid")) +
  scale_color_manual(values = palette, 
                     breaks = c("1e-10", "1e-09", "1e-08",
                                "1e-07", "1e-06"),
                     labels = c(bquote("1 \u00D7 " * 10^-10),
                                bquote("1 \u00D7 " * 10^-9),
                                bquote("1 \u00D7 " * 10^-8),
                                bquote("1 \u00D7 " * 10^-7),
                                bquote("1 \u00D7 " * 10^-6))) +
  coord_cartesian(ylim = c(0.1, 5e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(10, 10, 0),
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

all_results = data.frame()

for(tr_param in c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10)){
  
  parameters[["alpha"]] = tr_param
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  results  = results %>%
    select(time, Bet, Be, Bt) %>%
    mutate(tr_param = tr_param)
  
  all_results = rbind(all_results, results)
  
}

pb = ggplot(all_results) +
  geom_line(aes(time, Bet, colour = as.factor(tr_param), linetype = "Double-resistant"), size = 0.8) +
  geom_line(aes(time, Be+Bt, colour = as.factor(tr_param), linetype = "Single-resistant"),
            alpha = 0.5,size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("segment", x = 0+4, xend = 0, y = 1e4, yend = 1e4) +
  geom_label(x = 0+4, y = 4, label = "Phage +", size = 3) +
  geom_vline(xintercept = 10, linetype = "dashed") +
  annotate("segment", x = 10+5, xend = 10, y = 1e9, yend = 1e9) +
  geom_label(x = 10+5, y = 9, label = "Antibiotics +", size = 3) +
  theme_bw() +
  labs(x = "Time (hours)",
       y = "Double-resistant bacteria (cfu/mL)",
       colour = "Transducing phage probability:",
       linetype = "Bacteria:") +
  scale_x_continuous(breaks = seq(0,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_linetype_manual(breaks = c("Single-resistant", "Double-resistant"),
                        values = c("dashed", "solid")) +
  scale_color_manual(values = palette, 
                     breaks = c("1e-10", "1e-09", "1e-08",
                                "1e-07", "1e-06"),
                     labels = c(bquote("1 \u00D7 " * 10^-10),
                                bquote("1 \u00D7 " * 10^-9),
                                bquote("1 \u00D7 " * 10^-8),
                                bquote("1 \u00D7 " * 10^-7),
                                bquote("1 \u00D7 " * 10^-6))) +
  coord_cartesian(ylim = c(0.1, 5e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 0, 10),
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

all_results = data.frame()

for(tr_param in c(1e-6, 1e-7, 1e-8, 1e-9, 1e-10)){
  
  parameters[["alpha"]] = tr_param
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  results  = results %>%
    select(time, Bet, Be, Bt) %>%
    mutate(tr_param = tr_param)
  
  all_results = rbind(all_results, results)
  
}

pc = ggplot(all_results) +
  geom_line(aes(time, Bet, colour = as.factor(tr_param), linetype = "Double-resistant"), size = 0.8) +
  geom_line(aes(time, Be+Bt, colour = as.factor(tr_param), linetype = "Single-resistant"),
            alpha = 0.5,size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 10, linetype = "dashed") +
  annotate("segment", x = 10+4, xend = 10, y = 1e4, yend = 1e4) +
  geom_label(x = 10+4, y = 4, label = "Phage +", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("segment", x = 0+5, xend = 0, y = 1e9, yend = 1e9) +
  geom_label(x = 0+5, y = 9, label = "Antibiotics +", size = 3) +
  theme_bw() +
  labs(x = "Time (hours)",
       y = "Double-resistant bacteria (cfu/mL)",
       colour = "Transducing phage probability:",
       linetype = "Bacteria:") +
  scale_x_continuous(breaks = seq(0,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_linetype_manual(breaks = c("Single-resistant", "Double-resistant"),
                        values = c("dashed", "solid")) +
  scale_color_manual(values = palette, 
                     breaks = c("1e-10", "1e-09", "1e-08",
                                "1e-07", "1e-06"),
                     labels = c(bquote("1 \u00D7 " * 10^-10),
                                bquote("1 \u00D7 " * 10^-9),
                                bquote("1 \u00D7 " * 10^-8),
                                bquote("1 \u00D7 " * 10^-7),
                                bquote("1 \u00D7 " * 10^-6))) +
  coord_cartesian(ylim = c(0.1, 5e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  



## sensitivity ###############

all_results = data.frame(mu_e = runif(500, bac_params$mu_e[2], bac_params$mu_e[3]),
                         mu_t = runif(500, bac_params$mu_t[2], bac_params$mu_t[3]),
                         mu_et = runif(500, bac_params$mu_et[2], bac_params$mu_et[3]),
                         Nmax = bac_params$Nmax[1],
                         beta = runif(500, pha_params$beta_0.975, pha_params$beta_0.025),
                         L = runif(500, pha_params$L_0.025, pha_params$L_0.975),
                         tau = runif(500, pha_params$tau_0.025, pha_params$tau_0.975),
                         alpha = runif(500, pha_params$alpha_0.975, pha_params$alpha_0.025),
                         gamma = runif(500, 0, 0.1),
                         ery_kill_max_BE = abx_params$kmax[1],
                         ery_kill_max_BT = abx_params$kmax[3],
                         ery_kill_max_BET = abx_params$kmax[5],
                         tet_kill_max_BE = abx_params$kmax[2],
                         tet_kill_max_BT = abx_params$kmax[4],
                         tet_kill_max_BET = abx_params$kmax[6],
                         EC_ery_BE = abx_params$EC50[1],
                         EC_ery_BT = abx_params$EC50[3],
                         EC_ery_BET = abx_params$EC50[5],
                         EC_tet_BE = abx_params$EC50[2],
                         EC_tet_BT = abx_params$EC50[4],
                         EC_tet_BET = abx_params$EC50[6],
                         pow_ery_BE = abx_params$pow[1],
                         pow_ery_BT = abx_params$pow[3],
                         pow_ery_BET = abx_params$pow[5],
                         pow_tet_BE = abx_params$pow[2],
                         pow_tet_BT = abx_params$pow[4],
                         pow_tet_BET = abx_params$pow[6],
                         gamma_ery = runif(500, 0, 0.1),
                         gamma_tet = runif(500, 0, 0.1),
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
  
  if(i %% round(nrow(all_results)/10) == 0) cat(i/round(nrow(all_results))*100, "% done\n")
  
  parameters = as.vector(all_results[i,c(1:29)])
  
  event_dat = data.frame(var = c("ery", "tet", "Pl"),
                         time = c(0, 0, 0),
                         value = c(1, 1, 1e9),
                         method = c("add", "add", "add"))
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  all_results$max_bet[i] = max(max(results$Bet, na.rm = T), 0)
  all_results$end_bacteria[i] = max(tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1), 0)
  
}

all_results_m = all_results
all_results = all_results_m %>%
  select(mu_e, mu_t, mu_et, beta, L, tau,
         alpha, gamma, gamma_ery, gamma_tet, end_bacteria)

cor_end_bac = epi.prcc(all_results)
cor_end_bac$param = colnames(all_results)[-ncol(all_results)]
cor_end_bac$cor = "end_bac"

all_results = all_results_m %>%
  select(mu_e, mu_t, mu_et, beta, L, tau,
         alpha, gamma, gamma_ery, gamma_tet, max_bet)

cor_max_bet = epi.prcc(all_results)
cor_max_bet$param = colnames(all_results)[-ncol(all_results)]
cor_max_bet$cor = "max_bet"

all_results = rbind(cor_end_bac, cor_max_bet)
all_results$param = as.factor(all_results$param)
all_results$param = factor(all_results$param,
                           levels = levels(all_results$param)[c(2,6,10,1,3:5,7,9,8)])

pd = ggplot(all_results) +
  geom_pointrange(aes(x = param, y = est, group = cor, 
                      ymin = lower, ymax = upper, colour = cor),
                  position=position_dodge(width=0.2), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12)) +
  labs(colour = "", x = "\n", y = "Correlation coefficient") +
  scale_color_discrete(labels = c("Remaining bacteria", bquote("Maximum B"[ET]))) +
  scale_x_discrete(labels = c(bquote(beta),
                              bquote(delta*"max"),
                              bquote(tau),
                              bquote(alpha),
                              bquote(gamma[P]),
                              bquote(gamma[E]),
                              bquote(gamma[T]),
                              bquote(mu*"max"[E]),
                              bquote(mu*"max"[T]), 
                              bquote(mu*"max"[ET]))) +
  coord_cartesian(clip = "off", ylim = c(-1,1)) +
  annotate("text", x = 2.5, y = -1.37, label = "Phage parameters") +
  annotate("text", x = 6, y = -1.37, label = "Decay parameters") +
  annotate("text", x = 9, y = -1.37, label = "Bacteria parameters") +
  annotate("segment", x = 1, xend = 4, y = -1.30, yend = -1.30) +
  annotate("segment", x = 5, xend = 7, y = -1.30, yend = -1.30) +
  annotate("segment", x = 8, xend = 10, y = -1.30, yend = -1.30)


## final plot ############

plot_grid(plot_grid(plot_grid(pa + theme(legend.position = "none"),
                              pb + theme(legend.position = "none"),
                              pc + theme(legend.position = "none"),
                              ncol = 3,
                              labels = c("a)", "b)", "c)")),
                    get_legend(pa + guides(linetype = F) + theme(legend.position = "bottom")),
                    get_legend(pa + guides(colour = F) + theme(legend.position = "bottom",
                                                               legend.key.width = unit(2, "cm"))),
                    nrow = 3,
                    rel_heights = c(1,0.1,0.1)),
          pd, rel_heights = c(1,0.7), nrow = 2, labels = c("", "d)"))

ggsave(here::here("Figures", "fig6.png"), height = 10, width = 10)
