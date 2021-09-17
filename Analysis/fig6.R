
library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(RColorBrewer)

source(here::here("Model", "model.R"))

palette = c(brewer.pal(n = 9, name = "BuGn")[-1], "#000000")

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


event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 0, 0),
                       value = c(0, 0, 1e9),
                       method = c("add", "add", "add"))

all_results = data.frame()

for(tr_param in c(1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 5e-10, 1e-10)){
  
  parameters[["alpha"]] = tr_param
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  results  = results %>%
    select(time, Bet, Be, Bt) %>%
    mutate(tr_param = tr_param)
  
  all_results = rbind(all_results, results)
  
  # max_bet = max(max(results$Bet, na.rm = T), 0.01)
  # 
  # duration_bet = sum(results$Bet>1)
  # 
  # end_bacteria = max(tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1), 0.01)
  # 
  # all_results = rbind(all_results,
  #                     data.frame(tr = tr_param, max_bet = max_bet,
  #                                duration_bet = duration_bet, end_bacteria = end_bacteria))
  
}

pa = ggplot(all_results) +
  geom_line(aes(time, Bet, colour = as.factor(tr_param), linetype = "Double-resistant"), size = 0.8) +
  geom_line(aes(time, Be+Bt, colour = as.factor(tr_param), linetype = "Single-resistant"),
            alpha = 0.1,size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw() +
  labs(x = "Time (hours)",
       y = "Double-resistant bacteria (cfu/mL)",
       colour = "Transducing phage probability:",
       linetype = "Bacteria:") +
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_linetype_manual(breaks = c("Single-resistant", "Double-resistant"),
                        values = c("dashed", "solid")) +
  scale_color_manual(values = palette, 
                     breaks = c("1e-10", "5e-10", "1e-09", "5e-09", "1e-08",
                                "5e-08", "1e-07", "5e-07", "1e-06"),
                     labels = c(bquote("1 \u00D7 " * 10^-10),
                                bquote("5 \u00D7 " * 10^-10),
                                bquote("1 \u00D7 " * 10^-9),
                                bquote("5 \u00D7 " * 10^-9),
                                bquote("1 \u00D7 " * 10^-8),
                                bquote("5 \u00D7 " * 10^-8),
                                bquote("1 \u00D7 " * 10^-7),
                                bquote("5 \u00D7 " * 10^-7),
                                bquote("1 \u00D7 " * 10^-6))) +
  coord_cartesian(ylim = c(0.1, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 0, 0),
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

all_results = data.frame()

for(tr_param in c(1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9, 5e-10, 1e-10)){
  
  parameters[["alpha"]] = tr_param
  
  results = phage_tr_model(parameters, yinit, times, event_dat)
  
  results  = results %>%
    select(time, Bet, Be, Bt) %>%
    mutate(tr_param = tr_param)
  
  all_results = rbind(all_results, results)
  
  # max_bet = max(max(results$Bet, na.rm = T), 0.01)
  # 
  # duration_bet = sum(results$Bet>1)
  # 
  # end_bacteria = max(tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1), 0.01)
  # 
  # all_results = rbind(all_results,
  #                     data.frame(tr = tr_param, max_bet = max_bet,
  #                                duration_bet = duration_bet, end_bacteria = end_bacteria))
  
}

pb = ggplot(all_results) +
  geom_line(aes(time, Bet, colour = as.factor(tr_param), linetype = "Double-resistant"), size = 0.8) +
  geom_line(aes(time, Be+Bt, colour = as.factor(tr_param), linetype = "Single-resistant"),
            alpha = 0.1,size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme_bw() +
  labs(x = "Time (hours)",
       y = "Double-resistant bacteria (cfu/mL)",
       colour = "Transducing phage probability:",
       linetype = "Bacteria:") +
  scale_x_continuous(breaks = seq(0,24,2)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_linetype_manual(breaks = c("Single-resistant", "Double-resistant"),
                        values = c("dashed", "solid")) +
  scale_color_manual(values = palette, 
                     breaks = c("1e-10", "5e-10", "1e-09", "5e-09", "1e-08",
                                "5e-08", "1e-07", "5e-07", "1e-06"),
                     labels = c(bquote("1 \u00D7 " * 10^-10),
                                bquote("5 \u00D7 " * 10^-10),
                                bquote("1 \u00D7 " * 10^-9),
                                bquote("5 \u00D7 " * 10^-9),
                                bquote("1 \u00D7 " * 10^-8),
                                bquote("5 \u00D7 " * 10^-8),
                                bquote("1 \u00D7 " * 10^-7),
                                bquote("5 \u00D7 " * 10^-7),
                                bquote("1 \u00D7 " * 10^-6))) +
  coord_cartesian(ylim = c(0.1, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

plot_grid(plot_grid(pa + theme(legend.position = "none"),
                    pb + theme(legend.position = "none"),
                    ncol = 2,
                    labels = c("a)", "b)")),
          get_legend(pa + guides(linetype = F) + theme(legend.position = "bottom")),
          get_legend(pa + guides(colour = F) + theme(legend.position = "bottom")),
          nrow = 3,
          rel_heights = c(1,0.1,0.1))

ggsave(here::here("Figures", "fig6.png"), height = 6, width = 10)
