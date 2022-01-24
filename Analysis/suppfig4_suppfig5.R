

library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(RColorBrewer)
library(dplyr)

source(here::here("Model", "model.R"))

palette = rev(brewer.pal(n = 9, name = "RdBu"))[c(1,3,7,9)]

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

times = seq(0, 48, 0.1)

yinit = c(Be = 1e9,
          Bt = 1e6,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

all_results = data.frame()

for(pha_con in c(10^7, 10^8, 10^9, 10^10)){
  
  results_con = data.frame(abx_start = c(rep(0,25), 1:24),
                           phage_start = c(0:24, rep(0,24)),
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

#things to test
test_pha_con = 1e8
test_diffs = c(0,3,5,15)

p_duration = ggplot(all_results) +
  geom_line(aes(diff, duration_bet, colour = as.factor(pha_con)), size = 0.8) +
  geom_point(data = all_results %>% filter(pha_con == test_pha_con) %>% filter(diff %in% test_diffs),
             aes(diff, duration_bet, shape = as.factor(diff)), size = 3, fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Double-resistant bacteria presence time",
       colour = "Phage added\n(pfu/mL):",
       shape = "Condition:") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(limits = c(0, 47)) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+07", "1e+08", "1e+09", "1e+10"),
                     labels = c(bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10))) +
  scale_shape_manual(values = c(21,22,8,24),
                     breaks = as.factor(test_diffs),
                     labels = c("1.", "2.", "3.", "4.")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_max = ggplot(all_results) +
  geom_line(aes(diff, max_bet, colour = as.factor(pha_con)), size = 0.8) +
  geom_point(data = all_results %>% filter(pha_con == test_pha_con) %>% filter(diff %in% test_diffs),
             aes(diff, max_bet, shape = as.factor(diff)), size = 3, fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Maximum double-resistant bacteria",
       colour = "Phage added\n(pfu/mL):",
       shape = "Condition:") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+07", "1e+08", "1e+09", "1e+10"),
                     labels = c(bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10))) +
  scale_shape_manual(values = c(21,22,8,24),
                     breaks = as.factor(test_diffs),
                     labels = c("1.", "2.", "3.", "4.")) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_remain = ggplot(all_results) +
  geom_line(aes(diff, end_bacteria, colour = as.factor(pha_con)), size = 0.8) +
  geom_point(data = all_results %>% filter(pha_con == test_pha_con) %>% filter(diff %in% test_diffs),
             aes(diff, end_bacteria, shape = as.factor(diff)), size = 3, fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Total remaining bacteria after 48h",
       colour = "Phage added\n(pfu/mL):",
       shape = "Condition:") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+07", "1e+08", "1e+09", "1e+10"),
                     labels = c(bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10))) +
  scale_shape_manual(values = c(21,22,8,24),
                     breaks = as.factor(test_diffs),
                     labels = c("1.", "2.", "3.", "4.")) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


pb = plot_grid(p_remain + theme(legend.position = "none"),
               p_max + theme(legend.position = "none"),
               p_duration + theme(legend.position = "none"),
               get_legend(p_remain + guides(shape = F) + theme(legend.position = "bottom")),
               ncol = 1,
               rel_heights = c(1,1,1,0.15,0.15))


all_results_abx = data.frame()

for(abx_con in c(0.25, 0.5, 1, 2)){
  
  results_con = data.frame(abx_start = c(rep(0,25), 1:24),
                           phage_start = c(0:24, rep(0,24)),
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
  
  all_results_abx = rbind(all_results_abx, results_con)
  
}


p_duration = ggplot(all_results_abx) +
  geom_line(aes(diff, duration_bet, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Double-resistant bacteria presence time",
       colour = "Antibiotics added\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(limits = c(0, 47)) +
  scale_color_manual(values = palette) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_max = ggplot(all_results_abx) +
  geom_line(aes(diff, max_bet, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Maximum double-resistant bacteria",
       colour = "Antibiotics added\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_remain = ggplot(all_results_abx) +
  geom_line(aes(diff, end_bacteria, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Total remaining bacteria after 48h",
       colour = "Antibiotics added\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


pa = plot_grid(p_remain + theme(legend.position = "none"),
               p_max + theme(legend.position = "none"),
               p_duration + theme(legend.position = "none"),
               get_legend(p_remain + theme(legend.position = "bottom")),
               ncol = 1,
               rel_heights = c(1,1,1,0.15,0.15))



abx_time = test_diffs[1]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx1_pha1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  geom_point(aes(x=47, y=10^11.5), pch = 21, fill = "black", size = 5)


abx_time = test_diffs[2]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx6_pha1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  geom_point(aes(x=47, y=10^11.5), pch = 22, fill = "black", size = 5)



abx_time = test_diffs[3]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx8_pha1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  geom_point(aes(x=47, y=10^11.5), pch = 8, fill = "black", size = 5)



abx_time = test_diffs[4]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx16_pha1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  geom_point(aes(x=47, y=10^11.5), pch = 24, fill = "black", size = 5)


final_plot = plot_grid(pa,
                       NULL,
                       pb,
                       NULL,
                       plot_grid(
                         plot_grid(p_abx1_pha1 + theme(legend.position = "none"),
                                   p_abx6_pha1 + theme(legend.position = "none"),
                                   p_abx8_pha1 + theme(legend.position = "none"),
                                   p_abx16_pha1 + theme(legend.position = "none"),
                                   ncol = 2),
                         get_legend(p_abx1_pha1 + theme(legend.position = "bottom")),
                         ncol = 1,
                         rel_heights = c(1,0.05)),
                       ncol = 5,
                       rel_widths = c(0.6,0.05,0.6,0.05,1.1),
                       labels = c("a)", "", "b)", "", "c)"))

ggsave(here::here("Figures", "suppfig4.png"), final_plot, height = 12, width = 17)



## same, but with BT in minority ######################################

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

times = seq(0, 48, 0.1)

yinit = c(Be = 1e6,
          Bt = 1e9,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

all_results = data.frame()

for(pha_con in c(10^7, 10^8, 10^9, 10^10)){
  
  results_con = data.frame(abx_start = c(rep(0,25), 1:24),
                           phage_start = c(0:24, rep(0,24)),
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
  geom_point(data = all_results %>% filter(pha_con == test_pha_con) %>% filter(diff %in% test_diffs),
             aes(diff, duration_bet, shape = as.factor(diff)), size = 3, fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Double-resistant bacteria presence time",
       colour = "Phage added\n(pfu/mL):",
       shape = "Condition:") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(limits = c(0, 47)) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+07", "1e+08", "1e+09", "1e+10"),
                     labels = c(bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10))) +
  scale_shape_manual(values = c(21,22,8,24),
                     breaks = as.factor(test_diffs),
                     labels = c("1.", "2.", "3.", "4.")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_max = ggplot(all_results) +
  geom_line(aes(diff, max_bet, colour = as.factor(pha_con)), size = 0.8) +
  geom_point(data = all_results %>% filter(pha_con == test_pha_con) %>% filter(diff %in% test_diffs),
             aes(diff, max_bet, shape = as.factor(diff)), size = 3, fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Maximum double-resistant bacteria",
       colour = "Phage added\n(pfu/mL):",
       shape = "Condition:") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+07", "1e+08", "1e+09", "1e+10"),
                     labels = c(bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10))) +
  scale_shape_manual(values = c(21,22,8,24),
                     breaks = as.factor(test_diffs),
                     labels = c("1.", "2.", "3.", "4.")) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_remain = ggplot(all_results) +
  geom_line(aes(diff, end_bacteria, colour = as.factor(pha_con)), size = 0.8) +
  geom_point(data = all_results %>% filter(pha_con == test_pha_con) %>% filter(diff %in% test_diffs),
             aes(diff, end_bacteria, shape = as.factor(diff)), size = 3, fill = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Total remaining bacteria after 48h",
       colour = "Phage added\n(pfu/mL):",
       shape = "Condition:") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette,#c("grey", "grey", "grey", "grey", "grey", "grey"), 
                     breaks = c("1e+07", "1e+08", "1e+09", "1e+10"),
                     labels = c(bquote(10^7),
                                bquote(10^8),
                                bquote(10^9),
                                bquote(10^10))) +
  scale_shape_manual(values = c(21,22,8,24),
                     breaks = as.factor(test_diffs),
                     labels = c("1.", "2.", "3.", "4.")) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


pb = plot_grid(p_remain + theme(legend.position = "none"),
               p_max + theme(legend.position = "none"),
               p_duration + theme(legend.position = "none"),
               get_legend(p_remain + guides(shape = F) + theme(legend.position = "bottom")),
               ncol = 1,
               rel_heights = c(1,1,1,0.15,0.15))


all_results_abx = data.frame()

for(abx_con in c(0.25, 0.5, 1, 2)){
  
  results_con = data.frame(abx_start = c(rep(0,25), 1:24),
                           phage_start = c(0:24, rep(0,24)),
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
  
  all_results_abx = rbind(all_results_abx, results_con)
  
}


p_duration = ggplot(all_results_abx) +
  geom_line(aes(diff, duration_bet, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Double-resistant bacteria presence time",
       colour = "Antibiotics added\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(limits = c(0, 47)) +
  scale_color_manual(values = palette) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_max = ggplot(all_results_abx) +
  geom_line(aes(diff, max_bet, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Maximum double-resistant bacteria",
       colour = "Antibiotics added\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  

p_remain = ggplot(all_results_abx) +
  geom_line(aes(diff, end_bacteria, colour = as.factor(abx_con)), size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Antibiotics addition time, relative to phage addition time",
       y = "Total remaining bacteria after 48h",
       colour = "Antibiotics added\n(mg/L):") +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n = 6),
                     labels = c("x", "0", expression(10^0),
                                expression(10^2), expression(10^4), expression(10^6),
                                expression(10^8), "x")) +
  scale_color_manual(values = palette) +
  coord_cartesian(ylim = c(0.01, 1e9)) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12))  


pa = plot_grid(p_remain + theme(legend.position = "none"),
               p_max + theme(legend.position = "none"),
               p_duration + theme(legend.position = "none"),
               get_legend(p_remain + theme(legend.position = "bottom")),
               ncol = 1,
               rel_heights = c(1,1,1,0.15,0.15))



abx_time = test_diffs[1]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx1_pha1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  geom_point(aes(x=47, y=10^11.5), pch = 21, fill = "black", size = 5)


abx_time = test_diffs[2]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx6_pha1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  geom_point(aes(x=47, y=10^11.5), pch = 22, fill = "black", size = 5)



abx_time = test_diffs[3]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx8_pha1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  geom_point(aes(x=47, y=10^11.5), pch = 8, fill = "black", size = 5)



abx_time = test_diffs[4]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, test_pha_con),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$Pl[results$Pl == 0] = NA

p_abx16_pha1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10, yend = 10^10) +
  geom_label(x = abx_time+8, y = 10, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+7.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = sum(results[481,c(2:4)]), linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
  scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl"),
                      values = c("#685cc4","#6db356","#c2484d","#c88a33"),
                      labels = c(expression(B[E]),
                                 expression(B[T]),
                                 expression(B[ET]),
                                 expression(P[L]))) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold")) +
  geom_point(aes(x=47, y=10^11.5), pch = 24, fill = "black", size = 5)


final_plot = plot_grid(pa,
                       NULL,
                       pb,
                       NULL,
                       plot_grid(
                         plot_grid(p_abx1_pha1 + theme(legend.position = "none"),
                                   p_abx6_pha1 + theme(legend.position = "none"),
                                   p_abx8_pha1 + theme(legend.position = "none"),
                                   p_abx16_pha1 + theme(legend.position = "none"),
                                   ncol = 2),
                         get_legend(p_abx1_pha1 + theme(legend.position = "bottom")),
                         ncol = 1,
                         rel_heights = c(1,0.05)),
                       ncol = 5,
                       rel_widths = c(0.6,0.05,0.6,0.05,1.1),
                       labels = c("a)", "", "b)", "", "c)"))

ggsave(here::here("Figures", "suppfig5.png"), final_plot, height = 12, width = 17)
