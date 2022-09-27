

library(deSolve)
library(ggplot2)
library(ggtext)
library(scales)
library(cowplot)
library(RColorBrewer)
library(dplyr)

source(here::here("Model", "new_model.R"))

#things to test
test_pha_con = 1e8
test_abx_con = 1
test_diffs = c(0,3,5,15)

#palette = rev(brewer.pal(n = 9, name = "RdBu"))[c(1,3,7,9)]
palette_blues = brewer.pal(n = 7, name = "Blues")
palette_greens = brewer.pal(n = 7, name = "Greens")

abx_params = read.csv(here::here("Parameters", "abx_params.csv"))
pha_params = read.csv(here::here("Parameters", "pha_params_new.csv"))
bac_params = read.csv(here::here("Parameters", "bac_params.csv"))

parameters = c(mu_e = bac_params$mu_e[1],
               mu_t = bac_params$mu_t[1],
               mu_et = bac_params$mu_et[1],
               Nmax = bac_params$Nmax[1],
               beta = pha_params$beta,
               L = pha_params$L,
               tau = pha_params$tau,
               alpha = pha_params$alpha,
               P50 = pha_params$P50,
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
          Bt = 1e9,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

all_results_pha = data.frame()

for(pha_con in c(10^5, 5*10^5,
                 10^6, 5*10^6,
                 10^7, 5*10^7,
                 10^8, 5*10^8,
                 10^9, 5*10^9,
                 10^10)){
  
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
    
    results_con$max_bet[i] = max(results$Bet, na.rm = T)
    
    results_con$duration_bet[i] = sum(results$Bet>1)/10
    
    results_con$end_bacteria[i] = max(tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1),0)
    
  }
  
  results_con$diff = results_con$abx_start - results_con$phage_start
  results_con$pha_con = pha_con
  
  all_results_pha = rbind(all_results_pha, results_con)
  
}

all_results_pha_cat = all_results_pha %>%
  mutate(duration_bet_cat = cut(duration_bet, breaks = c(0, 0.99, 10, 20, 30, 40, 50),
                                labels = c("0", "1 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50"),
                                include.lowest = T)) %>%
  mutate(max_bet_cat = cut(max_bet, breaks = c(0, 10^0, 10^1, 10^2, 10^4, 10^6, 10^8, 10^10),
                           labels = c("0", "1 - 10", "10 - 10^2", "10^2 - 10^4", "10^4 - 10^6",
                                      "10^6 - 10^8", "10^8 - 10^10"),
                           include.lowest = T)) %>%
  mutate(end_bacteria_cat = cut(end_bacteria, breaks = c(0, 10^0, 10^1, 10^2, 10^4, 10^6, 10^8, 10^10),
                                labels = c("0", "1 - 10","10 - 10^2", "10^2 - 10^4", "10^4 - 10^6",
                                           "10^6 - 10^8", "10^8 - 10^10"),
                                include.lowest = T))

all_results_pha_cat = all_results_pha_cat %>%
  mutate(expanded = (pha_con == test_pha_con & diff %in% test_diffs)) %>%
  mutate(pha_con = gsub("1e+", "10^", pha_con, fixed = T)) %>%
  mutate(pha_con = gsub("5e+", "5*10^", pha_con, fixed = T)) %>%
  mutate(pha_con = gsub("^0", "^", pha_con, fixed = T)) %>%
  mutate(pha_con = factor(pha_con, levels = unique(pha_con)))



p_duration = ggplot(all_results_pha_cat[order(all_results_pha_cat$expanded),],
                    aes(x = diff, y = pha_con,
                        fill = as.factor(duration_bet_cat))) +
  geom_tile(aes(color = expanded), size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.y = element_markdown(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_fill_manual(values = palette_blues) +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  labs(x = "Antibiotics addition time, relative to phage",
       y = "Phage (pfu/mL)",
       fill = "Double-resistant bacteria\npresence time (h)") +
  scale_color_manual(guide = "none", values = c(`TRUE` = "black", `FALSE` = "NA")) +
  scale_y_discrete(labels = c("10^5", "",
                              "10^6", "",
                              "10^7", "",
                              "**10^8**", "",
                              "10^9", "",
                              "10^10"))

p_max = ggplot(all_results_pha_cat[order(all_results_pha_cat$expanded),],
               aes(x = diff, y = pha_con,
                   fill = as.factor(max_bet_cat))) +
  geom_tile(aes(color = expanded), size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.y = element_markdown(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_fill_manual(values = palette_blues) +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  labs(x = "Antibiotics addition time, relative to phage",
       y = "Phage (pfu/mL)",
       fill = "Maximum double-resistant\nbacteria (cfu/mL)") +
  scale_color_manual(guide = "none", values = c(`TRUE` = "black", `FALSE` = "NA")) +
  scale_y_discrete(labels = c("10^5", "",
                              "10^6", "",
                              "10^7", "",
                              "**10^8**", "",
                              "10^9", "",
                              "10^10"))

p_remain = ggplot(all_results_pha_cat[order(all_results_pha_cat$expanded),],
                  aes(x = diff, y = pha_con,
                      fill = as.factor(end_bacteria_cat))) +
  geom_tile(aes(color = expanded), size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.y = element_markdown(size=12),
        axis.text.x = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_fill_manual(values = palette_blues) +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  labs(x = "Antibiotics addition time, relative to phage",
       y = "Phage (pfu/mL)",
       fill = "Total bacteria\nremaining (cfu/mL)") +
  scale_color_manual(guide = "none", values = c(`TRUE` = "black", `FALSE` = "NA")) +
  scale_y_discrete(labels = c("10^5", "",
                              "10^6", "",
                              "10^7", "",
                              "**10^8**", "",
                              "10^9", "",
                              "10^10"))


pb = plot_grid(p_remain,
               p_max,
               p_duration,
               ncol = 1, align = "v")


all_results_abx = data.frame()

for(abx_con in seq(0.2, 2.2, 0.2)){
  
  results_con = data.frame(abx_start = c(rep(0,25), 1:24),
                           phage_start = c(0:24, rep(0,24)),
                           max_bet = 0,
                           duration_bet = 0,
                           end_bacteria = 0)
  
  for(i in 1:nrow(results_con)){
    
    event_dat = data.frame(var = c("ery", "tet", "Pl"),
                           time = c(results_con$abx_start[i], results_con$abx_start[i], results_con$phage_start[i]),
                           value = c(abx_con, abx_con, 1e8),
                           method = c("add", "add", "add"))
    
    #times = seq(0, max(24+results_con$abx_start[i], 24+results_con$phage_start[i]), 1)
    
    results = phage_tr_model(parameters, yinit, times, event_dat)
    
    results_con$max_bet[i] = max(results$Bet, na.rm = T)
    
    results_con$duration_bet[i] = sum(results$Bet>1)/10
    
    results_con$end_bacteria[i] = max(tail(results$Be,1) + tail(results$Bt,1) + tail(results$Bet,1),0)
    
  }
  
  results_con$diff = results_con$abx_start - results_con$phage_start
  results_con$abx_con = abx_con
  
  all_results_abx = rbind(all_results_abx, results_con)
  
}


all_results_abx_cat = all_results_abx %>%
  mutate(duration_bet_cat = cut(duration_bet, breaks = c(0, 0.99, 10, 20, 30, 40, 50),
                                labels = c("0", "1 - 10", "10 - 20", "20 - 30", "30 - 40", "40 - 50"),
                                include.lowest = T)) %>%
  mutate(max_bet_cat = cut(max_bet, breaks = c(0, 10^0, 10^1, 10^2, 10^4, 10^6, 10^8, 10^10),
                           labels = c("0", "1 - 10", "10 - 10^2", "10^2 - 10^4", "10^4 - 10^6",
                                      "10^6 - 10^8", "10^8 - 10^10"),
                           include.lowest = T)) %>%
  mutate(end_bacteria_cat = cut(end_bacteria, breaks = c(0, 10^0, 10^1, 10^2, 10^4, 10^6, 10^8, 10^10),
                                labels = c("0", "1 - 10", "10 - 10^2", "10^2 - 10^4", "10^4 - 10^6",
                                           "10^6 - 10^8", "10^8 - 10^10"),
                                include.lowest = T))

all_results_abx_cat = all_results_abx_cat %>%
  mutate(expanded = (abx_con == test_abx_con & diff %in% test_diffs)) %>%
  mutate(abx_con = factor(abx_con, levels = unique(abx_con)))


p_duration = ggplot(all_results_abx_cat[order(all_results_abx_cat$expanded),],
                    aes(x = diff, y = abx_con,
                        fill = as.factor(duration_bet_cat))) +
  geom_tile(aes(color = expanded), size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_markdown(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_fill_manual(values = palette_greens) +
  scale_color_manual(guide = "none", values = c(`TRUE` = "black", `FALSE` = "NA")) +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_discrete(labels = c("0.2", "", "0.6", "", "**1**",
                              "", "1.4", "", "1.8", "", "2.2")) +
  labs(x = "Antibiotics addition time, relative to phage",
       y = "Antibiotics (mg/L)",
       fill = "Double-resistant bacteria\npresence time (h)")

p_max = ggplot(all_results_abx_cat[order(all_results_abx_cat$expanded),],
               aes(x = diff, y = abx_con,
                   fill = as.factor(max_bet_cat))) +
  geom_tile(aes(color = expanded), size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_markdown(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_fill_manual(values = palette_greens) +
  scale_color_manual(guide = "none", values = c(`TRUE` = "black", `FALSE` = "NA")) +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_discrete(labels = c("0.2", "", "0.6", "", "**1**",
                              "", "1.4", "", "1.8", "", "2.2")) +
  labs(x = "Antibiotics addition time, relative to phage",
       y = "Antibiotics (mg/L)",
       fill = "Maximum double-resistant\nbacteria (cfu/mL)")

p_remain = ggplot(all_results_abx_cat[order(all_results_abx_cat$expanded),],
                  aes(x = diff, y = abx_con,
                      fill = as.factor(end_bacteria_cat))) +
  geom_tile(aes(color = expanded), size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_markdown(size=12),
        axis.title = element_text(size=12),
        legend.text = element_markdown(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  scale_fill_manual(values = palette_greens) +
  scale_color_manual(guide = "none", values = c(`TRUE` = "black", `FALSE` = "NA")) +
  scale_x_continuous(breaks = seq(-24,24,4)) +
  scale_y_discrete(labels = c("0.2", "", "0.6", "", "**1**",
                              "", "1.4", "", "1.8", "", "2.2")) +
  labs(x = "Antibiotics addition time, relative to phage",
       y = "Antibiotics (mg/L)",
       fill = "Total bacteria\nremaining (cfu/mL)")


pa = plot_grid(p_remain,
               p_max,
               p_duration,
               ncol = 1, align = "v")

ptop = plot_grid(pa, pb,
                 nrow = 1,
                 labels = c("a)", "b)"),
                 hjust = 0, vjust = 1)




#expanded view
abx_time = test_diffs[1]
pha_time = 0
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(test_abx_con, test_abx_con, test_pha_con),
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
  geom_label(x = pha_time+7, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_text(y = 4, x = 1.5, label = "0h", size = 4) +
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
        plot.title = element_text(face = "bold"))
#geom_point(aes(x=47, y=10^11.5), pch = 21, fill = "black", size = 5)


abx_time = test_diffs[2]
pha_time = 0 
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(test_abx_con, test_abx_con, test_pha_con),
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
  geom_label(x = pha_time+7, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_text(y = 4, x = 1.5, label = "3h", size = 4) +
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
        plot.title = element_text(face = "bold"))
#geom_point(aes(x=47, y=10^11.5), pch = 22, fill = "black", size = 5)



abx_time = test_diffs[3]
pha_time = 0 
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(test_abx_con, test_abx_con, test_pha_con),
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
  geom_label(x = pha_time+7, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_text(y = 4, x = 2.5, label = "5h", size = 4) +
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
        plot.title = element_text(face = "bold"))
#geom_point(aes(x=47, y=10^11.5), pch = 8, fill = "black", size = 5)



abx_time = test_diffs[4]
pha_time = 0 
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(test_abx_con, test_abx_con, test_pha_con),
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
  geom_label(x = pha_time+7, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_text(y = 4, x = 7.5, label = "15h", size = 4) +
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
        plot.title = element_text(face = "bold"))
#geom_point(aes(x=47, y=10^11.5), pch = 24, fill = "black", size = 5)


final_plot = plot_grid(ptop,
                       NULL,
                       plot_grid(
                         plot_grid(p_abx1_pha1 + theme(legend.position = "none"),
                                   p_abx6_pha1 + theme(legend.position = "none"),
                                   p_abx8_pha1 + theme(legend.position = "none"),
                                   p_abx16_pha1 + theme(legend.position = "none"),
                                   ncol = 2),
                         get_legend(p_abx1_pha1 + theme(legend.position = "right")),
                         nrow = 1,
                         rel_widths = c(1,0.1)),
                       ncol = 1,
                       rel_heights = c(1.3, 0.05, 1),
                       labels = c("", "" ,"c)"), hjust = 0)

ggsave(here::here("Figures", "fig5.png"), final_plot, height = 15, width = 13)

