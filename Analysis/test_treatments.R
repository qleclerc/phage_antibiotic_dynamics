
## SETUP ##########

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
               alpha = 0,
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

times = seq(0, 25, 1)

yinit = c(Be = 1e9,
          Bt = 1e8,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)


## ERY ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 100, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

p_abx_ery = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "Ery only") +
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
        strip.text.x = element_text(size=12))


## TET ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 0, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

p_abx_tet = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "Tet only") +
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
        strip.text.x = element_text(size=12))


## ERY AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 0, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

p_abx_both = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "Ery and tet") +
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
        strip.text.x = element_text(size=12))



## PHAGE ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 100, 0) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

p_phage = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "Phage only") +
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
        strip.text.x = element_text(size=12))


## PHAGE AND ABX ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 0, 0) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)


p_both = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "Phage, ery and tet") +
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
        strip.text.x = element_text(size=12))


## PHAGE ONLY, TRANSDUCTION ##########

parameters[["alpha"]] = pha_params$alpha

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 100, 0) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

p_phage_tr = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "Phage only + transduction") +
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
        strip.text.x = element_text(size=12))


## PHAGE AND ABX, TRANSDUCTION ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 0, 0) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

p_both_tr = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "Phage, ery and tet + transduction") +
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
        strip.text.x = element_text(size=12))


## FINAL PLOT ##########

legend = get_legend(p_both)
plot_grid(p_abx_ery + theme(legend.position = "none"),
          p_abx_tet + theme(legend.position = "none"),
          p_abx_both + theme(legend.position = "none"),
          p_phage + theme(legend.position = "none"),
          p_both + theme(legend.position = "none"),
          legend,
          p_phage_tr + theme(legend.position = "none"),
          p_both_tr + theme(legend.position = "none"))

ggsave(here::here("Figures", "treatments.png"))
