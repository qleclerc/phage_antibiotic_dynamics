
library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)

source(here::here("Model", "new_model.R"))

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
               alpha = 0,
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
                       time = c(0, 0, 0) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))
# value = c(5, 1.2, 1e9)
results = phage_tr_model(parameters, yinit, times, event_dat)
results$Pl[results$Pl == 0] = NA

p1 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 5e9)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "a)") +
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

optim_concentrations = function(concentrations, target){
  event_dat = data.frame(var = c("ery", "tet", "Pl"),
                         time = c(0, 0, 100) ,
                         value = c(concentrations[1], concentrations[2], 1e9),
                         method = c("add", "add", "add"))
  # value = c(5, 1.2, 1e9)
  results = phage_tr_model(parameters, yinit, times, event_dat) %>%
    filter(times > 6)

  sum((results$Be+results$Bt - target)^2)
  
}

results = results %>% filter(times > 6)
opti_concentrations = optim(c(1.1,1.1),
                            optim_concentrations,
                            target = results$Be+results$Bt,
                            lower = c(1,1), upper = c(2,2))


event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(0, 0, 100) ,
                       value = c(opti_concentrations$par[1], opti_concentrations$par[2], 1e9),
                       method = c("add", "add", "add"))

cat(opti_concentrations$par[1], opti_concentrations$par[2])

results = phage_tr_model(parameters, yinit, times, event_dat)

p2 = ggplot(results) +
  geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 5e9)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:", title = "b)") +
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

plot_grid(plot_grid(p1 + theme(legend.position = "none"),
                    p2 + theme(legend.position = "none")),
          plot_grid(NULL,get_legend(p1 + theme(legend.position = "bottom")), NULL,
                    nrow = 1),
          nrow = 2,
          rel_heights = c(1,0.1))

ggsave(here::here("Figures", "suppfig2.png"))
