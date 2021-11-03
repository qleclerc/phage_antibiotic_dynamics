
## SETUP ##########

library(deSolve)
library(ggplot2)
library(scales)
library(cowplot)
library(dplyr)
library(reshape2)

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
          Bt = 1e6,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

all_results = data.frame()

## NOTHING ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 100, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## PHAGE ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 100, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## ERY ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 100, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Erythromycin"
all_results = rbind(all_results, results)


## TET ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 1, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Tetracycline"
all_results = rbind(all_results, results)


## ERY AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 1, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Erythromycin + tetracycline"
all_results = rbind(all_results, results)


## PHAGE AND ERY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 100, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Erythromycin"
all_results = rbind(all_results, results)


## PHAGE AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 1, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Tetracycline"
all_results = rbind(all_results, results)


## PHAGE AND ERY AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 1, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Erythromycin + tetracycline"
all_results = rbind(all_results, results)



## PHAGE ONLY, TRANSDUCTION ##########

parameters[["alpha"]] = pha_params$alpha

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 100, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## PHAGE AND ERY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 100, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Erythromycin"
all_results = rbind(all_results, results)


## PHAGE AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 1, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Tetracycline"
all_results = rbind(all_results, results)


## PHAGE AND ABX, TRANSDUCTION ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 1, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Erythromycin + tetracycline"
all_results = rbind(all_results, results)


## FINAL PLOT ##########

all_results$abx = as.factor(all_results$abx)
all_results$abx = factor(all_results$abx, levels = levels(all_results$abx)[c(3,1,4,2)])

all_results$Pl[all_results$Pl == 0] = NA

pa = ggplot(all_results) +
  geom_line(aes(time, Be, colour = "Be"), size = 1, alpha = 0.6) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 1, alpha = 0.6) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 1) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  facet_grid(pha ~ abx) +
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
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  theme(legend.position = "bottom")


all_results2 = all_results %>%
  filter(pha == "Phage w/ transduction") %>%
  select(-c("Pe", "Pt", "ery", "tet")) %>%
  reshape2::melt(. , id.vars = c("time","pha","abx"))

bac_labs = c("Ery-resistant", "Tet-resistant", "Double-resistant", "Phage")
names(bac_labs) = c("Be", "Bt", "Bet", "Pl")

pb = ggplot(all_results2) +
  geom_line(aes(time, value, colour = abx), size = 1, alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  facet_grid(pha ~ variable, labeller = labeller(variable = bac_labs)) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Antibiotic:") +
  scale_colour_manual(values = c("black", "royalblue3","green3", "purple3")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  theme(legend.position = "bottom")


plot_grid(pa,
          NULL,
          pb,
          labels = c("a)","", "b)"),
          nrow = 3, rel_heights = c(1,0.05,0.5))

ggsave(here::here("Figures", "fig4.png"), width = 10, height = 11)



### REPEAT FOR OTHER STRAIN ###########

yinit = c(Be = 1e6,
          Bt = 1e9,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

all_results = data.frame()

## NOTHING ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 100, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## PHAGE ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 100, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## ERY ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 100, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Erythromycin"
all_results = rbind(all_results, results)


## TET ONLY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 1, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Tetracycline"
all_results = rbind(all_results, results)


## ERY AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 1, 100) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "No phage added"
results$abx = "Erythromycin + tetracycline"
all_results = rbind(all_results, results)


## PHAGE AND ERY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 100, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Erythromycin"
all_results = rbind(all_results, results)


## PHAGE AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 1, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Tetracycline"
all_results = rbind(all_results, results)


## PHAGE AND ERY AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 1, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage added"
results$abx = "Erythromycin + tetracycline"
all_results = rbind(all_results, results)



## PHAGE ONLY, TRANSDUCTION ##########

parameters[["alpha"]] = pha_params$alpha

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 100, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "No antibiotic"
all_results = rbind(all_results, results)


## PHAGE AND ERY ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 100, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Erythromycin"
all_results = rbind(all_results, results)


## PHAGE AND TET ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 1, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Tetracycline"
all_results = rbind(all_results, results)


## PHAGE AND ABX, TRANSDUCTION ##########

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(1, 1, 1) ,
                       value = c(1, 1, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

results$pha = "Phage w/ transduction"
results$abx = "Erythromycin + tetracycline"
all_results = rbind(all_results, results)


## FINAL PLOT ##########

all_results$abx = as.factor(all_results$abx)
all_results$abx = factor(all_results$abx, levels = levels(all_results$abx)[c(3,1,4,2)])

all_results$Pl[all_results$Pl == 0] = NA

pa = ggplot(all_results) +
  geom_line(aes(time, Be, colour = "Be"), size = 1, alpha = 0.6) +
  geom_line(aes(time, Bt, colour = "Bt"), size = 1, alpha = 0.6) +
  geom_line(aes(time, Bet, colour = "Bet"), size = 1) +
  geom_line(aes(time, Pl, colour = "Pl"), size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  facet_grid(pha ~ abx) +
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
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  theme(legend.position = "bottom")


all_results2 = all_results %>%
  filter(pha == "Phage w/ transduction") %>%
  select(-c("Pe", "Pt", "ery", "tet")) %>%
  reshape2::melt(. , id.vars = c("time","pha","abx"))

bac_labs = c("Ery-resistant", "Tet-resistant", "Double-resistant", "Phage")
names(bac_labs) = c("Be", "Bt", "Bet", "Pl")

pb = ggplot(all_results2) +
  geom_line(aes(time, value, colour = abx), size = 1, alpha = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  facet_grid(pha ~ variable, labeller = labeller(variable = bac_labs)) +
  coord_cartesian(ylim = c(0.1, 3e11)) +
  scale_x_continuous(breaks=seq(0,max(results$time),4))+
  theme_bw() +
  labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Antibiotic:") +
  scale_colour_manual(values = c("black", "royalblue3","green3", "purple3")) +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12),
        legend.title = element_text(size=12)) +
  theme(legend.position = "bottom")


plot_grid(pa,
          NULL,
          pb,
          labels = c("a)","", "b)"),
          nrow = 3, rel_heights = c(1,0.05,0.5))

ggsave(here::here("Figures", "supp_fig3.png"), width = 10, height = 11)