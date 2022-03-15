
library(deSolve)
library(ggplot2)
library(scales)

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

times = seq(0, 48, 0.1)

yinit = c(Be = 1e9,
          Bt = 1e9,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

abx_time = 5
pha_time = 0
extra_label = "" #"morep"
event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(abx_time, abx_time, pha_time),
                       value = c(1, 1, 1e8),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat, T)
results$Pl[results$Pl == 0] = NA

ggplot(results) +
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
        strip.text.x = element_text(size=12)) +
  geom_vline(xintercept = abx_time, linetype = "dashed") +
  annotate("segment", x = abx_time+5, xend = abx_time, y = 10^10.5, yend = 10^10.5) +
  geom_label(x = abx_time+6, y = 10.5, label = "Antibiotics +", size = 3) +
  geom_vline(xintercept = pha_time, linetype = "dashed") +
  annotate("segment", x = pha_time+4.5, xend = pha_time, y = 10^11.5, yend = 10^11.5) +
  geom_label(x = pha_time+5.5, y = 11.5, label = "Phage +", size = 3) +
  geom_hline(yintercept = max(results$Bet, na.rm = T), linetype = "dotted") +
  annotate("segment", y = 10^(log10(max(results$Bet, na.rm = T))+1),
           yend = max(results$Bet), x = 40, xend = 40) +
  geom_label(y = log10(max(results$Bet, na.rm = T))+1,
             x = 40, label = "Max DRP", size = 3) +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey", size = 0.8)
  
ggsave(here::here("test_plots", paste0("a", abx_time, "p", pha_time, extra_label,".png")))

# ggplot(results) +
#   geom_line(aes(time, Be, colour = "Be"), size = 0.8) +
#   geom_line(aes(time, Bt, colour = "Bt"), size = 0.8) +
#   geom_line(aes(time, Bet, colour = "Bet"), size = 0.8) +
#   geom_line(aes(time, Pl, colour = "Pl"), size = 0.8) +
#   geom_line(aes(time, Bet_lysis, colour = "Bet_l"), size = 0.8) +
#   geom_line(aes(time, Bet_new, colour = "Bet_n"), size = 0.8) +
#   scale_y_continuous(trans=log10_trans(),
#                      breaks=trans_breaks("log10", function(x) 10^x),
#                      labels=trans_format("log10", math_format(10^.x))) +
#   coord_cartesian(ylim = c(0.1, 3e11)) +
#   scale_x_continuous(breaks=seq(0,max(results$time),4))+
#   theme_bw() +
#   labs(y = "cfu or pfu per mL", x = "Time (hours)", colour = "Organism:") +
#   scale_colour_manual(breaks = c("Be", "Bt", "Bet", "Pl","Bet_l", "Bet_n"),
#                       values = c("#685cc4","#6db356","#c2484d","#c88a33", "black", "blue"),
#                       labels = c(expression(B[E]),
#                                  expression(B[T]),
#                                  expression(B[ET]),
#                                  expression(P[L]),
#                                  "Bet_l", "Bet_n")) +
#   theme(axis.text.x = element_text(size=12),
#         axis.title.x = element_text(size=12),
#         axis.text.y = element_text(size=12),
#         axis.title.y = element_text(size=12),
#         legend.text = element_text(size=12),
#         legend.title = element_text(size=12),
#         strip.text.x = element_text(size=12))
