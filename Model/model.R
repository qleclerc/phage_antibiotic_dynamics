
library(deSolve)
library(ggplot2)
library(scales)

phage_tr_model <- function(parameters, init.state, times, event_dat) {
  
  model_dde <- function(time, state, parameters) {
    
    mu_e = parameters[["mu_e"]]
    mu_t = parameters[["mu_t"]]
    mu_et = parameters[["mu_et"]]
    Nmax = parameters[["Nmax"]]
    
    beta = parameters[["beta"]]
    L = parameters[["L"]]
    gamma = parameters[["gamma"]]
    alpha = parameters[["alpha"]]
    tau = parameters[["tau"]]
    
    ery_kill_max = parameters[["ery_kill_max"]]
    tet_kill_max = parameters[["tet_kill_max"]]
    EC_ery = parameters[["EC_ery"]]
    EC_tet = parameters[["EC_tet"]]
    gamma_tet = parameters[["gamma_tet"]]
    gamma_ery = parameters[["gamma_ery"]]
    
    Be = state[["Be"]]
    Bt = state[["Bt"]]
    Bet = state[["Bet"]]
    Pl = state[["Pl"]]
    Pe = state[["Pe"]]
    Pt = state[["Pt"]]
    
    ery = state[["ery"]]
    tet = state[["tet"]]
    
    N = Be + Bt + Bet
    
    if(time <= tau){
      Be_past = 0
      Bt_past = 0
      Bet_past = 0
      Pl_past = 0
      Pe_past = 0
      Pt_past = 0
      N_past = 1
      ery_past = 0
      tet_past = 0
    } else {
      Be_past = lagvalue(time - tau, 1)
      Bt_past = lagvalue(time - tau, 2)
      Bet_past = lagvalue(time - tau, 3)
      Pl_past = lagvalue(time - tau, 4)
      Pe_past = lagvalue(time - tau, 5)
      Pt_past = lagvalue(time - tau, 6)
      N_past = Be_past + Bt_past + Bet_past
      ery_past = lagvalue(time - tau, 7)
      tet_past = lagvalue(time - tau, 8)
    }
    
    link = (1 - N/Nmax)
    L = L * link + 1
    
    lambda = (1 - exp(-beta * N))
    phi_Pl = (1 - exp(-lambda * Pl/N))
    phi_Pe = (1 - exp(-lambda * Pe/N))
    phi_Pt = (1 - exp(-lambda * Pt/N))
    
    lambda_past = (1 - exp(-beta * N_past))
    phi_Pl_past = (1 - exp(-lambda_past * Pl_past/N_past))
    phi_Pe_past = (1 - exp(-lambda_past * Pe_past/N_past))
    phi_Pt_past = (1 - exp(-lambda_past * Pt_past/N_past))
    
    ery_effect = ery_kill_max * ery/(EC_ery+ery)
    ery_effect = 1 - ery_effect/mu_t
    tet_effect = tet_kill_max * tet/(EC_tet+tet)
    
    ery_effect_past = ery_kill_max * ery_past/(EC_ery+ery_past)
    ery_effect_past = 1 - ery_effect_past/mu_t
    tet_effect_past = tet_kill_max * tet_past/(EC_tet+tet_past)
    
    dBe = mu_e * link * (Be - ((phi_Pl + phi_Pt) * Be) ) - (phi_Pl + phi_Pt) * Be - tet_effect*Be
    dBt = (mu_t * link * ery_effect) * (Bt - ((phi_Pl + phi_Pe) * Bt) ) - (phi_Pl + phi_Pe) * Bt
    dBet = mu_et * link * (Bet - (phi_Pl*Bet) ) - phi_Pl * Bet +
      phi_Pe * Bt + phi_Pt * (Be - tet_effect*Be)
    
    dPl = phi_Pl_past * L * (1-alpha) * Be_past +
      phi_Pl_past * L * ery_effect_past * (1-alpha) * Bt_past +
      phi_Pl_past * L * (1-2*alpha) * Bet_past -
      lambda * Pl - gamma * Pl
    dPe = phi_Pl_past * L * alpha * Be_past +
      phi_Pl_past * L * alpha * Bet_past -
      lambda * Pe - gamma * Pe 
    dPt = phi_Pl_past * L * ery_effect_past * alpha * Bt_past +
      phi_Pl_past * L * alpha * Bet_past -
      lambda * Pt - gamma * Pt 
    
    dery = -ery*gamma_ery
    dtet = -tet*gamma_tet
    
    return(list(c(dBe, dBt, dBet, dPl, dPe, dPt, dery, dtet)))
  }
  
  trajectory <- data.frame(dede(y = init.state,
                                times = times,
                                func = model_dde,
                                parms = parameters,
                                events = list(data = event_dat)))
  
  return(trajectory)
  
}

parameters = c(mu_e = 1.609,
               mu_t = 1.513,
               mu_et = 1.445,
               Nmax = 2.76e9,
               beta = 1e-10,
               L = 50,
               tau = 0.67,
               alpha = 1e-6,
               gamma = 0.005,
               ery_kill_max = 1.513,
               tet_kill_max = 3,
               EC_ery = 0.8,
               EC_tet = 0.8,
               gamma_ery = 0.1,
               gamma_tet = 0.1)

times = seq(0, 24, 1)

yinit = c(Be = 1e8,
          Bt = 1e8,
          Bet = 0,
          Pl = 0,
          Pe = 0,
          Pt = 0,
          ery = 0,
          tet = 0)

event_dat = data.frame(var = c("ery", "tet", "Pl"),
                       time = c(100, 10, 10) ,
                       value = c(10, 10, 1e9),
                       method = c("add", "add", "add"))

results = phage_tr_model(parameters, yinit, times, event_dat)

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
        strip.text.x = element_text(size=12))
                      

