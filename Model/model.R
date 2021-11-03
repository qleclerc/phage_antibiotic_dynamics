
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
    
    ery_kill_max_BE = parameters[["ery_kill_max_BE"]]
    tet_kill_max_BE = parameters[["tet_kill_max_BE"]]
    EC_ery_BE = parameters[["EC_ery_BE"]]
    EC_tet_BE = parameters[["EC_tet_BE"]]
    pow_ery_BE = parameters[["pow_ery_BE"]]
    pow_tet_BE = parameters[["pow_tet_BE"]]
    
    ery_kill_max_BT = parameters[["ery_kill_max_BT"]]
    tet_kill_max_BT = parameters[["tet_kill_max_BT"]]
    EC_ery_BT = parameters[["EC_ery_BT"]]
    EC_tet_BT = parameters[["EC_tet_BT"]]
    pow_ery_BT = parameters[["pow_ery_BT"]]
    pow_tet_BT = parameters[["pow_tet_BT"]]
    
    ery_kill_max_BET = parameters[["ery_kill_max_BET"]]
    tet_kill_max_BET = parameters[["tet_kill_max_BET"]]
    EC_ery_BET = parameters[["EC_ery_BET"]]
    EC_tet_BET = parameters[["EC_tet_BET"]]
    pow_ery_BET = parameters[["pow_ery_BET"]]
    pow_tet_BET = parameters[["pow_tet_BET"]]
    
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
    
    lambda = (1 - exp(-beta * N))
    phi_Pl = (1 - exp(-lambda * Pl/N))
    phi_Pe = (1 - exp(-lambda * Pe/N))
    phi_Pt = (1 - exp(-lambda * Pt/N))
    
    lambda_past = (1 - exp(-beta * N_past))
    phi_Pl_past = (1 - exp(-lambda_past * Pl_past/N_past))
    phi_Pe_past = (1 - exp(-lambda_past * Pe_past/N_past))
    phi_Pt_past = (1 - exp(-lambda_past * Pt_past/N_past))
    
    ery_effect_BE = ery_kill_max_BE * ery^pow_ery_BE/(EC_ery_BE^pow_ery_BE+ery^pow_ery_BE)
    tet_effect_BE = tet_kill_max_BE * tet^pow_tet_BE/(EC_tet_BE^pow_tet_BE+tet^pow_tet_BE)
    ery_effect_BT = ery_kill_max_BT * ery^pow_ery_BT/(EC_ery_BT^pow_ery_BT+ery^pow_ery_BT)
    tet_effect_BT = tet_kill_max_BT * tet^pow_tet_BT/(EC_tet_BT^pow_tet_BT+tet^pow_tet_BT)
    ery_effect_BET = ery_kill_max_BET * ery^pow_ery_BET/(EC_ery_BET^pow_ery_BET+ery^pow_ery_BET)
    tet_effect_BET = tet_kill_max_BET * tet^pow_tet_BET/(EC_tet_BET^pow_tet_BET+tet^pow_tet_BET)
    
    ery_effect_past_BE = ery_kill_max_BE * ery_past^pow_ery_BE/(EC_ery_BE^pow_ery_BE+ery_past^pow_ery_BE)
    tet_effect_past_BE = tet_kill_max_BE * tet_past^pow_tet_BE/(EC_tet_BE^pow_tet_BE+tet_past^pow_tet_BE)
    ery_effect_past_BT = ery_kill_max_BT * ery_past^pow_ery_BT/(EC_ery_BT^pow_ery_BT+ery_past^pow_ery_BT)
    tet_effect_past_BT = tet_kill_max_BT * tet_past^pow_tet_BT/(EC_tet_BT^pow_tet_BT+tet_past^pow_tet_BT)
    ery_effect_past_BET = ery_kill_max_BET * ery_past^pow_ery_BET/(EC_ery_BET^pow_ery_BET+ery_past^pow_ery_BET)
    tet_effect_past_BET = tet_kill_max_BET * tet_past^pow_tet_BET/(EC_tet_BET^pow_tet_BET+tet_past^pow_tet_BET)
    
    L_ET = L * max(0, (link - tet_effect_BET - ery_effect_BET)) + 1
    L_E = L * max(0, (link - tet_effect_BE - ery_effect_BE)) + 1
    L_T = L * max(0, (link - tet_effect_BT - ery_effect_BT)) + 1
    
    
    dBe = (mu_e * link) * (Be - ((phi_Pl + phi_Pt) * Be) ) -
      (phi_Pl + phi_Pt) * Be -  (tet_effect_BE + ery_effect_BE)*mu_e*Be
    dBt = (mu_t * link) * (Bt - ((phi_Pl + phi_Pe) * Bt) ) -
      (phi_Pl + phi_Pe) * Bt - (ery_effect_BT + tet_effect_BT)*mu_t*Bt
    dBet = mu_et * link * (Bet - (phi_Pl*Bet) ) - phi_Pl * Bet +
      phi_Pe * Bt + phi_Pt * Be - (ery_effect_BET + tet_effect_BET)*mu_et*Bet
    
    # if(time %% 1 == 0) cat(time, lambda, "\n")
    
    dPl = phi_Pl_past * L_E * (1-alpha) * Be_past +
      phi_Pl_past * L_T * (1-alpha) * Bt_past +
      phi_Pl_past * L_ET * (1-2*alpha) * Bet_past -
      lambda * Pl - gamma * Pl
    dPe = phi_Pl_past * L_E * alpha * Be_past +
      phi_Pl_past * L_ET * alpha * Bet_past -
      lambda * Pe - gamma * Pe 
    dPt = phi_Pl_past * L_T * alpha * Bt_past +
      phi_Pl_past * L_ET * alpha * Bet_past -
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

