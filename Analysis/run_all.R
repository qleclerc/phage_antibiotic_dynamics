
all_params = data.frame()

for(bac in c("327", "201kt7")){
  for(abx in c("ery", "tet")){
    
    source(here::here("Analysis", "clean_plot_data.R"))
    source(here::here("Analysis", "estimate_parameters.R"))
    
    p = data.frame(bacteria = bac,
                   antibiotic = abx,
                   mu = growth_par[1],
                   Bmax = growth_par[2],
                   kmax = abx_par[1],
                   EC50 = abx_par[2],
                   pow = abx_par[3])
    
    all_params = rbind(all_params, p)
    
  }
}

write.csv(all_params, "all_params.csv", row.names = F)
