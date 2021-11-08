library(deSolve)
library(dplyr)

data = data %>%
  filter(time < 8)

quick_model <- function(parameters, init.state, times) {
  
  model_ode <- function(time, state, parameters) {
    
    mu = parameters[["mu"]]
    Bmax = parameters[["Bmax"]]
    
    abx_effect = parameters[["abx_effect"]]
    
    B = state[["B"]]
    
    link = (1 - B/Bmax)
    
    dB = mu * link * B - abx_effect * mu * B
    
    return(list(dB))
  }
  
  trajectory <- data.frame(ode(y = init.state,
                               times = times,
                               func = model_ode,
                               parms = parameters))
  
  return(trajectory)
  
}


#fit to 0
data_0 = data %>%
  filter(Concentration == 0)

fit_0 = function(params, data_0){
  
  params = c(mu = params[1], Bmax = params[2], abx_effect = 0)
  results = quick_model(params, c(B=data_0$cfu[data_0$time==0]), data_0$time)
  
  results = merge(data_0, results, by = "time")
  
  sum((results$cfu-results$B)^2)
  
}

growth_par = optim(c(1.5, 1e9), fit_0, data_0 = data_0)$par


#fit to all
fit_all = function(abx_effect, growth_par, data_c){

  params = c(mu = growth_par[1], Bmax = growth_par[2], abx_effect = abx_effect)
  results = quick_model(params, c(B=data_c$cfu[data_c$time==0]), data_c$time)

  for(i in 1:nrow(results)){
    
    results$B[i] = results$B[i]/(10^(max(floor(log10(data_c$cfu[i])),1)-1))
    data_c$cfu[i] = data_c$cfu[i]/(10^(max(floor(log10(data_c$cfu[i])),1)-1))
    
    
  }
  
  results = merge(data_c, results, by = "time")
  
  sum((results$cfu-results$B)^2, na.rm = T)
  
}

summary_effect = data.frame(Concentration = unique(data$Concentration), effect = 0)

for(i in unique(data$Concentration)){
  
  if(i == 0) next
  
  data_c = data %>%
    filter(Concentration == i)
  
  c_par = optimise(fit_all,interval = c(0.01,10), 
                growth_par = growth_par, data_c = data_c)$minimum

  summary_effect$effect[summary_effect$Concentration == i] = c_par
  
}


#fit abx data
fit_abx = function(params, summary_effect){
  
  results = params[1] * (summary_effect$Concentration^params[3]/(summary_effect$Concentration^params[3]+params[2]^params[3]))
  
  results = cbind(summary_effect, results)
  
  sum((results$effect-results$results)^2)
  
}

abx_par = optim(c(2, 0.1, 1), fit_abx, summary_effect = summary_effect)$par

summary_effect$fitted_effect = abx_par[1] * (summary_effect$Concentration^abx_par[3]/(summary_effect$Concentration^abx_par[3]+abx_par[2]^abx_par[3]))

ggplot(summary_effect) +
  geom_line(aes(Concentration, effect, colour = "Data"), size = 1) +
  geom_point(aes(Concentration, effect, colour = "Data"), size = 3) +
  geom_line(aes(Concentration, fitted_effect, colour = "Fitted"), size = 1) +
  geom_point(aes(Concentration, fitted_effect, colour = "Fitted"), size = 3) +
  labs(x = "Concentration (mg/L)", y = "Kill rate per hour", colour = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12),
        legend.title = element_text(size=12))

ggsave(here::here("Figures", paste0(bac, "_", abx, "_effect.png")))

model_results = c()
#summary_effect$fitted_effect = c(0,1.55,1.93,2.41,2.95,3.54,4.15,4.77,5.36)
for(i in unique(summary_effect$Concentration)){
  
  data_c = data %>%
    filter(Concentration == i)
  
  params = c(mu = growth_par[1],
             Bmax = growth_par[2],
             abx_effect = summary_effect$fitted_effect[summary_effect$Concentration == i])
  
  results = quick_model(params, c(B=data_c$cfu[data_c$time==0]), data_c$time)
  
  results = cbind(rep(i, nrow(results)),
                  results,
                  rep("Model", nrow(results)))
  
  colnames(results) = c("Concentration", "time", "cfu", "source")
  
  model_results = rbind(model_results, results)
  
}

data$source = "Data"
model_results = rbind(data, model_results)

ggplot(model_results, aes(time, cfu, colour = as.factor(Concentration), linetype = source)) +
  geom_line(size=1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmax(0,cfu-se), ymax=cfu+se), width = 0.2, size = 1) +
  theme_bw() +
  labs(x = "Time (hours)", y = "cfu per mL", colour = "Antibiotic concentration (mg/L):", linetype = "Source:") +
  scale_color_manual(values=palette)+
  coord_cartesian(ylim=c(1e1,1e10))+
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0,24,2)) +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12),
        legend.title = element_text(size=12))

ggsave(here::here("Figures", paste0(bac, "_", abx, "_fitted",".png")))
