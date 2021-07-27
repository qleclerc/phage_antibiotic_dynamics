model_name <- "transduction model"
model_state.names <- c("Be","Bt","Bet","Pl", "Pe", "Pt")
model_theta.names <- c("beta", "L", "gamma", "alpha", "tau")

model_simulateDeterministic <- function(theta,init.state,times) {
  
}


## function to compute log-prior
model_prior <- function(theta, log = TRUE) {
  
  #log.prior.L <- dunif(theta[["L"]], min = 1, max = 500, log = TRUE)
  log.prior.L <- dnorm(theta[["L"]], mean = 40, sd = 7, log = TRUE)
  log.prior.beta <- dunif(theta[["beta"]], min = 1, max = 1e20, log = TRUE)
  log.prior.gamma <- dunif(theta[["gamma"]], min = 1, max = 1e10, log = TRUE)
  log.prior.alpha <- dunif(theta[["alpha"]], min = 1, max = 1e10, log = TRUE)
  
  #log.prior.tau <- dunif(theta[["tau"]], min = 0.01, max = 0.8, log = TRUE)
  log.prior.tau <- dnorm(theta[["tau"]], mean = 0.67, sd = 0.07, log = TRUE)
  
  log.sum <- log.prior.L + log.prior.beta + log.prior.gamma + log.prior.alpha + log.prior.tau
  
  return(ifelse(log, log.sum, exp(log.sum)))
}

## function to compute the likelihood of one data point
model_pointLike <- function(data.point, model.point, theta, log = FALSE){
  
  # dpoisBe = dpois(x = data.point[["Be"]],
  #                 lambda = model.point[["Be"]],
  #                 log = log)
  # 
  # dpoisBt = dpois(x = data.point[["Bt"]],
  #                 lambda = model.point[["Bt"]],
  #                 log = log)
  dpoisBe = dpoisBt = dpoisPe = dpoisPt = 0
  
  dpoisBet = dpois(x = data.point[["Bet"]],
                   lambda = model.point[["Bet"]],
                   log = log)
  if(is.infinite(dpoisBet)) dpoisBet = -1e7
  
  
  dpoisPl = dpois(x = round(data.point[["P"]]/(10^(max(nchar(as.character(round(model.point[["Pl"]]))),2)-2))),
                  lambda = model.point[["Pl"]]/(10^(max(nchar(as.character(round(model.point[["Pl"]]))),2)-2)),
                  log = log)
  
  ## the prevalence is observed through a Poisson process
  return(sum(dpoisBe, dpoisBt, dpoisBet, dpoisPl, dpoisPe, dpoisPt))
}

## function to generate observation from a model simulation
# phagebac_genObsPoint <- function(model.point, theta){
#   
#   ## the prevalence is observed through a Poisson process
#   obs.point <- rpois(n = 1, lambda = model.point[["I"]])
#   
#   return(c(obs = obs.point))
# }

## create deterministic SIR fitmodel
model <- fitR::fitmodel(
  name = model_name,
  state.names = model_state.names,
  theta.names = model_theta.names,
  simulate = model_simulateDeterministic,
  dprior = model_prior,
  dPointObs = model_pointLike)

saveRDS(model, here::here("Model", "transduction_model.rds"))
