

# WARNING
# Solving differential equations is not a clean thing!
# Computers struggle with really tiny numbers, so you might see weird things like -1e-20 bacteria or phage
# Generally, you can treat anything with an absolute value smaller than 1e-15 as equal to 0


#packages required
#if anything doesn't load, package is not installed!
#in that case, use install.package(that_package_name) once, then it'll always work later
library(deSolve)
library(ggplot2)
library(scales)


#function containing AND solving the model
phage_model <- function(parameters,               #the parameter values
                        init.state,               #initial numbers of phage and bacteria
                        times,                    #the times for which to solve the model
                        frequency = FALSE,        #if the model should be frequency dependent (is density otherwise)
                        link_adsorption = FALSE,  #if adsorption rate should be linked to growth rate
                        link_burst = FALSE,       #if burst size should be linked to growth rate
                        tell_me_more = FALSE) {   #should the model spit out extra info at each time step
  
  
  #function containing the model
  model_dde <- function(time,             #time step
                        state,            #states (ie compartments, ie bacteria and phage)
                        parameters) {     #parameter values
  
    #first we recover each parameter value stored in the "parameters" vector
    mu_max = parameters[["mu_max"]]
    Nmax = parameters[["Nmax"]]
    beta = parameters[["beta"]]
    delta = parameters[["delta"]]
    gamma = parameters[["gamma"]]
    tau = parameters[["tau"]]

    
    #then we recover the current bacteria and phage numbers stored in the "state" vector
    B = state[["B"]]
    P = state[["P"]]

    
    #remember, there's a latent period tau between phage infection and burst!
    #so, at each time t, the number of bacteria bursting depends on infections that happened at time t minus tau
    #to recover this, we need to recover the number of bacteria and phage we had at time t-tau
    if(time <= tau){
      #if time<tau, there has not yet been any phage infection cycle!
      B_past = 0
      P_past = 0
    } else {
      B_past = lagvalue(time - tau, 1)
      P_past = lagvalue(time - tau, 2)
    }
    
    
    #here we work out the link function to decrease parameters as bacteria reach stationary phase
    link = (1 - B/Nmax)
    #we also have to work out this value for time t-tau !
    #hopefully you can understand why, but let me know if that's not clear otherwise :)
    link_past = (1 - B_past/Nmax)

    
    #if we link adsorption to growth rate, the value of beta should decrease as bacteria reach stationary phase
    if(link_adsorption){
      beta = beta * link
      beta_past = beta * link_past
    } else beta_past = beta
    
    
    #if we link burst to growth rate, the value of delta should decrease as bacteria reach stationary phase
    #can you figure out why we don't care about delta_past though? ;)
    if(link_burst) delta = delta * link
    
    
    #here we define lambda and phi depending on whether we use frequency or density dependent interaction
    #here I define:
    # - lambda: NUMBER of phage binding to a bacteria
    # - phi: NUMBER of bacteria being infected by at least one phage
    if(frequency){
      
      #first work out number of phage binding to any bacteria
      lambda = (1 - exp(-beta * B)) * P
      #then number of bacteria infected by at least one phage
      phi = (1 - exp(-lambda/B)) * B
      
      #however, because of latent period, we need to recover phi at t-tau
      #since phi depends on lambda, we also have to recover lambda at t-tau
      lambda_past = (1 - exp(-beta_past * B)) * P
      phi_past = (1 - exp(-lambda_past/B)) * B
      
      #finally, we set omega to 1 to apply the growth correction I explained
      omega = 1
      
    } else {
      
      #in density dependent interaction, lambda and phi are the same!
      lambda = phi = beta * P * B
      #however, because of latent period, we need to recover phi at t-tau
      phi_past = beta_past * P_past * B_past
      
      #no need for growth correction here!
      omega = 0
    }
    
    
    #actual equations here
    dB = mu_max * link * (B - omega*phi) - phi
    #change in B = growth - infection at time t
    
    dP = delta * phi_past - lambda - gamma * P
    #change in P = burst from infections at time t-tau - infection at time t - decay
    

    #bit of code to spit out info if you're curious (see below)
    if(tell_me_more & time%%1 == 0){
      cat("\n\nTime is", time, ". Current values:",
          "\n- adsorption rate:", beta,
          "\n- burst size:", delta,
          "\n- number of phage:", format(P, scientific=T, digits = 3),
          "\n- number of bacteria:", format(B, scientific=T, digits = 3),
          "\n- number of new infections:", format(phi, scientific=T, digits = 3))
    }
    
    #return the change in B and P for that time step
    return(list(c(dB, dP)))
    
  }
  
  #function solving the model
  trajectory <- data.frame(dede(y = init.state,
                                times = times,
                                func = model_dde,
                                parms = parameters))
  
  #return the model results
  return(trajectory)
  
}


#okay now we have the model function, we need the parameters, times, and starting values!


#parameter values:
#you can experiment with these!
#test your understanding by first deciding what you think should happen if you increase something,
# then see what happens when you do it!
parameters = c(mu_max = 1.5, #max growth rate
               Nmax = 1e9,   #carrying capacity
               beta = 1e-10,  #phage adsorption rate
               delta = 100,  #phage burst size
               tau = 0.67,   #phage latent period (in hours)
               gamma = 0.003)  #phage decay rate


#times (in hours):
#again you can change this for longer periods of time
#seq(starting time, end time, time-step)
times = seq(0, 24, 1)


#starting values:
#again, feel free to change these!
yinit = c(B = 1e4,
          P = 1e4)


#aand run the model:
results = phage_model(parameters, yinit, times,
                      frequency = FALSE,
                      link_adsorption = FALSE,
                      link_burst = FALSE,
                      tell_me_more = FALSE)

# ! ABOUT TELL_ME_MORE !
#if you set that to TRUE, model will spit out some estimates for each time step
#because of weird computer stuff, it spits out info multiple time per time-step, that's normal

#example use case: with the default parameters + density model with nothing linked, can you find 
# the time when number of infections is magically HIGHER than number of bacteria? ;)

#(remember, ignore anything <1e-15, it's just noise)


#plot the results!
#don't worry about the warnings!
ggplot(results) +
  geom_line(aes(time, B, colour = "B")) +
  geom_line(aes(time, P, colour = "P")) +
  scale_y_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x),
                     labels=trans_format("log10", math_format(10^.x))) +
  coord_cartesian(ylim = c(0.1, 1e11)) +
  labs(x = "Time (hours)", y = "cfu or pfu per mL", colour = "Organism:") +
  theme_bw()
