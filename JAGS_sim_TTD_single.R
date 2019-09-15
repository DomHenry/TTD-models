# ====================================================================================
#
#   JAGS_sim_TTD.R
#
#   JAGS formulation for hierarchical single-survey time to detection occupancy model used 
#   in simulations (sim_runs.R)
#  
#   Henry et al. 2019. 
#
# ====================================================================================

sink("ttd_single_sim.txt")
cat("
    model {

    # Priors
    int.psi ~ dunif(0, 1)               # Intercept occupancy on prob. scale
    int.lambda ~ dgamma(0.0001, 0.0001) # Poisson rate parameter

    # Likelihood
    for (i in 1:M){
      
      # Model for occurrence
      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- logit(int.psi) #+ beta1 * covB[i]
        
      # Observation model
      # Exponential model for time to detection ignoring censoring
      ttd[i] ~ dexp(lambda[i])
      log(lambda[i]) <- log(int.lambda) 
      
      # Model for censoring due to species absence and ttd>=Tmax
      d[i] ~ dbern(theta[i])
      theta[i] <- z[i] * step(ttd[i] - Tmax) + (1 - z[i])
  }
}
    ",fill = TRUE)
sink()



