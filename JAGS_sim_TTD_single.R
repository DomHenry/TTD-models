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
    int.psi ~ dunif(0, 1)               # Occupancy intecept
    int.lambda ~ dgamma(0.0001, 0.0001) # Detection rate intercept

    # Likelihood
    for (i in 1:M){
      
      # Ecological model for true occurrence
      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- logit(int.psi) 
      
      # Exponential model for time to detection
      ttd[i] ~ dexp(lambda[i])
      log(lambda[i]) <- log(int.lambda) 
      
      # Model for censoring due to species absence and ttd >= Tmax
      d[i] ~ dbern(theta[i])
      theta[i] <- z[i] * step(ttd[i] - Tmax) + (1 - z[i])
  }
}
    ",fill = TRUE)
sink()



