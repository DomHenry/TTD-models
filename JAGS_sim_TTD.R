# ====================================================================================
#
#   JAGS_sim_TTD.R
#
#   JAGS formulation for hierarchical multi-survey time to detection occupancy model used 
#   in simulations (sim_runs.R)
#  
#   Henry et al. 2019. 
#
# ====================================================================================

sink("ttd_multi_sim.txt")
cat("
    model {
    
    # Priors
    int.psi ~ dunif(0, 1)               # Intercept occupancy on prob. scale
    int.lambda ~ dgamma(0.0001, 0.0001) # Poisson rate parameter

    # Likelihood
    # Ecological model for true occurrence
    for (i in 1:M) {
      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- logit(int.psi) 
    }

    # Observation model
    # Exponential model for time to detection ignoring censoring
    for (i in 1:M) {
      for (j in 1:n_surveys){ 
        ttd[i,j] ~ dexp(lambda[i,j])
        log(lambda[i,j]) <- log(int.lambda)
        # Model for censoring due to species absence and ttd>=Tmax
        d[i,j] ~ dbern(theta[i,j])
        theta[i,j] <- z[i] * step(ttd[i,j] - Tmax) + (1 - z[i])
      }
    }
  }
    ",fill = TRUE)
sink()
