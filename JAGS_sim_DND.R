sink("dnd_sim.txt")
cat("
    model {
    
    # Priors
    int.psi ~ dunif(0, 1)               # Occupancy intercept
    int.p ~ dunif(0, 1)                 # Detection intercept

    # Likelihood
    # Ecological model for true occurrence
    for (i in 1:M) {
      z[i] ~ dbern(psi[i])
      logit(psi[i]) <- logit(int.psi) 
    }

    # Observation model for replicated detection/nondetection observations
    for (i in 1:M) {
      for (j in 1:n_surveys){ 
        logit(p[i,j]) <-  logit(int.p) 
        mup[i,j] <- z[i] * p[i,j]
        dnd[i,j] ~ dbern(mup[i,j])
      }
    }
  }
    ",fill = TRUE)
sink()
