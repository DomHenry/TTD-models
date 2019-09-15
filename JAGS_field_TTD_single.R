# ====================================================================================
#
#   JAGS_field_TTD_single.R
#
#   JAGS formulation for hierarchical single-survey time to detection occupancy 
#   model used to analyse bird community data from South Africa.
#  
#   Henry et al. 2019. 
#
#   Notes:  These models are based on those developed in Kéry & Royle (2016) Applied 
#           Hierarchical Modeling in Ecology - Chapter 10 and by the models in the 
#           supplementary material of Bornand, C. N., Kéry, M., Bueche, L., & 
#           Fischer, M. (2014). Hide-and-seek in vegetation: Time-to-detection is 
#           an efficient design for estimating detectability and 
#           occurrence. Methods in Ecology and Evolution, 5(5), 433–442.
#
# ====================================================================================

sink("TTDmod_single.txt")
cat("
    model {
    
    # Priors for species-specific effects 
    for(k in 1:n_spp){	                                 
      logLambda[k] ~ dnorm(mu.logLambda, tau.logLambda) 
      lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)                     
      beta1[k] ~ dnorm(mu.beta1, tau.beta1)                # NDVI
      beta2[k] ~ dnorm(mu.beta2, tau.beta2)                # Precipitation Concentration Index
      beta3[k] ~ dnorm(mu.beta3, tau.beta3)                # Elevation   
      beta4[k] ~ dnorm(mu.beta4, tau.beta4)                # Terrain Ruggedness Index
    }	
    
    # Hyperpriors for model of occupancy #
    # ********************************** #    
    mu.psi <- ilogit(mu.lpsi)
    mu.lpsi ~ dt(muT,tauT,kT)    
    
    muT <- 0
    tauT <- pow(sdT, -2)
    sdT <- 1.566
    kT <- 7.763
    
    # The Cauchy distribution is a special case of 
    # the t distribution, with 1 degree of freedom. The half-Cauchy prior for 
    # standard deviation (not variance) can be coded by truncating the t distribution:
    # dt(mu, tau, 1)T(0,) - truncating the positve values
    
    sd.lpsi ~ dt(0,1,1)T(0,)  # Half-Cauchy distribution
    tau.lpsi <- pow(sd.lpsi, -2)
    
    mu.beta1 ~ dnorm(0, 0.1)
    sd.beta1 ~ dt(0,1,1)T(0,)  
    tau.beta1 <- pow(sd.beta1, -2)
    
    mu.beta2 ~ dnorm(0, 0.1)
    sd.beta2 ~ dt(0,1,1)T(0,)
    tau.beta2 <- pow(sd.beta2, -2)
    
    mu.beta3 ~ dnorm(0, 0.1)
    sd.beta3 ~ dt(0,1,1)T(0,)
    tau.beta3 <- pow(sd.beta3, -2)
    
    mu.beta4 ~ dnorm(0, 0.1)
    sd.beta4 ~ dt(0,1,1)T(0,)
    tau.beta4 <- pow(sd.beta4, -2)
    
    # mu.beta5 ~ dnorm(0, 0.1)
    # sd.beta5 ~ dt(0,1,1)T(0,)
    # tau.beta5 <- pow(sd.beta5, -2)
    
    # Hyperpriors for model of detection #
    # ********************************** #
    
    mu.logLambda <- log(lambdaP)
    lambdaP ~ dgamma(0.0001, 0.0001)
    
    sd.logLambda ~ dt(0,1,1)T(0,) 
    tau.logLambda <- pow(sd.logLambda, -2)
    
    # Detection covariate priors #
    # ************************** #
    
    gamma1 ~ dnorm(0, 0.1)       # Temperature
    gamma2 ~ dnorm(0, 0.1)       # Wind speed
    gamma3 ~ dnorm(0, 0.1)       # Cloud
    gamma4 ~ dnorm(0, 0.1)       # Time of day (mins)
    gamma5 ~ dnorm(0, 0.1)       # mins (quad)
    gamma6 ~ dnorm(0, 0.1)       # Julian day
    
    # Ecological model for latent occurrence z #
    # **************************************** #
    
    for(k in 1:n_spp){ 
      for (i in 1:n_pen) {
        ## Ecological model of true occurence [Occ z at site i for spp k]
        z[i,k] ~ dbern(psi[i,k])                           
        logit(psi[i,k]) <- lpsi[k] + beta1[k]*ndvi[i] + beta2[k]*map_ctn[i] + 
                                     beta3[k]*elev[i] + beta4[k]*tri[i]
        }
    }
    
    # Time-to-detection observation model #
    # *********************************** #
    for(k in 1:n_spp){
    for (i in 1:n_pen) {
  
      Y[i,k] ~ dexp(lambda[i,k])
      
        ## Linear model for log(rate)
        log(lambda[i,k]) <- logLambda[k] + gamma1*temp[i,k] + gamma2*wind[i,k] +
                              gamma3*cloud[i,k] + gamma4*mins[i,k] + 
                              gamma5*pow(mins[i,k],2) + gamma6*jday[i,k]
        
        ## Accomodation of z=0 and censoring
        d[i,k] ~ dbern(theta[i,k])
        theta[i,k] <- z[i,k] * step(Y[i,k] -  tmax) + (1 - z[i,k])
        
        }
    }
    
    
    # Derived quantities #
    # ****************** #
    
    for(k in 1:n_spp){
      Nocc.fs[k] <- sum(z[,k])	                  # Number of occupied sites
      lam.est.sp[k] <- exp(logLambda[k])
      p.sp[k]<- 1-exp(-lam.est.sp[k]*time)        # Mean species-specific detec prob after 700 seconds (T)
    }
    
    for (i in 1:n_pen) {                               
     Nsite[i] <- sum(z[i,])                      # Number of occurring species at each site
    } 
    
    # The Cauchy distribution is a special case of 
    # the t distribution, with 1 degree of freedom. The half-Cauchy prior for 
    # standard deviation (not variance) can be coded by truncating the t distribution:
    # dt(mu, tau, 1)T(0,) - truncating the positve values
    
    }
    ",fill=TRUE)
sink()
