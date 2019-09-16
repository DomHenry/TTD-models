# ====================================================================================
#
#   JAGS_field_TTD.R
#
#   JAGS formulation for a hierarchical multi-species, multi-survey time-to-detection 
#   occupancy model used to analyse bird community data from South Africa.
#  
#   Henry et al. 2019. 
#
#   Notes:  These models are based on those developed in Kéry & Royle (2016) Applied 
#           Hierarchical Modeling in Ecology - Chapter 10 and by the models in the 
#           Supporting Information of Bornand, C. N., Kéry, M., Bueche, L., & 
#           Fischer, M. (2014). Hide-and-seek in vegetation: Time-to-detection is 
#           an efficient design for estimating detectability and 
#           occurrence. Methods in Ecology and Evolution, 5(5), 433–442.
#
# ====================================================================================

sink("TTDmod.txt")
cat("
    model {
    
    # Priors for species-specific effects #
    # *********************************** #
    for(k in 1:n_spp){                                     # Loop over species 
      logLambda[k] ~ dnorm(mu.logLambda, tau.logLambda)    # Detection rate intercept
      lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)                   # Occupancy probability intercept  
      beta1[k] ~ dnorm(mu.beta1, tau.beta1)                # NDVI
      beta2[k] ~ dnorm(mu.beta2, tau.beta2)                # Precipitation Concentration Index
      beta3[k] ~ dnorm(mu.beta3, tau.beta3)                # Elevation   
      beta4[k] ~ dnorm(mu.beta4, tau.beta4)                # Terrain Ruggedness Index 
    }

    mu.psi <- ilogit(mu.lpsi)
    mu.lpsi ~ dt(muT,tauT,kT)    
    
    muT <- 0
    tauT <- pow(sdT, -2)
    sdT <- 1.566
    kT <- 7.763
    
    sd.lpsi ~ dt(0,1,1)T(0,)          # Half-Cauchy distribution
    tau.lpsi <- pow(sd.lpsi, -2)
    
    mu.logLambda <- log(lambdaP)
    lambdaP ~ dgamma(0.0001, 0.0001)
    
    sd.logLambda ~ dt(0,1,1)T(0,) 
    tau.logLambda <- pow(sd.logLambda, -2)
    
    # Occupancy covariate priors #
    # ************************** #
    
    mu.beta1 ~ dnorm(0, 0.1)          # NDVI
    sd.beta1 ~ dt(0,1,1)T(0,)  
    tau.beta1 <- pow(sd.beta1, -2)
    
    mu.beta2 ~ dnorm(0, 0.1)          # Precipitation Concentration Index
    sd.beta2 ~ dt(0,1,1)T(0,)
    tau.beta2 <- pow(sd.beta2, -2)
    
    mu.beta3 ~ dnorm(0, 0.1)          # Elevation
    sd.beta3 ~ dt(0,1,1)T(0,)
    tau.beta3 <- pow(sd.beta3, -2)
    
    mu.beta4 ~ dnorm(0, 0.1)          # Terrain Ruggedness Index 
    sd.beta4 ~ dt(0,1,1)T(0,)
    tau.beta4 <- pow(sd.beta4, -2)

    # Detection covariate priors #
    # ************************** #
    gamma1 ~ dnorm(0, 0.1)     # Temperature
    gamma2 ~ dnorm(0, 0.1)     # Wind speed
    gamma3 ~ dnorm(0, 0.1)     # Cloud
    gamma4 ~ dnorm(0, 0.1)     # Time of day (mins)
    gamma5 ~ dnorm(0, 0.1)     # mins (quad)
    gamma6 ~ dnorm(0, 0.1)     # Julian day
    
    # Ecological model for latent occurrence z #
    # **************************************** #
    
    for(k in 1:n_spp){                # Loop over species
      for (i in 1:n_pen) {            # Loop over sites
        z[i,k] ~ dbern(psi[i,k])                           
        logit(psi[i,k]) <- lpsi[k] + beta1[k]*ndvi[i] + beta2[k]*map_ctn[i] + 
                                     beta3[k]*elev[i] + beta4[k]*tri[i] 
      }
    }
    
    # Time-to-detection observation model #
    # *********************************** #
    for(k in 1:n_spp){                                    # Loop over species
      for (i in 1:n_pen) {                                # Loop over sites    
        for (j in 1:transvec[i]){                         # Loop over site survey replicates
        
          Y[i,j,k] ~ dexp(lambda[i,j,k])
          
          ## Linear model for log(detection rate)
          log(lambda[i,j,k]) <- logLambda[k] + gamma1*temp[i,j,k] + gamma2*wind[i,j,k] +
                                               gamma3*cloud[i,j,k] + gamma4*mins[i,j,k] + 
                                               gamma5*pow(mins[i,j,k],2) + gamma6*jday[i,j,k]
          
          ## Accomodation of z=0 and censoring
          d[i,j,k] ~ dbern(theta[i,j,k])
          theta[i,j,k] <- z[i,k] * step(Y[i,j,k] -  tmax) + (1 - z[i,k])
          
        }
      }
    }
    
    # Species richness #
    # **************** #
    for (i in 1:n_pen) {          # Loop over sites
      Nsite[i] <- sum(z[i,])      # Add up number of speices occurring at each site
    }
  }
    ",fill = TRUE)
sink()