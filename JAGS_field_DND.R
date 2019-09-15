# ====================================================================================
#
#   JAGS_field_DND.R
#
#   JAGS formulation for hierarchical multi-survey detection/non-detection occupancy 
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

sink("DNDmod.txt")
cat("
    model {
    
    # Priors for species-specific effects 
    for(k in 1:n_spp){                    # Loop over species
      lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
      lp[k] ~ dnorm(mu.lp, tau.lp)
      beta1[k] ~ dnorm(mu.beta1, tau.beta1)                # NDVI
      beta2[k] ~ dnorm(mu.beta2, tau.beta2)                # Precipitation Concentration Index
      beta3[k] ~ dnorm(mu.beta3, tau.beta3)                # Elevation   
      beta4[k] ~ dnorm(mu.beta4, tau.beta4)                # Terrain Ruggedness Index 
    }
    
    # Hyperpriors for model of occupancy #
    # ********************************** #
    mu.psi ~ dunif(0,1)
    mu.lpsi <- logit(mu.psi)
    sd.lpsi ~ dunif(0,5)
    tau.lpsi <- pow(sd.lpsi, -2)
    
    mu.p ~ dunif(0,1)
    mu.lp <- logit(mu.p)
    sd.lp ~ dunif(0,5)
    tau.lp <- pow(sd.lp, -2)
    
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
    
    # Detection covariate priors #
    # ************************** #
    
    alpha1 ~ dnorm(0, 0.1)       # Temperature
    alpha2 ~ dnorm(0, 0.1)       # Wind speed
    alpha3 ~ dnorm(0, 0.1)       # Cloud
    alpha4 ~ dnorm(0, 0.1)       # Time of day (mins)
    alpha5 ~ dnorm(0, 0.1)       # mins (quad)
    alpha6 ~ dnorm(0, 0.1)       # Julian day  
    
    # Ecological model for latent occurrence z #
    # **************************************** #
    
    for(k in 1:n_spp){                # Loop over species
      for (i in 1:n_pen) {            # Loop over sites
        z[i,k] ~ dbern(psi[i,k])
        logit(psi[i,k]) <- lpsi[k] + beta1[k]*ndvi[i] + beta2[k]*map_ctn[i] + 
                                     beta3[k]*elev[i] + beta4[k]*tri[i] 
      }
    }
    
    
    # Detection/non-detection observation model #
    # ***************************************** #
    for(k in 1:n_spp){ 
      for (i in 1:n_pen) {
        for (j in 1:rep_vec[i]){ 
          logit(p[i,j,k]) <-lp[k] +                                            
                            alpha1*temp[i,j] +  alpha2*wind[i,j] +
                            alpha3*cloud[i,j] + alpha4*mins[i,j] + 
                            alpha5*pow(mins[i,j],2) + alpha6*jday[i,j]
          
          mup[i,j,k] <- z[i,k] * p[i,j,k]
          Y[i,j,k] ~ dbern(mup[i,j,k])
        }
      }
    }
    
    # Derived quantities #
    # ****************** #
    for(k in 1:n_spp){            # Loop over species
      Nocc.fs[k] <- sum(z[,k])    # Add up number of occupied sites among the list of species
      }
    
    for (i in 1:n_pen) {          # Loop over sites
      Nsite[i] <- sum(z[i,])      # Add up number of occurring species at each site
      }
    }
    ",fill = TRUE)
sink()