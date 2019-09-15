# ====================================================================================
#
#   sim_runs.R
#
#   Code used to run simulations comparing time to detection and 
#   detection/non-detection occupancy models under different survey replicates
#   used in simulations (sim_runs.R)
#  
#   Henry et al. 2019. 
#
# ====================================================================================

# install.packages("AHMbook")
# install.packages("tidyverse")
# install.packages("jagsUI")

library(AHMbook)
library(jagsUI)
source("sim_ttd_data.R")


# Define initial parameters -----------------------------------------------
nsims <- 1000   # Number of simulations
nvars <- 3      # Number of variables in arrays (posterior mean, posterior SD and MSE)
psi_levels <- 2 # Number of simulated occupancy levels (psi = 0.3 & 0.6)
p_levels <- 3   # Number of simulated detection levels (p = 0.1, 0.3 & 0.6)
ints <- 2       # Number of intercepts (int.psi and int.p)

# Define dimension names --------------------------------------------------
dim_n1 <- c("int.psi","int.p")
dim_n2 <- c("post_mean","post_sd","MSE")
dim_n3 <- NULL
dim_n4 <- c("psi1","psi2")
dim_n5 <- c("p1","p2", "p3")

# Create output arrays ----------------------------------------------------
occ_sim_results1 <-   # 1. TTD1
  occ_sim_results2 <- # 2. TTD2
  occ_sim_results3 <- # 3. TTD4
  occ_sim_results4 <- # 4. TTD8
  occ_sim_results5 <- # 5. DND2
  occ_sim_results6 <- # 6. DND4
  occ_sim_results7 <- # 7. DND8
  array(dim=c(ints,nvars,nsims,psi_levels,p_levels)) 

dimnames(occ_sim_results1) <- dimnames(occ_sim_results2) <- dimnames(occ_sim_results3) <- 
  dimnames(occ_sim_results4) <- dimnames(occ_sim_results5) <- dimnames(occ_sim_results6) <- 
  dimnames(occ_sim_results7) <- list(dim_n1,dim_n2,dim_n3,dim_n4,dim_n5)

# Combine to list of TTD arrays
arr_list_ttd <- list(occ_sim_results2,occ_sim_results3,occ_sim_results4)
## Note that TTD 1 is modelled and stored separately (not in the list above) ##

# Combine to list of DND arrays
arr_list_dnd <- list(occ_sim_results5,occ_sim_results6,occ_sim_results7)


# Create occupancy and detection combinations -----------------------------
psi_vec <- c(0.3,0.6)
p_vec <- c(0.1,0.3,0.6)

combos <- expand.grid(dim_n4,dim_n5)
combo_list <- paste(combos[,1],combos[,2],sep= ".")


# Assign number of survey replicates in each for each protocol ------------
n_survs_ttd <- c(2,4,8) 
n_surs_dnd <- c(2,4,8)

# Assign JAGS chain values --------------------------------------------------

# Full run
ni <- 80000 # iterations
nb <- 40000 # burn-in
nc <- 3     # chains
nt <- 20    # thinning rate

# Testing phase
ni <- 10000
nb <- 3000
nc <- 3
nt <- 20

# Start simulations -------------------------------------------------------
for(i in 1:nsims){
  for(s in seq_along(psi_vec)){
    for(p in seq_along(p_vec)){
      cat("\n sim number", i)
      cat("\n psi level",psi_vec[s])
      cat("\n p level ", p_vec[p], "\n")
      
      ## TIME TO DETECTION MODELS ##
      ## -------------------------##
      
      # Simulate dataset for 8 TTD surveys
      sim_data <- sim_ttd_data(M = 50, mean.psi = psi_vec[s],mean.lambda = p_vec[p], Tmax = 10, n_surveys = 8)
      
      # Create subsets of sample data without replacement
      ran8 <- c(1:8)         # TTD8
      ran4 <- sample(ran8,4) # TTD4
      ran2 <- sample(ran4,2) # TTD2
      ran1 <- sample(ran2,1) # TTD1
      
      # Single survey simulation model
      source("JAGS_sim_TTD_single.R")
      sim_datasub <- sim_data
      sim_datasub$ttd_uncen <- sim_datasub$ttd_uncen[,ran1]
      sim_datasub$ttd <- sim_datasub$ttd[,ran1]
      sim_datasub$d <- sim_datasub$d[,ran1]
      sim_datasub$n_surveys <- length(ran1) 
      
      zst <- rep(1, length(sim_datasub$ttd))
      ttdst <-rep(sim_datasub$Tmax+1, sim_datasub$M)
      ttdst[sim_datasub$d == 0] <- NA
      inits <- function(){list(z =zst, ttd = ttdst, int.psi = runif(1), int.lambda = runif(1))}
      
      params <- c("int.psi", "int.lambda", "n.occ")
      out_ttd_single <- jags(sim_datasub, inits, params, "ttd_single_sim.txt", n.chains=nc, n.iter=ni,
                             n.burn = nb, n.thin=nt, parallel = T, verbose = FALSE)
      
      occ_sim_results1[1,1,i,s,p] <- out_ttd_single$summary[1,1] # psi posterior mean
      occ_sim_results1[2,1,i,s,p] <- out_ttd_single$summary[2,1] # p posterior mean
      occ_sim_results1[1,2,i,s,p] <- out_ttd_single$summary[1,2] # psi posterior SD
      occ_sim_results1[2,2,i,s,p] <- out_ttd_single$summary[2,2] # p posterior SD
      
      # MSE = SD^2 + bias^2
      occ_sim_results1[1,3,i,s,p] <- out_ttd_single$summary[1,2]^2 + (out_ttd_single$summary[1,1] - psi_vec[s])^ 2  # psi MSE
      occ_sim_results1[2,3,i,s,p] <- out_ttd_single$summary[2,2]^2 + (out_ttd_single$summary[2,1] - p_vec[p])^2     # p MSE
      
      # Multiple TTD survey simulations
      surv_list <- list(ran2,ran4,ran8)
      
      for(n in 1:length(surv_list)){
        
        sim_datasub <- sim_data
        sim_datasub$ttd_uncen <- sim_datasub$ttd_uncen[,surv_list[[n]]]
        sim_datasub$ttd <- sim_datasub$ttd[,surv_list[[n]]]
        sim_datasub$d <- sim_datasub$d[,surv_list[[n]]]
        sim_datasub$n_surveys <- length(surv_list[[n]]) 
        
        source("JAGS_sim_TTD.R")
        zst <- rep(1, length(sim_datasub$ttd[,1]))
        ttdst <- sim_datasub$ttd
        I <- which(is.na(sim_datasub$ttd))
        ttdst[I] <- sim_datasub$Tmax
        ttdst[-I] <- NA
        inits <- function(){list(z =zst, ttd = ttdst, int.psi = runif(1), int.lambda = runif(1))}
        
        params <- c("int.psi", "int.lambda", "n.occ")
        out_ttd_multi <- jags(sim_datasub, inits, params, "ttd_multi_sim.txt", n.chains=nc, n.iter=ni, 
                              n.burn = nb, n.thin=nt, parallel = T, verbose = FALSE)
        
        arr_list_ttd[[n]][1,1,i,s,p] <- out_ttd_multi$summary[1,1]
        arr_list_ttd[[n]][2,1,i,s,p] <- out_ttd_multi$summary[2,1]
        arr_list_ttd[[n]][1,2,i,s,p] <- out_ttd_multi$summary[1,2]
        arr_list_ttd[[n]][2,2,i,s,p] <- out_ttd_multi$summary[2,2]
        
        # MSE = SD^2 + bias^2
        arr_list_ttd[[n]][1,3,i,s,p] <-  out_ttd_multi$summary[1,2]^2 + (out_ttd_multi$summary[1,1] - psi_vec[s])^2
        arr_list_ttd[[n]][2,3,i,s,p] <- out_ttd_multi$summary[2,2]^2 + (out_ttd_multi$summary[2,1] - p_vec[p])^2
        
      }
      
      ## DETECTION/NON-DETECTION MODELS ##
      ## ------------------------------ ##
      
      # Simulate dataset for 8 DND surveys
      sim_data <- AHMbook::simOcc(M = 50, J = 8, mean.occupancy = psi_vec[s], mean.detection = p_vec[p],
                                  beta1 = 0, beta2 = 0, beta3 = 0,time.effects = c(0, 0), # Supress covariate effects
                                  alpha1 = 0, alpha2 = 0, alpha3 = 0, sd.lp = 0,          # Supress covariate effects
                                  b = 0, show.plot = F)
      
      # Need to create subsets of sample data without replacement
      ran8 <- c(1:8)         # DND8
      ran4 <- sample(ran8,4) # DND4
      ran2 <- sample(ran4,2) # DND2
      
      # Multiple TTD survey simulations
      surv_list <- list(ran2,ran4,ran8)
      
      for(n in 1:length(surv_list)){
        
        sim_datasub <- list(M = sim_data$M, n_surveys = length(surv_list[[n]]), dnd = sim_data$y[,surv_list[[n]]])  
        
        source("JAGS_sim_DND.R")
        zst <- rep(1, length(sim_datasub$dnd[,1]))
        zst[is.na(zst)] <- 1
        inits <- function() list(z = zst)
        
        params <- c("int.psi", "int.p", "n.occ")
        out_dnd <- jags(sim_datasub, inits, params, "dnd_sim.txt", n.chains=nc, n.iter=ni, 
                        n.burn = nb, n.thin=nt, parallel = T,verbose = FALSE)
        
        arr_list_dnd[[n]][1,1,i,s,p] <- out_dnd$summary[1,1]
        arr_list_dnd[[n]][2,1,i,s,p] <- out_dnd$summary[2,1]
        arr_list_dnd[[n]][1,2,i,s,p] <- out_dnd$summary[1,2]
        arr_list_dnd[[n]][2,2,i,s,p] <- out_dnd$summary[2,2]
       
        # MSE = SD^2 + bias^2 #
        arr_list_dnd[[n]][1,3,i,s,p] <- out_dnd$summary[1,2]^2 + (out_dnd$summary[1,1] - psi_vec[s])^2  
        arr_list_dnd[[n]][2,3,i,s,p] <- out_dnd$summary[2,2]^2 + (out_dnd$summary[2,1] - p_vec[p])^2     
        
      }
    }
  }
}


save.image("simulation_results_test.RData") # Save data for plotting
