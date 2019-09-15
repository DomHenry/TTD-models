# ====================================================================================
#
#   sim_ttd_data.R
#
#   Function to generate time to detection data used in simulations (sim_runs.R)
#  
#   Henry et al. 2019. 
#
#   Note: This occupancy data-generation function is based on the 
#         AHMbook::simOccttd function.
#
#   Citation: Marc Kery, Andy Royle and Mike Meredith (2017). 
#             AHMbook: Functions and Data for the Book 'Applied Hierarchical 
#             Modeling in Ecology'. R package version 0.1.4. 
#             https://CRAN.R-project.org/package=AHMbook)
#
# ====================================================================================

sim_ttd_data <- function(M, mean.psi, mean.lambda, Tmax, n_surveys) {
  
  z_vec <- rbinom(M, 1, mean.psi)
  lambda <- rep(mean.lambda,M)
  
  # Create uncensored TTD matrix
  ttd_uncen <- matrix(nrow = M, ncol= n_surveys) 
  for(i in 1:n_surveys){
    ttd_uncen[,i] <- rexp(M, lambda)
  }
  
  # Create censored TTD matrix
  ttd <- ttd_uncen
  ttd[z_vec == 0,] <- NA            # Censor on occurence
  ttd[ttd_uncen >= Tmax] <- NA      # Censor on max survey time
  
  # Create indicator matrix
  d <- matrix(as.numeric(is.na(ttd)), nrow=dim(ttd)[1])
  
  return(list(M = M, mean.psi = mean.psi, mean.lambda = mean.lambda, 
              Tmax = Tmax, z_vec = z_vec, n_surveys = n_surveys,
              ttd_uncen = ttd_uncen, 
              ttd = ttd, d = d))
}
