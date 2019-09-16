# ====================================================================================
#
#   array_to_df.R
#
#   Function to convert simulation output arrays to data frames 
#   for plotting in sim_result_plots.R                 
#  
#   Henry et al. 2019. 
#
# ====================================================================================

arr2df <- function(combo_list,simarray,int,coln,measure) {
  
  dflist <- list()
  
  for(i in 1:length(combo_list)){
    
    ref <- unlist(strsplit(combo_list[i],"[.]"))
    s <- as.numeric(gsub("\\D", "", ref))[1] 
    p <- as.numeric(gsub("\\D", "", ref))[2]
    
    df <- data.frame(name = as.vector(simarray[int,measure,,s,p])) 
    colnames(df) <- combo_list[[i]]
    dflist[[i]] <- df
    
  }  
  
  dfall <- do.call(cbind, dflist) 
  plotdf <- dfall %>% gather(combo,!!coln) %>% as_tibble()
  return(plotdf)
}