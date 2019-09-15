# ====================================================================================
#
#   sim_results_plots.R
#
#   Code used to plot results of simulations (Figs. 2, 3, S1, S2, 
#   S3 and S4 in manuscript)
#  
#   Henry et al. 2019. 
#
# ====================================================================================

# install.packages("tidyverse")
library(tidyverse)

# Load data ---------------------------------------------------------------
source("array_to_df.R")
load("simulation_results.RData")

# Create occupancy plotting data frames -----------------------------------
combos <- list(combo_list)
measures <- list("MSE","post_mean","post_sd")
int <- list("int.psi")

ttd1_psi <- pmap(.l = list(combos,list(occ_sim_results1),int,list("ttd1_psi_mse","ttd1_psi_mean","ttd1_psi_sd"),measures),
                 .f = arr2df) 

ttd2_psi <- pmap(.l = list(combos,list(arr_list_ttd[[1]]),int,list("ttd2_psi_mse","ttd2_psi_mean","ttd2_psi_sd"),measures),
                 .f = arr2df) 

ttd4_psi <- pmap(.l = list(combos,list(arr_list_ttd[[2]]),int,list("ttd4_psi_mse","ttd4_psi_mean","ttd4_psi_sd"),measures),
                 .f = arr2df) 

ttd8_psi <- pmap(.l = list(combos,list(arr_list_ttd[[3]]),int,list("ttd8_psi_mse","ttd8_psi_mean","ttd8_psi_sd"),measures),
                 .f = arr2df) 

dnd2_psi <- pmap(.l = list(combos,list(arr_list_dnd[[1]]),int,list("dnd2_psi_mse","dnd2_psi_mean","dnd2_psi_sd"),measures),
                 .f = arr2df) 

dnd4_psi <- pmap(.l = list(combos,list(arr_list_dnd[[2]]),int,list("dnd4_psi_mse","dnd4_psi_mean","dnd4_psi_sd"),measures),
                 .f = arr2df) 

dnd8_psi <- pmap(.l = list(combos,list(arr_list_dnd[[3]]),int,list("dnd8_psi_mse","dnd8_psi_mean","dnd8_psi_sd"),measures),
                 .f = arr2df) 

plotdf_psi <- bind_cols(ttd1_psi,ttd2_psi,ttd4_psi,ttd8_psi,
                        dnd2_psi,dnd4_psi,dnd8_psi) %>% 
  select(-contains("com")) %>% 
  mutate(combo = ttd1_psi[[1]]$combo) %>% 
  gather(protocol,value,-combo) %>% 
  separate(protocol, c("protocol","param","measure"))

# Create detection plotting data frames -----------------------------------
int <- list("int.p")

ttd1_p <- pmap(.l = list(combos,list(occ_sim_results1),int,list("ttd1_p_mse","ttd1_p_mean","ttd1_p_sd"),measures),
               .f = arr2df) 

ttd2_p <- pmap(.l = list(combos,list(arr_list_ttd[[1]]),int,list("ttd2_p_mse","ttd2_p_mean","ttd2_p_sd"),measures),
               .f = arr2df) 

ttd4_p <- pmap(.l = list(combos,list(arr_list_ttd[[2]]),int,list("ttd4_p_mse","ttd4_p_mean","ttd4_p_sd"),measures),
               .f = arr2df) 

ttd8_p <- pmap(.l = list(combos,list(arr_list_ttd[[3]]),int,list("ttd8_p_mse","ttd8_p_mean","ttd8_p_sd"),measures),
               .f = arr2df) 

dnd2_p <- pmap(.l = list(combos,list(arr_list_dnd[[1]]),int,list("dnd2_p_mse","dnd2_p_mean","dnd2_p_sd"),measures),
               .f = arr2df) 

dnd4_p <- pmap(.l = list(combos,list(arr_list_dnd[[2]]),int,list("dnd4_p_mse","dnd4_p_mean","dnd4_p_sd"),measures),
               .f = arr2df) 

dnd8_p <- pmap(.l = list(combos,list(arr_list_dnd[[3]]),int,list("dnd8_p_mse","dnd8_p_mean","dnd8_p_sd"),measures),
               .f = arr2df) 

plotdf_p <- bind_cols(ttd1_p,ttd2_p,ttd4_p,ttd8_p,
                      dnd2_p,dnd4_p,dnd8_p) %>% 
  select(-contains("com")) %>% 
  mutate(combo = ttd1_p[[1]]$combo) %>% 
  gather(protocol,value,-combo) %>% 
  separate(protocol, c("protocol","param","measure"))

plotdf <- bind_rows(plotdf_psi,plotdf_p)

labwidth <- 20

combo_names <- list(
  'psi1.p1' = str_wrap("Spp A: rare & very inconspicuous",width = labwidth),
  'psi1.p2' = str_wrap("Spp B: rare & inconspicuous",width = labwidth),
  'psi1.p3' = str_wrap("Spp C: rare & conspicuous", width = labwidth),
  'psi2.p1' = str_wrap("Spp D: widespread & very inconspicuous",width = labwidth),
  'psi2.p2' = str_wrap("Spp E: widespread & inconspicuous",width = labwidth),
  'psi2.p3' = str_wrap("Spp F: widespread & conspicuous", width = labwidth)
)

combo_labeller <- function(variable,value){
  return(combo_names[value])
}


# Create geometery for the shaded panels
shading <- tibble(min = seq(from = 0.5, 
                            to = max(as.numeric(as.factor(plotdf$protocol))), 
                            by = 1),
                  max = seq(from = 1.5, to = max(as.numeric(as.factor(plotdf$protocol))) + 0.5, 
                            by = 1),
                  colshade = as_factor(c(0,0,1,1,0,0,0)),
                  protocol = c("ttd1", "dnd2", "ttd2","dnd4", "ttd4","dnd8", "ttd8"))

plotdf <- plotdf %>% 
  left_join(shading, by = "protocol") %>% 
  mutate(combo = as.factor(combo)) %>%
  mutate(protocol = fct_relevel(protocol, 
                                "ttd1", "dnd2", "ttd2","dnd4", "ttd4","dnd8", "ttd8"))

plotsims1 <- function(meas, para, ylims, label){
  
  p <- plotdf %>% 
    filter(measure == meas & param == para) %>% 
    ggplot(aes(y=value,x=protocol))+
    labs(title = label, x = "", y = "")+
    geom_blank()+
    geom_rect(aes(xmin = min, xmax = max,
                  ymin = -Inf, ymax = Inf,
                  fill = colshade),
              alpha = 0.009)+
    geom_boxplot(aes(y=value,
                     x=protocol))+
    scale_fill_manual(values = c("gray63","white"),
                      labels = element_blank())+
    theme_classic()+
    theme(legend.title = element_blank(),
          axis.text = element_text(size = 13),
          strip.text.x =  element_text(size = 13),
          legend.text = element_text(size = 16),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 18))+
    scale_y_continuous(limits = ylims)+
    facet_wrap(~combo, labeller = combo_labeller,scales = "free")
  return(p)
}

# Fig. 2
plotsims1(meas = "mse", para = "p", ylims = c(0,0.2),label = "Mean square error (detection)")

# Fig. 3
plotsims1(meas = "mse", para = "psi", ylims = c(0,0.2),label = "Mean square error (occupancy)")

# Fig. S1
plotsims1(meas = "sd", para = "p", ylims = c(0,0.4),label = "Parameter standard deviation (detection)")

# Fig. S2
plotsims1(meas = "sd", para = "psi", ylims = c(0,0.4),label = "Parameter standard deviation (occupancy)")



plotsims2 <- function(para, label){
  
  p <- plotdf %>% 
    filter(measure == "mean" & param == para) %>% 
    ggplot(aes(y=value,x=protocol))+
    labs(title = label, x = "", y = "")+
    geom_blank()+
    geom_rect(aes(xmin = min, xmax = max,
                  ymin = -Inf, ymax = Inf,
                  fill = colshade),
              alpha = 0.009)+
    geom_boxplot(aes(y=value,
                     x=protocol))+
    scale_fill_manual(values = c("gray63","white"),
                      labels = element_blank())+
    theme_classic()+
    theme(legend.title = element_blank(),
          axis.text = element_text(size = 13),
          strip.text.x =  element_text(size = 13),
          legend.text = element_text(size = 16),
          panel.border = element_rect(color = "black", fill = NA),
          plot.title = element_text(size = 18))+
    scale_y_continuous(breaks = c(0,0.3,0.6,0.9), limits = c(0,1))
  
  if(para == "psi") {
    
    p <- p +
      geom_hline(data=filter(plotdf_p, combo=="psi1.p1"), aes(yintercept=0.3), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi1.p2"), aes(yintercept=0.3), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi1.p3"), aes(yintercept=0.3), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi2.p1"), aes(yintercept=0.6), linetype = "dashed", col = "red", size = 1.0) +
      geom_hline(data=filter(plotdf_p, combo=="psi2.p2"), aes(yintercept=0.6), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi2.p3"), aes(yintercept=0.6), linetype = "dashed", col = "red", size = 1.0) + 
      facet_wrap(~combo, labeller = combo_labeller,scales = "free")
  } else {
    p <- p +
      geom_hline(data=filter(plotdf_p, combo=="psi1.p1"), aes(yintercept=0.1), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi1.p2"), aes(yintercept=0.3), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi1.p3"), aes(yintercept=0.6), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi2.p1"), aes(yintercept=0.1), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi2.p2"), aes(yintercept=0.3), linetype = "dashed", col = "red", size = 1.0) + 
      geom_hline(data=filter(plotdf_p, combo=="psi2.p3"), aes(yintercept=0.6), linetype = "dashed", col = "red", size = 1.0) + 
      facet_wrap(~combo, labeller = combo_labeller,scales = "free")
  }
  return(p)  
}

# Fig. S3
plotsims2(para = "p", label = "Parameter mean (detection)")

# Fig. S4
plotsims2(para = "psi", label = "Parameter mean (occupancy)")
