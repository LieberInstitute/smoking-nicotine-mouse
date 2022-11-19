
## 1.2 Comparison of DEG 

load(here("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
     


### 1.2.1 T-stats plots

## Function to add DE info of exons in both groups of samples
add_DE_info <-function(top_exons1, top_exons2, name_1, name_2) {
  DE<-vector()
  for (i in 1:dim(top_exons1)[1]) {
    ## DE exons in both groups
    if (top_exons1$adj.P.Val[i]<0.05 && top_exons2$adj.P.Val[i]<0.05) {
      DE<-append(DE, "sig Both")
    }
    ## DE exons in only one of the groups
    else if (top_exons1$adj.P.Val[i]<0.05 && !top_exons2$adj.P.Val[i]<0.05) {
      DE<-append(DE, paste("sig", name_1))
    }
    
    else if (top_exons2$adj.P.Val[i]<0.05 && !top_exons1$adj.P.Val[i]<0.05) {
      DE<-append(DE, paste("sig",name_2))
    }
    ## No DE genes in neither group
    else {
      DE<-append(DE, "None")
    }
  }
  return(DE)
}



## Compare t-stats of exons from different groups of samples
t_stat_plot <- function(top_exons1, top_exons2, name_1, name_2, title){
  
  ## Correlation coeff
  rho <- cor(top_exons1$t, top_exons2$t, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Merge data
  t_stats<-data.frame(t1=top_exons1$t, t2=top_exons2$t)
  ## Add DE info for both groups
  t_stats$DE<-add_DE_info(top_exons1, top_exons2, name_1, name_2)
  
  cols <- c("red", "#ffad73","#26b3ff", "dark grey") 
  names(cols)<-c("sig Both", paste0("sig ", name_1), paste0("sig ", name_2), "None")
  alphas <- c( 1, 1, 1,0.5)  
  names(alphas)<-c("sig Both", paste0("sig ", name_1), paste0("sig ", name_2), "None")
  
  plot <- ggplot(t_stats, aes(x = t1, y = t2, color=DE, alpha=DE)) +
    geom_point(size = 1) +
    labs(x = paste("t-stats", name_1), 
         y = paste("t-stats", name_2),
         title = title, 
         subtitle = rho_anno, 
         parse = T) +
    theme_bw() +
    scale_color_manual(values = cols) + 
    scale_alpha_manual(values = alphas)
  
  plot
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Exon_analysis/t_stats_", gsub(" ", "_", title), 
                        ".pdf", sep=""), height = 20, width = 25, units = "cm")
}



## Function to create multiple t-stats plots
tstats_plots<-function(top_exons_pairs, name_1, name_2, titles){
  
  plots<-list()
  for (i in 1:length(top_exons_pairs)){
    p<-t_stat_plot(top_exons_pairs[[i]][[1]], top_genes_pairs[[i]][[2]], name_1, name_2, titles[i])
    plots[[i]]<-p
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=2)
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/Exon_analysis/t_stats_",gsub(" ", "_", name_1), "_VS_",
                        gsub(" ", "_", name_2), ".pdf", sep=""), 
         height = 20, width = 25, units = "cm")
}





################################
# Smoking vs nicotine (exons)
################################
t_stat_plot(top_exons_nic, top_exons_smo, "Nicotine pups", "Smoking pups", "Nic vs Smo DE exons")



##########################################
# Smoking genes vs smoking exons' genes 
##########################################
t_stat_plot(top_exons_nic, top_exons_smo, "Nicotine pups", "Smoking pups", "Nic vs Smo DE exons")






top_exons_pairs<-list(list(top_genes_blood_smoking_naive, top_genes_adults_smoking_naive), 
                      list(top_genes_blood_smoking_fitted, top_genes_adults_smoking_fitted),
                      list(top_genes_blood_smoking_interaction, top_genes_adults_smoking_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking blood", "Smoking brain", models)


