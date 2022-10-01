
## 1.2 Comparison of DEG 

load(here("processed-data/04_DEA/top_genes_blood_smoking_naive.Rdata"))
load(here("processed-data/04_DEA/top_genes_blood_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_blood_smoking_interaction.Rdata"))

load(here("processed-data/04_DEA/top_genes_adults_nicotine_naive.Rdata"))
load(here("processed-data/04_DEA/top_genes_adults_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_adults_nicotine_interaction.Rdata"))

load(here("processed-data/04_DEA/top_genes_adults_smoking_naive.Rdata"))
load(here("processed-data/04_DEA/top_genes_adults_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_adults_smoking_interaction.Rdata"))

load(here("processed-data/04_DEA/de_genes_pups_nicotine_naive.Rdata"))
load(here("processed-data/04_DEA/top_genes_pups_nicotine_naive.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_pups_nicotine_interaction.Rdata"))

load(here("processed-data/04_DEA/de_genes_pups_smoking_naive.Rdata"))
load(here("processed-data/04_DEA/top_genes_pups_smoking_naive.Rdata"))
load(here("processed-data/04_DEA/de_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/04_DEA/top_genes_pups_smoking_interaction.Rdata"))



### 1.2.1 T-stats plots

## Function to add DE info of genes in both groups
add_DE_info <-function(top_genes1, top_genes2, name_1, name_2) {
  DE<-vector()
  for (i in 1:dim(top_genes1)[1]) {
    ## DE genes in both groups
    if (top_genes1$adj.P.Val[i]<0.05 && top_genes2$adj.P.Val[i]<0.05) {
      DE<-append(DE, "sig Both")
    }
    ## DE genes in only one of the groups
    else if (top_genes1$adj.P.Val[i]<0.05 && !top_genes2$adj.P.Val[i]<0.05) {
      DE<-append(DE, paste("sig", name_1))
    }
    
    else if (top_genes2$adj.P.Val[i]<0.05 && !top_genes1$adj.P.Val[i]<0.05) {
      DE<-append(DE, paste("sig",name_2))
    }
    ## No DE genes in neither group
    else {
      DE<-append(DE, "None")
    }
  }
  return(DE)
}


## Compare t-stats of genes from different samples' groups
t_stat_plot <- function(top_genes1, top_genes2, name_1, name_2, model_name){
  
  ## Correlation coeff
  rho <- cor(top_genes1$t, top_genes2$t, method = "spearman")
  rho_anno = paste0("rho = ", format(round(rho, 2), nsmall = 2))
  
  ## Merge data
  t_stats<-data.frame(t1=top_genes1$t, t2=top_genes2$t)
  ## Add DEG info for both groups
  t_stats$DEG<-add_DE_info(top_genes1, top_genes2, name_1, name_2)
  
  cols <- c("red", "#ffad73","#26b3ff", "dark grey") 
  names(cols)<-c("sig Both", paste0("sig ", name_1), paste0("sig ", name_2), "None")
  alphas <- c( 1, 1, 1,0.5)  
  names(alphas)<-c("sig Both", paste0("sig ", name_1), paste0("sig ", name_2), "None")
   
  plot <- ggplot(t_stats, aes(x = t1, y = t2, color=DEG, alpha=DEG)) +
     geom_point(size = 1) +
     labs(x = paste("t-stats", name_1), 
          y = paste("t-stats", name_2),
          title = model_name, 
          subtitle = rho_anno, 
          parse = T) +
     theme_bw() +
     scale_color_manual(values = cols)+ 
     scale_alpha_manual(values = alphas)

 return(plot)
}


## Function to create multiple t-stats plots
tstats_plots<-function(top_genes_pairs, name_1, name_2, models){
  
  plots<-list()
  for (i in 1:length(top_genes_pairs)){
    p<-t_stat_plot(top_genes_pairs[[i]][[1]], top_genes_pairs[[i]][[2]], name_1, name_2, models[i])
    plots[[i]]<-p
  }
  plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=2)
  ggsave(filename=paste("plots/04_DEA/02_Comparisons/t_stats_",gsub(" ", "_", name_1), "_VS_",
                        gsub(" ", "_", name_2), ".pdf", sep=""), 
                        height = 20, width = 25, units = "cm")
}




####################################
# Blood vs brain adults (smoking)
####################################

top_genes_pairs<-list(list(top_genes_blood_smoking_naive, top_genes_adults_smoking_naive), 
                      list(top_genes_blood_smoking_fitted, top_genes_adults_smoking_fitted),
                      list(top_genes_blood_smoking_interaction, top_genes_adults_smoking_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking blood", "Smoking brain", models)


#########################################
# Smoking adults vs nicotine adults 
#########################################

top_genes_pairs<-list(list(top_genes_adults_smoking_naive, top_genes_adults_nicotine_naive), 
                      list(top_genes_adults_smoking_fitted, top_genes_adults_nicotine_fitted),
                      list(top_genes_adults_smoking_interaction, top_genes_adults_nicotine_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking adults", "Nicotine adults", models)


#####################################
# Smoking pups vs nicotine pups 
#####################################

top_genes_pairs<-list(list(top_genes_pups_smoking_naive, top_genes_pups_nicotine_naive), 
                      list(top_genes_pups_smoking_fitted, top_genes_pups_nicotine_fitted),
                      list(top_genes_pups_smoking_interaction, top_genes_pups_nicotine_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking pups", "Nicotine pups", models)


#####################################
# Smoking pups vs smoking adults
#####################################

top_genes_pairs<-list(list(top_genes_pups_smoking_naive, top_genes_adults_smoking_naive), 
                      list(top_genes_pups_smoking_fitted, top_genes_adults_smoking_fitted),
                      list(top_genes_pups_smoking_interaction, top_genes_adults_smoking_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Smoking pups", "Smoking adults", models)


#####################################
# Nicotine pups vs nicotine adults
#####################################

top_genes_pairs<-list(list(top_genes_pups_nicotine_naive, top_genes_adults_nicotine_naive), 
                      list(top_genes_pups_nicotine_fitted, top_genes_adults_nicotine_fitted),
                      list(top_genes_pups_nicotine_interaction, top_genes_adults_nicotine_interaction))
models<-c("Naive model", "Fitted model", "Interaction model")
tstats_plots(top_genes_pairs,  "Nicotine pups", "Nicotine adults", models)


####################################
# Naive vs fitted nicotine pups
####################################

t<-t_stat_plot(top_genes_pups_nicotine_naive, top_genes_pups_nicotine_fitted, 
               "Naive model", "Fitted model", "Nicotine pups")
ggsave("plots/04_DEA/02_Comparisons/t_stats_Naive_VS_Fitted_Nicotine.pdf", t, 
       height = 10, width = 12, units = "cm")


####################################
# Naive vs fitted smoking pups
####################################

t<-t_stat_plot(top_genes_pups_smoking_naive, top_genes_pups_smoking_fitted, 
               "Naive model", "Fitted model", "Smoking pups")
ggsave("plots/04_DEA/02_Comparisons/t_stats_Naive_VS_Fitted_Smoking.pdf", t, 
       height = 10, width = 12, units = "cm")










### 1.2.3 Venn diagrams

## Function to create multiple Venn diagrams
venn_plot<-function(DEG_lists, colors, filename){
  plots<-list()
  pdf(file = paste("plots/04_DEA/02_Comparisons/Venn_", filename, ".pdf", sep=""))
  for (i in 1:length(DEG_lists)){
     v<-venn.diagram(DEG_lists[[i]], fill=colors[[i]], alpha = rep(0.5, length(DEG_lists[[i]])), 
                     lwd =0, margin=0.2, cat.cex=0.6, cex=0.6, height = 35, width = 35, units = "cm", 
                     cat.dist=rep(0.09, length(DEG_lists[[i]])),filename=NULL)
     plots[[i]]<-v
  }
  
  if (i==4){
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol=2)
    dev.off()
  }
  else{
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol=2)
    dev.off()
  }
}



################################################################################

## Compare naive VS fitted models DEG 
## Nicotine DEG
DEG_nic_naive_vs_fitted <-list(
  "Nicotine naive"=de_genes_pups_nicotine_naive$Symbol,
  "Nicotine fitted"=de_genes_pups_nicotine_fitted$Symbol)
## Smoking DEG
DEG_smo_naive_vs_fitted <-list(
  "Smoking naive"=de_genes_pups_smoking_naive$Symbol,
  "Smoking fitted"=de_genes_pups_smoking_fitted$Symbol)

## Compare smoking VS nicotine DEG
## Naive model DEG
DEG_naive_smo_vs_nic <-list(
  "Nicotine naive"=de_genes_pups_nicotine_naive$Symbol,
  "Smoking naive"=de_genes_pups_smoking_naive$Symbol)
## Fitted model DEG
DEG_fitted_smo_vs_nic <-list(
  "Nicotine fitted"=de_genes_pups_nicotine_fitted$Symbol,
  "Smoking fitted"=de_genes_pups_smoking_fitted$Symbol)


DEG_lists<-list(DEG_nic_naive_vs_fitted, DEG_smo_naive_vs_fitted, 
                DEG_naive_smo_vs_nic, DEG_fitted_smo_vs_nic  )
colors<-list(c("palegreen2", "yellow3"), c("lightsalmon", "slategray3"), 
             c("palegreen2", "lightsalmon"), c("yellow3", "slategray3"))
venn_plot(DEG_lists, colors, "model_smo_and_nic")


################################################################################

## Venn diagrams for Up and Down regulated DEG
## Datasets
nic_naive_up<-de_genes_pups_nicotine_naive[which(de_genes_pups_nicotine_naive$logFC>0),"Symbol"]
nic_naive_down<-de_genes_pups_nicotine_naive[which(de_genes_pups_nicotine_naive$logFC<0),"Symbol"]
nic_fitted_up<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC>0),"Symbol"]
nic_fitted_down<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC<0),"Symbol"]
smo_naive_up<-de_genes_pups_smoking_naive[which(de_genes_pups_smoking_naive$logFC>0),"Symbol"]
smo_naive_down<-de_genes_pups_smoking_naive[which(de_genes_pups_smoking_naive$logFC<0),"Symbol"]
smo_fitted_up<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC>0),"Symbol"]
smo_fitted_down<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC<0),"Symbol"]



## Compare Up and down Regulated nicotine/smoking genes 
## Naive VS fitted models

## Nicotine DEG
DEG_nic_up<-list(
  "Nicotine naive up"=nic_naive_up,
  "Nicotine fitted up"=nic_fitted_up)
DEG_nic_down<-list(
  "Nicotine naive down"=nic_naive_down,
  "Nicotine fitted down"=nic_fitted_down)

## Smoking DEG
DEG_smo_up<-list(
  "Smoking naive up"=smo_naive_up,
  "Smoking fitted up"=smo_fitted_up)
DEG_smo_down<-list(
  "Smoking naive down"=smo_naive_down,
  "Smoking fitted down"=smo_fitted_down)

DEG_lists<-list(DEG_nic_up, DEG_nic_down, DEG_smo_up, DEG_smo_down)
colors<-list(c("firebrick3", "brown1"), c("dodgerblue2", "deepskyblue1"), 
             c("coral1", "darksalmon"), c("darkturquoise", "darkslategray2"))
venn_plot(DEG_lists, colors, "nic_and_smo_Up_and_Down")




## Compare smoking VS nicotine Up and Down regulated genes
## For genes from either naive or fitted model
DEG_smo_vs_nic_Up<-list(
  "Smoking up"=union(smo_fitted_up, smo_naive_up),
  "Nicotine up"=union(nic_fitted_up, nic_naive_up))

DEG_smo_vs_nic_Down<-list(
  "Smoking down"=union(smo_fitted_down, smo_naive_down),
  "Nicotine down"=union(nic_fitted_down, nic_naive_down))

DEG_smoUp_vs_nicDown<-list(
  "Smoking up"=union(smo_fitted_up, smo_naive_up),
  "Nicotine down"=union(nic_fitted_down, nic_naive_down))

DEG_smoDown_vs_nicUp<-list(
  "Smoking down"=union(smo_fitted_down, smo_naive_down),
  "Nicotine up"=union(nic_fitted_up, nic_naive_up))

DEG_lists<-list(DEG_smo_vs_nic_Up, DEG_smo_vs_nic_Down, DEG_smoUp_vs_nicDown, DEG_smoDown_vs_nicUp)
colors<-list(c("coral", "brown3"), c("cyan4", "cyan3"), c("coral", "cyan3"), c("cyan4", "brown3"))
venn_plot(DEG_lists, colors, "smo_VS_nic_Up_and_Down")




## Compare smoking VS nicotine Up and Down regulated genes 
## For naive/fitted genes only

## Naive model DEG
DEG_naive_smo_vs_nic_up<-list(
  "Smoking naive up"=smo_naive_up,
  "Nicotine naive up"=nic_naive_up)
DEG_naive_smo_vs_nic_down<-list(
  "Smoking naive down"=smo_naive_down,
  "Nicotine naive down"=nic_naive_down)
DEG_naive_smoUp_nicDown<-list(
  "Smoking naive up"=smo_naive_up,
  "Nicotine naive down"=nic_naive_down)
DEG_naive_smoDown_nicUp<-list(
  "Smoking naive down"=smo_naive_down,
  "Nicotine naive up"=nic_naive_up)

DEG_lists<-list(DEG_naive_smo_vs_nic_up, DEG_naive_smo_vs_nic_down, 
                DEG_naive_smoUp_nicDown, DEG_naive_smoDown_nicUp)
colors=list(c("coral1", "firebrick3"), c("darkturquoise", "dodgerblue2"), 
            c("coral1", "dodgerblue2"), c("darkturquoise", "firebrick3"))
venn_plot(DEG_lists, colors, "Naive_smo_VS_nic_Up_and_Down")


## Fitted model DEG
DEG_fitted_smo_vs_nic_up<-list(
  "Smoking fitted up"=smo_fitted_up,
  "Nicotine fitted up"=nic_fitted_up)
save(DEG_fitted_smo_vs_nic_up, file="processed-data/04_DEA/DEG_fitted_smo_vs_nic_up.Rdata")

DEG_fitted_smo_vs_nic_down<-list(
  "Smoking fitted down"=smo_fitted_down,
  "Nicotine fitted down"=nic_fitted_down)
save(DEG_fitted_smo_vs_nic_down, file="processed-data/04_DEA/DEG_fitted_smo_vs_nic_down.Rdata")

DEG_fitted_smoUp_nicDown<-list(
  "Smoking fitted up"=smo_fitted_up,
  "Nicotine fitted down"=nic_fitted_down)
save(DEG_fitted_smoUp_nicDown, file = "processed-data/04_DEA/DEG_fitted_smoUp_nicDown.Rdata")

DEG_fitted_smoDown_nicUp<-list(
  "Smoking fitted down"=smo_fitted_down,
  "Nicotine fitted up"=nic_fitted_up)
save(DEG_fitted_smoDown_nicUp, file="processed-data/04_DEA/DEG_fitted_smoDown_nicUp.Rdata")

DEG_lists<-list(DEG_fitted_smo_vs_nic_up, DEG_fitted_smo_vs_nic_down, 
                DEG_fitted_smoUp_nicDown, DEG_fitted_smoDown_nicUp)
colors=list(c("darksalmon", "brown1"), c("darkslategray2", "deepskyblue1"),
            c("darksalmon", "deepskyblue1"), c("darkslategray2", "brown1"))
venn_plot(DEG_lists, colors, "Fitted_smo_VS_nic_Up_and_Down")


###############################################################################

## Compare all 4 groups of DEG (by expt and model)
DEG_all<-list(
  "Nicotine naive"=de_genes_pups_nicotine_naive$Symbol,
  "Nicotine fitted"=de_genes_pups_nicotine_fitted$Symbol,
  "Smoking naive"=de_genes_pups_smoking_naive$Symbol,
  "Smoking fitted"=de_genes_pups_smoking_fitted$Symbol)

## Compare DEG from the fitted model
## Compare all Up/Down regulated genes
DEG_all_Up_and_Down<-list(
  "Nicotine up"=nic_fitted_up,
  "Nicotine down"=nic_fitted_down,
  "Smoking up"=smo_fitted_up, 
  "Smoking down"=smo_fitted_down)

## Compare all smoking VS all nicotine DEG
DEG_all_smo_vs_all_nic<-list(
  "Nicotine"=union(de_genes_pups_nicotine_naive$Symbol, de_genes_pups_nicotine_fitted$Symbol),
  "Smoking"=union(de_genes_pups_smoking_naive$Symbol, de_genes_pups_smoking_fitted$Symbol))


DEG_lists<-list(DEG_all, DEG_all_Up_and_Down, DEG_all_smo_vs_all_nic)
colors<-list(c("palegreen2", "yellow3", "lightsalmon", "slategray3"), 
             c("brown3", "cyan3", "coral", "cyan4"), c("olivedrab1", "rosybrown2"))
venn_plot(DEG_lists, colors, "all_DEG")









## Reproducibility information

options(width = 120)
session_info()
# setting  value
# version  R version 4.2.0 (2022-04-22 ucrt)
# os       Windows 10 x64 (build 19044)
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  Spanish_Mexico.utf8
# ctype    Spanish_Mexico.utf8
# tz       America/Mexico_City
# date     2022-09-22
# rstudio  2022.02.3+492 Prairie Trillium (desktop)
# pandoc   2.17.1.1 @ C:/Program Files/RStudio/bin/quarto/bin/ (via rmarkdown)






