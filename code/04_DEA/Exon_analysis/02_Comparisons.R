
## 1.2 Comparison of DE exons

load(here("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_nic.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Exon_analysis/de_exons_smo.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_nicotine_fitted.Rdata"))
#load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
load(here("processed-data/04_DEA/Gene_analysis/de_genes_pups_smoking_fitted.Rdata"))
#load(here("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))
load(here("processed-data/05_GO_KEGG/Gene_analysis/intersections.Rdata"))     


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



################################
# Smoking vs nicotine (exons)
################################
t_stat_plot(top_exons_nic, top_exons_smo, "Nicotine pups", "Smoking pups", "Nic vs Smo DE exons")



##########################################
# Smoking genes vs smoking exons' genes ????
##########################################





### 1.2.3 Venn diagrams

## Function to create multiple Venn diagrams
venn_plot<-function(DE_lists, colors, name, titles){
  
  if (name=="smo_VS_nic_DE_exons"){
    height=10
    width=15
    margin=0.2
    dist=0.1
    cat=0.75
    cex=0.8
  }
  
  else if (name=="DEG_VS_exons_genes"){
    height=12
    width=16
    margin=5
    dist=0.07
    cat=0.75
    cex=0.8
  }
  
  else if (name=="smo_VS_nic_DE_exons_genes"){
    height=12
    width=16
    margin=0.2
    dist=0.09
    cat=0.75
    cex=0.8
  }
  
  else if (name=="intersections_DEG_VS_exons_genes"){
    height=10
    width=20
    margin=5
    dist=0.06
    cat=0.9
    cex=1
  }
  
  plots<-list()
  pdf(file = paste("plots/04_DEA/02_Comparisons/Exon_analysis/Venn_", name, ".pdf", sep=""), 
      height = height, width = width)
  for (i in 1:length(DE_lists)){
    v<-venn.diagram(DE_lists[[i]], fill=colors[[i]], alpha = rep(0.5, length(DE_lists[[i]])), 
                    lwd =0, margin=margin, cat.cex=cat, cex=cex, height = 10, width = 12.5, units = "cm", 
                    cat.dist=rep(dist, length(DE_lists[[i]])), filename=NULL, main = titles[i], 
                    main.pos = c(0.5,0.476), disable.logging=TRUE)
    plots[[i]]<-v
  }
  
  if (i==6){
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], 
                            plots[[6]], ncol=3)
    dev.off()
  }
  else if (i==9) {
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
                            plots[[6]], plots[[7]], plots[[8]], plots[[9]], ncol=4)
    dev.off()
  }
  else {
    gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]],
                            plots[[6]], plots[[7]], plots[[8]], ncol=4)
    dev.off()    
  }
}




## Define groups of DE exons
nic_up<-de_exons_nic[which(de_exons_nic$logFC>0),"exon_libdID"]
smo_up<-de_exons_smo[which(de_exons_smo$logFC>0),"exon_libdID"]
nic_down<-de_exons_nic[which(de_exons_nic$logFC<0),"exon_libdID"]
smo_down<-de_exons_smo[which(de_exons_smo$logFC<0),"exon_libdID"]


## Define groups of DE exons' genes
nic_up_genes<-unique(de_exons_nic[which(de_exons_nic$logFC>0),"Symbol"])
smo_up_genes<-unique(de_exons_smo[which(de_exons_smo$logFC>0),"Symbol"])
nic_down_genes<-unique(de_exons_nic[which(de_exons_nic$logFC<0),"Symbol"])
smo_down_genes<-unique(de_exons_smo[which(de_exons_smo$logFC<0),"Symbol"])

smoUp_nicUp_genes<-intersect(nic_up_genes, smo_up_genes)
smoDown_nicDown_genes<-intersect(nic_down_genes, smo_down_genes)
smoUp_nicDown_genes<-intersect(nic_down_genes, smo_up_genes)
smoDown_nicUp_genes<-intersect(nic_up_genes, smo_down_genes)
only_up_nic_genes<-nic_up_genes[which(! (nic_up_genes %in% smo_up_genes | 
                            nic_up_genes %in% smo_down_genes | 
                            nic_up_genes %in% nic_down_genes))]
only_up_smo_genes<-smo_up_genes[which(! (smo_up_genes %in% smo_down_genes | 
                            smo_up_genes %in% nic_down_genes | 
                            smo_up_genes %in% nic_up_genes))]
only_down_nic_genes<-nic_down_genes[which(! (nic_down_genes %in% smo_up_genes | 
                              nic_down_genes %in% smo_down_genes | 
                              nic_down_genes %in% nic_up_genes))]
only_down_smo_genes<-smo_down_genes[which(! (smo_down_genes %in% smo_up_genes | 
                              smo_down_genes %in% nic_down_genes | 
                              smo_down_genes %in% nic_up_genes))]


## Define groups of DEG
nic_DEG_up<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC>0),"Symbol"]
nic_DEG_down<-de_genes_pups_nicotine_fitted[which(de_genes_pups_nicotine_fitted$logFC<0),"Symbol"]
smo_DEG_up<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC>0),"Symbol"]
smo_DEG_down<-de_genes_pups_smoking_fitted[which(de_genes_pups_smoking_fitted$logFC<0),"Symbol"]

only_up_nic_DEG<-intersections[[1]]$Symbol
only_up_smo_DEG<-intersections[[2]]$Symbol
only_down_nic_DEG<-intersections[[3]]$Symbol
only_down_smo_DEG<-intersections[[4]]$Symbol
smoUp_nicUp_DEG<-intersections[[5]]$Symbol
smoDown_nicDown_DEG<-intersections[[6]]$Symbol
smoUp_nicDown_DEG<-intersections[[7]]$Symbol
smoDown_nicUp_DEG<-intersections[[8]]$Symbol





## Compare smoking VS nicotine DE exons

## Smo vs nic DE exons
DE_exons_smo_vs_nic <-list(
  "Nicotine"=de_exons_nic$exon_libdID,
  "Smoking"=de_exons_smo$exon_libdID)

## Smo vs nic Up DE exons
DE_exons_smo_vs_nic_Up <-list(
  "Nicotine up"=nic_up,
  "Smoking up"=smo_up)

## Smo vs nic Down DE exons
DE_exons_smo_vs_nic_Down <-list(
  "Nicotine down"=nic_down,
  "Smoking down"=smo_down)

## Smo Up vs nic Down DE exons
DE_exons_smoUp_vs_nicDown <-list(
  "Nicotine down"=nic_down,
  "Smoking up"=smo_up)

## Smo Down vs nic Up DE exons
DE_exons_smoDown_vs_nicUp <-list(
  "Nicotine up"=nic_up,
  "Smoking down"=smo_down)

## All DE exons
DE_exons_all<-list(
  "Nicotine up"=nic_up,
  "Nicotine down"=nic_down,
  "Smoking up"=smo_up, 
  "Smoking down"=smo_down)



## Venn diagrams
DE_lists<-list(DE_exons_smo_vs_nic, DE_exons_smo_vs_nic_Up, DE_exons_smo_vs_nic_Down, 
              DE_exons_smoUp_vs_nicDown, DE_exons_smoDown_vs_nicUp, DE_exons_all)
colors<-list(c("olivedrab1", "rosybrown2"), c("brown3", "coral"), c("cyan3", "cyan4"), 
             c("cyan3", "coral"), c("brown3", "cyan4"), c("brown3", "cyan3", "coral", "cyan4"))
venn_plot(DE_lists, colors, "smo_VS_nic_DE_exons", NULL)





## Compare DE exons' genes 

## Smo vs nic genes
DE_exons_genes_smo_vs_nic <-list(
  "Nicotine"=unique(de_exons_nic$Symbol),
  "Smoking"=unique(de_exons_smo$Symbol))

## Smo vs nic Up genes
DE_exons_genes_smo_vs_nic_Up <-list(
  "Nicotine up"=nic_up_genes,
  "Smoking up"=smo_up_genes)

## Smo vs nic Down genes
DE_exons_genes_smo_vs_nic_Down <-list(
  "Nicotine down"=nic_down_genes,
  "Smoking down"=smo_down_genes)

## Smo Up vs nic Down genes
DE_exons_genes_smoUp_vs_nicDown <-list(
  "Nicotine down"=nic_down_genes,
  "Smoking up"=smo_up_genes)

## Smo Down vs nic Up genes
DE_exons_genes_smoDown_vs_nicUp <-list(
  "Nicotine up"=nic_up_genes,
  "Smoking down"=smo_down_genes)

## Up vs Down genes
DE_exons_genes_Up_vs_Down <-list(
  "Up"=union(nic_up_genes, smo_up_genes),
  "Down"=union(nic_down_genes, smo_down_genes))

## Nic Up vs nic Down genes
DE_exons_genes_nicUp_vs_nicDown <-list(
  "Nicotine up"=nic_up_genes,
  "Nicotine down"=nic_down_genes)

## Smo Up vs smo Down genes
DE_exons_genes_smoUp_vs_smoDown <-list(
  "Smoking up"=smo_up_genes,
  "Smoking down"=smo_down_genes)

## All DE exons
DE_exons_genes_all<-list(
  "Nicotine up"=nic_up_genes,
  "Nicotine down"=nic_down_genes,
  "Smoking up"=smo_up_genes, 
  "Smoking down"=smo_down_genes)



## Venn diagrams
DE_lists<-list(DE_exons_genes_smo_vs_nic, DE_exons_genes_smo_vs_nic_Up, 
               DE_exons_genes_smo_vs_nic_Down, DE_exons_genes_smoUp_vs_nicDown, 
               DE_exons_genes_smoDown_vs_nicUp, DE_exons_genes_Up_vs_Down,
               DE_exons_genes_nicUp_vs_nicDown, DE_exons_genes_smoUp_vs_smoDown, 
               DE_exons_genes_all)
colors<-list(c("olivedrab1", "rosybrown2"), c("brown3", "coral"), 
             c("cyan3", "cyan4"), c("cyan3", "coral"), c("brown3", "cyan4"), 
             c("firebrick", "dodgerblue"),c("brown3", "cyan3"), c("coral", "cyan4"),
             c("brown3", "cyan3", "coral", "cyan4"))
venn_plot(DE_lists, colors, "smo_VS_nic_DE_exons_genes", NULL)





## Compare DE exons' genes with DEG

## All DEG vs all exons' genes
DEG_vs_Exons_all <-list(
  "DEG"= union(union(nic_DEG_down, nic_DEG_up), union(smo_DEG_down, smo_DEG_up)),
  "DE exons' genes"=union(union(nic_down_genes, nic_up_genes), union(smo_down_genes, smo_up_genes))
)

## Nic DEGs vs nic exons' genes 
DEG_nic_vs_Exons_nic <-list(
  "DEG"= union(nic_DEG_down, nic_DEG_up),
  "DE exons' genes"=union(nic_down_genes, nic_up_genes)
)

## Smo DEGs vs mo exons' genes
DEG_smo_vs_Exons_smo <-list(
  "DEG"= union(smo_DEG_down, smo_DEG_up),
  "DE exons' genes"=union(smo_down_genes, smo_up_genes)
)

## Up DEGs vs up exons' genes
DEG_up_vs_Exons_up <-list(
  "DEG"= union(nic_DEG_up, smo_DEG_up),
  "DE exons' genes"=union(nic_up_genes, smo_up_genes)
)

## Down DEGs vs down exons' genes
DEG_down_vs_Exons_down <-list(
  "DEG"= union(nic_DEG_down, smo_DEG_down),
  "DE exons' genes"=union(nic_down_genes, smo_down_genes)
)

## Nic up DEGs vs nic up exons' genes
DEG_nicUp_vs_Exons_nicUp <-list(
  "DEG"=nic_up_genes,
  "DE exons' genes"=nic_DEG_up
)

## Nic down DEGs vs nic down exons' genes
DEG_nicDown_vs_Exons_nicDown <-list(
  "DEG"=nic_down_genes,
  "DE exons' genes"=nic_DEG_down
)

## Smo up DEGs vs smo up exons' genes
DEG_smoUp_vs_Exons_smoUp <-list(
  "DEG"=smo_up_genes,
  "DE exons' genes"=smo_DEG_up
)

## Smo down DEGs vs smo down exons' genes
DEG_smoDown_vs_Exons_smoDown <-list(
  "DEG"=smo_down_genes,
  "DE exons' genes"=smo_DEG_down
)



## Venn diagrams
DE_lists<-list(DEG_vs_Exons_all, DEG_nic_vs_Exons_nic, DEG_smo_vs_Exons_smo,
               DEG_up_vs_Exons_up, DEG_down_vs_Exons_down, DEG_nicUp_vs_Exons_nicUp,
               DEG_nicDown_vs_Exons_nicDown, DEG_smoUp_vs_Exons_smoUp, DEG_smoDown_vs_Exons_smoDown)
colors<-list(c("rosybrown2", "navajowhite2"), c("plum", "lightgoldenrod2"), 
             c("pink1", "khaki2"), c("rosybrown4", "navajowhite4"), c("rosybrown1", "navajowhite1"), 
             c("plum3", "lightgoldenrod4"),c("plum1", "lightgoldenrod1"), c("pink3", "khaki3"),
             c("pink", "khaki1"))
venn_plot(DE_lists, colors, "DEG_VS_exons_genes", c("All", "Nicotine", "Smoking", "Up", "Down", 
                                                    "Nicotine Up", "Nicotine Down", "Smoking Up", "Smoking Down"))





## Compare DEG intersections with those of exons' genes

## Only up nic DEGs vs exons' genes 
DEG_vs_Exons_only_up_nic <-list(
  "DEG"= only_up_nic_DEG,
  "DE exons' genes"= only_up_nic_genes
)

## Only down nic DEGs vs exons' genes 
DEG_vs_Exons_only_down_nic <-list(
  "DEG"= only_down_nic_DEG,
  "DE exons' genes"= only_down_nic_genes
)

## Only up smo DEGs vs exons' genes 
DEG_vs_Exons_only_up_smo <-list(
  "DEG"= only_up_smo_DEG,
  "DE exons' genes"= only_up_smo_genes
)

## Only down smo DEGs vs exons' genes 
DEG_vs_Exons_only_down_smo <-list(
  "DEG"= only_down_smo_DEG,
  "DE exons' genes"= only_down_smo_genes
)

## Smo and nic up DEGs vs exons' genes
DEG_vs_Exons_smoUp_nicUp <-list(
  "DEG"=smoUp_nicUp_DEG,
  "DE exon's genes"=smoUp_nicUp_genes
)

## Smo and nic down DEGs vs exons' genes
DEG_vs_Exons_smoDown_nicDown <-list(
  "DEG"=smoDown_nicDown_DEG,
  "DE exon's genes"=smoDown_nicDown_genes
)

## Smo up and nic down DEGs vs exons' genes
DEG_vs_Exons_smoUp_nicDown <-list(
  "DEG"=smoUp_nicDown_DEG,
  "DE exon's genes"=smoUp_nicDown_genes
)

## Smo down and nic up DEGs vs exons' genes
DEG_vs_Exons_smoDown_nicUp <-list(
  "DEG"=smoDown_nicUp_DEG,
  "DE exon's genes"=smoDown_nicUp_genes
)



## Venn diagrams
DE_lists<-list(DEG_vs_Exons_only_up_nic, DEG_vs_Exons_only_down_nic, DEG_vs_Exons_only_up_smo,
               DEG_vs_Exons_only_down_smo, DEG_vs_Exons_smoUp_nicUp, DEG_vs_Exons_smoDown_nicDown,
               DEG_vs_Exons_smoUp_nicDown, DEG_vs_Exons_smoDown_nicUp)
colors<-list(c("darkslategray3", "lightgoldenrod3"), c("darkslategray1", "lightgoldenrod1"), 
             c("lightpink3", "olivedrab3"), c("pink", "olivedrab1"), c("palegreen3", "navajowhite3"),
             c("palegreen", "navajowhite"), c("palegreen","navajowhite3"), c("palegreen3", "navajowhite"))
venn_plot(DE_lists, colors, "intersections_DEG_VS_exons_genes", c("Only up in nic", "Only down in nic",
                                                                  "Only up in smo", "Only down in smo", 
                                                                  "Smo up, nic up", "Smo down, nic down",
                                                                  "Smo up, nic down", "Smo down, nic up"))





