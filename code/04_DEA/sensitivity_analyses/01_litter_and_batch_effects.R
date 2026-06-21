# 5. Sensitivity analyses to assess hidden litter and batch-related effects

library(here)
library(tidyr)
library(dplyr)
library(SummarizedExperiment)
library(pheatmap)
library(VennDiagram) 
library(ggplot2)
library(cowplot)
library(limma)
library(edgeR)
library(sva)
library(sessioninfo)


## Load data at 4 expr levels including filtered pup samples 
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_nicotine.Rdata"), verbose = T)
load(here("processed-data/03_EDA/04_Expl_Var_partition/rse_gene_brain_pups_smoking.Rdata"))

load(here("processed-data/03_EDA/03_PCA/rse_gene_brain_pups_qc_afterPCA.Rdata"), verbose = T)
load(here("processed-data/03_EDA/03_PCA/rse_tx_brain_pups_qc_afterPCA.Rdata"), verbose = T)
load(here("processed-data/03_EDA/03_PCA/rse_exon_brain_pups_qc_afterPCA.Rdata"), verbose = T)
load(here("processed-data/03_EDA/03_PCA/rse_jx_brain_pups_qc_afterPCA.Rdata"), verbose = T)


## Split by expt and group (independent litters)
rse_gene_nic_ctrl <- rse_gene_brain_pups_qc_afterPCA[, which(rse_gene_brain_pups_qc_afterPCA$Expt == "Nicotine" &
                                                             rse_gene_brain_pups_qc_afterPCA$Group == "Control")]
rse_gene_nic_expt <- rse_gene_brain_pups_qc_afterPCA[, which(rse_gene_brain_pups_qc_afterPCA$Expt == "Nicotine" &
                                                             rse_gene_brain_pups_qc_afterPCA$Group == "Experimental")]
rse_gene_smo_ctrl <- rse_gene_brain_pups_qc_afterPCA[, which(rse_gene_brain_pups_qc_afterPCA$Expt == "Smoking" & 
                                                             rse_gene_brain_pups_qc_afterPCA$Group == "Control")]
rse_gene_smo_expt <- rse_gene_brain_pups_qc_afterPCA[, which(rse_gene_brain_pups_qc_afterPCA$Expt == "Smoking" &
                                                             rse_gene_brain_pups_qc_afterPCA$Group == "Experimental")]

## For DTE, DEE and DJE
rse_tx_nic <- rse_tx_brain_pups_qc_afterPCA[, which(rse_tx_brain_pups_qc_afterPCA$Expt == "Nicotine")]
rse_tx_smo <- rse_tx_brain_pups_qc_afterPCA[, which(rse_tx_brain_pups_qc_afterPCA$Expt == "Smoking")]

rse_exon_nic <- rse_exon_brain_pups_qc_afterPCA[, which(rse_exon_brain_pups_qc_afterPCA$Expt == "Nicotine")]
rse_exon_smo <- rse_exon_brain_pups_qc_afterPCA[, which(rse_exon_brain_pups_qc_afterPCA$Expt == "Smoking")]

rse_jx_nic <- rse_jx_brain_pups_qc_afterPCA[, which(rse_jx_brain_pups_qc_afterPCA$Expt == "Nicotine")]
rse_jx_smo <- rse_jx_brain_pups_qc_afterPCA[, which(rse_jx_brain_pups_qc_afterPCA$Expt == "Smoking")]



# ------------------------------------------------------------------------------
#              5.1 Correlation in gene expression between pups
# ------------------------------------------------------------------------------

plot_cor_heatmap <- function(expt, group){
  
  rse <- get(paste0("rse_gene_", expt, "_", group))
  cor_mat <- cor(assays(rse)$logcounts, method = "pearson")
  dist_mat <- as.dist(1 - abs(cor_mat))
  colnames(cor_mat) <- rownames(cor_mat) <- rse$SAMPLE_ID
  
  ann_cols <- colData(rse)[, c("Group", "Sex")]
  rownames(ann_cols) <- rse$SAMPLE_ID
  colnames(ann_cols) <- c("Group", "Sex")
  ann_cols  <- as.data.frame(ann_cols)
  
  ann_colors = list()
  ann_colors[["Sex"]]=c("F"="hotpink1", "M"="dodgerblue")
  ann_colors[["Group"]]=c("Control"="seashell3", "Experimental"="orange3")
  
  pheatmap(
    cor_mat, 
    colorRampPalette(c("#FFE4C4", "#FF7F00", "#EE7600"))(5), 
    border_color = NA, 
    clustering_distance_rows = dist_mat, 
    clustering_distance_cols = dist_mat, 
    clustering_method = "average",
    show_rownames = FALSE,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE, 
    # cutree_rows = 3, 
    # cutree_cols = 3, 
    treeheight_col = 10, 
    treeheight_row = 10,
    annotation_col = ann_cols, 
    annotation_colors = ann_colors, 
    cellwidth = 5,
    cellheight = 5,
    width = 5,
    height = 5,
    filename = paste("plots/04_DEA/sensitivity_analyses/Corr_heatmap_", expt, "_", group, ".pdf", sep="")
  )
  
}

plot_cor_heatmap("nic", "ctrl")
plot_cor_heatmap("nic", "expt")
plot_cor_heatmap("smo", "ctrl")
plot_cor_heatmap("smo", "expt")



# ------------------------------------------------------------------------------
#        5.2 PCA and pup K-means clustering based on gene expression
# ------------------------------------------------------------------------------

colors = list("Group"=c("Control" = "seashell3", "Experimental" = "orange3"),
              "Sex"=c("F" = "hotpink1", "M" = "dodgerblue"),
              "plate"=c("Plate1" = "darkorange", "Plate2" = "lightskyblue", "Plate3" = "deeppink1"),
              "flowcell"=c("HKCG7DSXX" = "chartreuse2", "HKCMHDSXX" = "magenta",  
                           "HKCNKDSXX" = "turquoise3", "HKCTMDSXX" = "tomato", 
                           "HK7JHDSXX"="seagreen3","HKCJCDSXX"="palevioletred2"),
              "cluster" = c("1" = "red", "2" = "#EEEE00", "3" = "#7FFFD4", 
                            "4" = "#7FFF00", "5" = "#AB82FF", "6" = "#00C5CD", 
                            "7" = "#436EEE", "8" = "#FFBBFF"))

shapes = list("Group" = c("Experimental" = 16, "Control" = 1),
              "Sex" = c("F" = 5, "M" = 4),
              "plate" = c("Plate1" = 0, "Plate2" = 11, "Plate3" = 2),
              "flowcell" = c("HKCG7DSXX" = 7, "HKCMHDSXX" = 6, "HKCNKDSXX" = 3,
                           "HKCTMDSXX" = 1, "HK7JHDSXX" = 10, "HKCJCDSXX" = 12),
              "cluster" = c("1" = 15, "2" = 6, "3" = 8, 
                            "4" = 14, "5" = 11, "6" = 10, "7" = 16, "8" = 5))

## Max number of litters (clusters) per group (i.e. number of mothers)
num_centers = c("nic_ctrl" = 3,
                "nic_expt" = 3,
                "smo_ctrl" = 7,
                "smo_expt" = 8)

pca_kmeans_clust_resid <- function(expt, group, PCx, PCy, color_var, shape_var){
  
  rse <- get(paste0("rse_gene_", expt, "_", group))
  
  ## Compute TMM factors
  rse_norm <- calcNormFactors(rse, method = "TMM")
  
  ## Regress out cov (all as fixed effects) -- no Group effect
  f <- ~ Sex + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + 
    ERCCsumLogErr + mitoRate
  ## Model matrix
  model = model.matrix(f, data = rse_norm$samples)

  ## voom log-cpm and variance weights for lm fit
  v = voom(rse_norm, design = model, plot = TRUE)
  
  ## Fit linear model for each gene
  fitGene = lmFit(v)
  
  ## Residuals (include intercept = baseline gene expr)
  X_cov <- model[, -1, drop = FALSE]   # remove intercept column
  cov_effects <- fitGene$coefficients[, -1, drop = FALSE] %*% t(X_cov)
  residual_expr <- v$E - cov_effects
  
  ## PCs on residuals
  pca <- prcomp(t(residual_expr), center = TRUE, scale. = FALSE)
  var_per <- pca$sdev^2 / sum(pca$sdev^2)
  var_exp <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  ## PCs explaining 80% of variance
  nPC <- which(var_exp > 0.8)[1]
  pc_mat <- pca$x[,1:nPC]
  
  ## K mean clustering
  set.seed(123)
  km <- kmeans(pc_mat, centers = num_centers[paste0(expt, "_", group)], iter.max = 100)
  
  ## Add sample data and cluster
  pc_mat <- cbind(pc_mat, colData(rse))
  pc_mat$cluster = as.character(km$cluster)
  pc_mat <- as.data.frame(pc_mat)
  
  clusters <- as.character(km$cluster)
  
  plot = ggplot(data = pc_mat, 
                aes(x = get(PCx), y = get(PCy))) + 
    geom_point(aes(color = get(color_var), shape = get(shape_var)), 
               size = 2, stroke = 1) + 
    scale_color_manual(values = colors[[color_var]]) +
    scale_shape_manual(values = shapes[[shape_var]]) +
    theme_classic() + 
    labs(x = paste0(PCx, " (", signif(var_per[strtoi(gsub("PC","", PCx))], 2)*100, "%)"),
         y = paste0(PCy, " (", signif(var_per[strtoi(gsub("PC","", PCy))], 2)*100, "%)"),
         color = 'K-means cluster', shape = shape_var) +
    theme(legend.position = "right",
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          axis.title = element_text(size = 10)) 
  
  return(list(plot, clusters))
}
 

for(expt in c("nic", "smo")){
  for(group in c("expt", "ctrl")){
    for(cov in c("cluster", "Sex", "plate", "flowcell")){
      p1 <- pca_kmeans_clust_resid(expt, group, "PC1", "PC2", cov, "Sex")[[1]]
      p2 <- pca_kmeans_clust_resid(expt, group, "PC1", "PC3", cov, "Sex")[[1]]
      p3 <- pca_kmeans_clust_resid(expt, group, "PC2", "PC3", cov, "Sex")[[1]]
      p4 <- pca_kmeans_clust_resid(expt, group, "PC4", "PC5", cov, "Sex")[[1]]
      plot_grid(p1, p2, p3, p4, ncol = 2)
      ggsave(paste("plots/04_DEA/sensitivity_analyses/Resid_PCs_clusters_", 
                   cov, "_", expt, "_", group,".pdf", sep=""), width = 8, height = 5.5)
    } 
  }
}

## Add cluster info
rse_gene_nic_ctrl$cluster = paste0(pca_kmeans_clust_resid("nic", "ctrl", "PC1", "PC2", "cluster", "Sex")[[2]], "_nic_ctrl")
rse_gene_nic_expt$cluster = paste0(pca_kmeans_clust_resid("nic", "expt", "PC1", "PC2", "cluster", "Sex")[[2]], "_nic_expt")
rse_gene_smo_ctrl$cluster = paste0(pca_kmeans_clust_resid("smo", "ctrl", "PC1", "PC2", "cluster", "Sex")[[2]], "_smo_ctrl")
rse_gene_smo_expt$cluster = paste0(pca_kmeans_clust_resid("smo", "expt", "PC1", "PC2", "cluster", "Sex")[[2]], "_smo_expt")

## Bind ctrls and expt in each group
identical(rownames(rse_gene_nic_ctrl), rownames(rse_gene_nic_expt))
identical(colnames(colData(rse_gene_nic_ctrl)), colnames(colData(rse_gene_nic_expt)))
rse_gene_nic <- cbind(rse_gene_nic_ctrl, rse_gene_nic_expt)

identical(rownames(rse_gene_smo_ctrl), rownames(rse_gene_smo_expt))
identical(colnames(rse_gene_smo_ctrl), colnames(rse_gene_smo_expt))
rse_gene_smo <- cbind(rse_gene_smo_ctrl, rse_gene_smo_expt)

all_clusters <- c(unique(rse_gene_nic$cluster), unique(rse_gene_smo$cluster))
cluster_letters <- LETTERS[1:length(all_clusters)]
names(cluster_letters) <- all_clusters

rse_gene_nic$cluster <- cluster_letters[rse_gene_nic$cluster]
rse_gene_smo$cluster <- cluster_letters[rse_gene_smo$cluster]

colors_clusters <- c("#FF4500", "#9AFF9A", "#9B30FF", "#87CEEB", "#EEEE00", 
                     "#FFBBFF", "#2E8B57", "#B3EE3A", "#EE00EE", "#EEAD0E", 
                     "#FF7F00", "#76EEC6","#0000FF", "#FFAEB9", "#EECFA1", 
                     "#00E5EE", "#CDBA96", "#473C8B", "#FFA07A", "#66CD00", 
                     "#A52A2A")

names(colors_clusters) <- cluster_letters
colors$cluster <- colors_clusters
cluster_letter_shapes <- rep(16, length(cluster_letters))
names(cluster_letter_shapes) <- cluster_letters
shapes$cluster <- cluster_letter_shapes

## Plot PCs on all expt pups
plot_PCs_all_pups <- function(expt, PCx, PCy, color_var, shape_var){
  
  rse <- get(paste0("rse_gene_", expt))
  expr <- assays(rse)$logcounts
    
  ## PCs on lognorm counts
  pca <- prcomp(t(expr), center = TRUE, scale. = FALSE)
  var_per <- pca$sdev^2 / sum(pca$sdev^2)
  var_exp <- cumsum(pca$sdev^2 / sum(pca$sdev^2))

  ## Add sample data 
  pc_mat <- cbind(pca$x, colData(rse))
  pc_mat <- as.data.frame(pc_mat)
  
  plot = ggplot(data = pc_mat, 
                aes(x = get(PCx), y = get(PCy))) + 
    geom_point(aes(color = get(color_var), shape = get(shape_var)), 
               size = 2, stroke = 1) + 
    scale_color_manual(values = colors[[color_var]]) +
    scale_shape_manual(values = shapes[[shape_var]]) +
    theme_classic() + 
    labs(x = paste0(PCx, " (", signif(var_per[strtoi(gsub("PC","", PCx))], 2)*100, "%)"),
         y = paste0(PCy, " (", signif(var_per[strtoi(gsub("PC","", PCy))], 2)*100, "%)"),
         color = 'K-means cluster', shape = shape_var) +
    theme(legend.position = "right",
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          legend.key.height = unit(4, units = "mm"),
          axis.title = element_text(size = 10)) 
  
  if(shape_var == "cluster"){
    plot <- plot + guides(shape = "none")
  }
  
  return(plot)
  
}

## PNE PCs colored by cluster
for(cov in c("cluster", "Group", "Sex", "plate", "flowcell")){
  plots <- list()
  c = 1
  for(i in 1:5){
    for(j in (i+1):6){
      
      PCx = paste0("PC", i)
      PCy = paste0("PC", j)
      
      plots[[c]] <- plot_PCs_all_pups("nic", PCx, PCy, "cluster", cov)
      c = c + 1
    }
  }
  
  plot_grid(plotlist = plots, ncol = 5, align = "vh")
  ggsave(filename = paste0("plots/04_DEA/sensitivity_analyses/PNE_PCs_cluster_and_",
                           cov, ".pdf"), width = 21, height = 8)
}

## MSDP PCs colored by cluster
for(cov in c("cluster", "Group", "Sex", "plate", "flowcell")){
  plots <- list()
  c = 1
  for(i in 1:5){
    for(j in (i+1):6){
      
      PCx = paste0("PC", i)
      PCy = paste0("PC", j)
      
      plots[[c]] <- plot_PCs_all_pups("smo", PCx, PCy, "cluster", cov)
      c = c + 1
    }
  }
  
  plot_grid(plotlist = plots, ncol = 5, align = "vh")
  ggsave(filename = paste0("plots/04_DEA/sensitivity_analyses/MSDP_PCs_cluster_and_",
                           cov, ".pdf"), width = 21, height = 8)
}


## Fit random-effects model to assess DGE using clusters as proxy for litter
fit_lmm <- function(expt){
  
  rse <- get(paste0("rse_gene_", expt))

  ## TMM norm factors
  rse_norm <- calcNormFactors(rse, method = "TMM")
  
  ## Same model for PNE and MSDP as in fixed-effects model (04_DEA)
  f <- ~ Group + Sex + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate

  model = model.matrix(f, data = rse_norm$samples)
  
  ## log-cpm and inverse variance weights
  v = voom(rse_norm, design = model, plot=TRUE)
  
  ## Intra-cluster corr based on lognorm expr
  cor = duplicateCorrelation(v, design = model, block = rse_norm$samples$cluster)
  
  ## Re-compute voom weights accounting for intra-cluster corr
  v2 = voom(rse_norm, design = model, plot=TRUE, block = rse_norm$samples$cluster, 
            correlation = cor$consensus)
  
  ## Corr based on corrected expression
  cor2 = duplicateCorrelation(v2, design = model, block = rse_norm$samples$cluster)
  
  ## Fit linear model
  fit = lmFit(v2, design = model, block = rse_norm$samples$cluster, 
              correlation = cor2$consensus)
  eBGene = eBayes(fit)
  
  top_genes = topTable(eBGene, coef = "GroupExperimental", p.value = 1, 
                       number = nrow(rse), sort.by = "none")
  
  ## DEGs?
  de_genes <- subset(top_genes, adj.P.Val<0.05)
  de_genes <- de_genes[order(de_genes$adj.P.Val, decreasing = F), ]
  message("DGEs:")
  print(paste(dim(de_genes)[1], "DEGs:", 
              dim(subset(de_genes, logFC>0))[1], "up-regulated and", 
              dim(subset(de_genes, logFC<0))[1], "down-regulated"))
  
  return(top_genes)
}  

## DEGs
top_genes_nic_cluster_adjusted <- fit_lmm("nic")
# [1] "424 DEGs: 303 up-regulated and 121 down-regulated"
top_genes_smo_cluster_adjusted <- fit_lmm("smo")
# [1] "1032 DEGs: 380 up-regulated and 652 down-regulated"

de_genes_nic_cluster_adjusted <- subset(top_genes_nic_cluster_adjusted, adj.P.Val<0.05)
de_genes_nic_cluster_adjusted <- de_genes_nic_cluster_adjusted[order(de_genes_nic_cluster_adjusted$adj.P.Val, decreasing = F), ]

de_genes_smo_cluster_adjusted <- subset(top_genes_smo_cluster_adjusted, adj.P.Val<0.05)
de_genes_smo_cluster_adjusted <- de_genes_smo_cluster_adjusted[order(de_genes_smo_cluster_adjusted$adj.P.Val, decreasing = F), ]

save(top_genes_nic_cluster_adjusted, file="processed-data/04_DEA/sensitivity_analyses/top_genes_nic_cluster_adjusted.Rdata")
save(top_genes_smo_cluster_adjusted, file="processed-data/04_DEA/sensitivity_analyses/top_genes_smo_cluster_adjusted.Rdata")
save(de_genes_nic_cluster_adjusted, file="processed-data/04_DEA/sensitivity_analyses/de_genes_nic_cluster_adjusted.Rdata")
save(de_genes_smo_cluster_adjusted, file="processed-data/04_DEA/sensitivity_analyses/de_genes_smo_cluster_adjusted.Rdata")


## Check overlap with found DEGs
top_genes_nic_unadjusted <- get(load("processed-data/04_DEA/Gene_analysis/top_genes_pups_nicotine_fitted.Rdata"))
top_genes_smo_unadjusted <- get(load("processed-data/04_DEA/Gene_analysis/top_genes_pups_smoking_fitted.Rdata"))

de_genes_nic_unadjusted <- subset(top_genes_nic_unadjusted, adj.P.Val<0.05)
de_genes_smo_unadjusted <- subset(top_genes_smo_unadjusted, adj.P.Val<0.05)
nom_genes_nic_unadjusted <- subset(top_genes_nic_unadjusted, P.Value<0.05)
nom_genes_smo_unadjusted <- subset(top_genes_smo_unadjusted, P.Value<0.05)

DEG_nic<-list(
  "PNE adjusted DEGs" = de_genes_nic_cluster_adjusted$ensemblID,
  "PNE unadjusted DEGs"= de_genes_nic_unadjusted$ensemblID)
venn.diagram(DEG_nic, disable.logging = T, filename = NULL)

DEG_smo <- list(
  "MSDP adjusted DEGs" = de_genes_smo_cluster_adjusted$ensemblID,
  "MSDP unadjusted DEGs"= de_genes_smo_unadjusted$ensemblID)
venn.diagram(DEG_smo, disable.logging = T, filename = NULL)


## Compare t-stats
t_stat_plot <- function(expt, features, model1, model2){
  
  top_genes1 <- get(paste0("top_", features, "_", expt, "_", model1))
  top_genes2 <- get(paste0("top_", features, "_", expt, "_", model2))

  top_genes1$de <- top_genes1$adj.P.Val<0.05
  top_genes2$de <- top_genes2$adj.P.Val<0.05
    
  ## Model name
  name <- c("unadjusted" = "unadjusted", 
            "cluster_adjusted" = "cluster adjusted",
            "SV_adjusted" = "SV adjusted")
  
  ## Spearman's correlation coeff and pval
  res = cor.test(top_genes1$t, top_genes2$t, method="spearman", exact = T)
  res = data.frame(rho = res$estimate, rho_p = res$p.value)
  ## If p == 0, it's < 2.2e-16
  res$rho_p  <- ifelse(res$rho_p == 0, 2.2e-16, res$rho_p)

  ## Merge data
  t_stats <- data.frame(t1 = top_genes1$t, t2 = top_genes2$t, 
                        de1 = top_genes1$de, de2 = top_genes2$de)
  t_stats$de <- case_when(t_stats$de1 == TRUE & t_stats$de2 == TRUE ~ "sig. both",
                          t_stats$de1 == TRUE & t_stats$de2 == FALSE ~ "sig. 1 only",
                          t_stats$de1 == FALSE & t_stats$de2 == TRUE ~ "sig. 2 only",
                          t_stats$de1 == FALSE & t_stats$de2 == FALSE ~ "n.s.")
  
  ## Colors and transparency
  cols <- c("#8B8B00", "#CDCD00","#EED8AE", "gray80") 
  alphas <- c( 1, 1, 1,0.5)  
  names(cols) <- names(alphas) <- c("sig. both", "sig. 1 only", "sig. 2 only", "n.s.")
  
  plot <- ggplot(t_stats, aes(x = t1, y = t2, color=de, alpha=de)) +
    geom_point(size = 0.9) +
    scale_color_manual(values = cols, drop = T) + 
    scale_alpha_manual(values = alphas, drop=T) +
    labs(x = paste("t-stats", name[model1]), 
         y = paste("t-stats", name[model2]),
         subtitle = as.expression(bquote(~ rho  == .(signif(res$rho, 3)) ~ ", " ~ italic(.("p")) == .(signif(res$rho_p, 2)))), 
         color = "Differential expression",
         parse = T) +
    guides(alpha = 'none', color = guide_legend(override.aes = list(size=2))) + 
    theme_bw() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size=11),
          legend.title = element_text(size=12))
  
  return(plot)

}

p1 <- t_stat_plot("nic", "genes", "unadjusted", "cluster_adjusted")
p2 <- t_stat_plot("smo", "genes", "unadjusted", "cluster_adjusted")

plot_grid(p1, p2, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_genes_cluster_adjusted_vs_unadjusted.pdf", width = 9.5, height = 3.2)


## Compare gene rankings
rank_nic_cluster_adjusted <- rank(top_genes_nic_cluster_adjusted$P.Value)
rank_nic_unadjusted <- rank(top_genes_nic_unadjusted$P.Value)
plot(rank_nic_cluster_adjusted, rank_nic_unadjusted, pch=16, cex=.3)
abline(0,1,col=2)

rank_smo_cluster_adjusted <- rank(top_genes_smo_cluster_adjusted$P.Value)
rank_smo_unadjusted <- rank(top_genes_smo_unadjusted$P.Value)
plot(rank_smo_cluster_adjusted, rank_smo_unadjusted, pch=16, cex=.3)
abline(0,1,col=2)



# ------------------------------------------------------------------------------
#                     5.3 Surrogate Variable Analysis (SVA)
# ------------------------------------------------------------------------------
## SVs computed based on gene expression residuals
f <- ~ Group + Sex + plate + flowcell + rRNA_rate + overallMapRate + totalAssignedGene + ERCCsumLogErr + mitoRate

## SVs in PNE:
model = model.matrix(f, data = colData(rse_gene_nic))
# 1. Compute number of latent factors to estimate according to signif eigengenes (PCs)
n.sv = num.sv(assays(rse_gene_nic)$logcounts, model, method = "be", B = 100, seed = 123)
# 2. Estimate n.sv surrogate variables 
svatwostep <- twostepsva.build(assays(rse_gene_nic)$logcounts, model, n.sv)
SVs <- svatwostep$sv %>% as.data.frame()
colnames(SVs) <- paste0("SV", 1:svatwostep$n.sv)
## Add SVs to colData
colData(rse_gene_nic) <- cbind(colData(rse_gene_nic) , SVs)

## Add to tx, exon and jx rse objects (make sure samples are in same order as in gene rse)
colnames(rse_tx_nic) <- rse_tx_nic$SAMPLE_ID
rse_tx_nic <- rse_tx_nic[, rse_gene_nic$SAMPLE_ID]
colnames(rse_exon_nic) <- rse_exon_nic$SAMPLE_ID
rse_exon_nic <- rse_exon_nic[, rse_gene_nic$SAMPLE_ID]
colnames(rse_jx_nic) <- rse_jx_nic$SAMPLE_ID
rse_jx_nic <- rse_jx_nic[, rse_gene_nic$SAMPLE_ID]

identical(rse_gene_nic$SAMPLE_ID, rse_tx_nic$SAMPLE_ID)
identical(rse_gene_nic$SAMPLE_ID, rse_exon_nic$SAMPLE_ID)
identical(rse_gene_nic$SAMPLE_ID, rse_jx_nic$SAMPLE_ID)
# [1] TRUE

rse_tx_nic$cluster <- rse_exon_nic$cluster <- rse_jx_nic$cluster <- rse_gene_nic$cluster
colData(rse_tx_nic)[, paste0("SV", 1:4)] <- colData(rse_exon_nic)[, paste0("SV", 1:4)] <- colData(rse_jx_nic)[, paste0("SV", 1:4)] <- colData(rse_gene_nic)[, paste0("SV", 1:4)]

save(rse_gene_nic, file = paste0("processed-data/04_DEA/sensitivity_analyses/", 
                                 "rse_gene_brain_pups_nicotine_Kmeans_cluster_SV.Rdata"))
save(rse_tx_nic, file = paste0("processed-data/04_DEA/sensitivity_analyses/", 
                               "rse_tx_brain_pups_nicotine_Kmeans_cluster_SV.Rdata"))
save(rse_exon_nic, file = paste0("processed-data/04_DEA/sensitivity_analyses/", 
                                 "rse_exon_brain_pups_nicotine_Kmeans_cluster_SV.Rdata"))
save(rse_jx_nic, file = paste0("processed-data/04_DEA/sensitivity_analyses/", 
                                 "rse_jx_brain_pups_nicotine_Kmeans_cluster_SV.Rdata"))

## SVs in MSDP:
model = model.matrix(f, data = colData(rse_gene_smo))
n.sv = num.sv(assays(rse_gene_smo)$logcounts, model, method = "be", B = 100, seed = 123)
svatwostep <- twostepsva.build(assays(rse_gene_smo)$logcounts, model, n.sv)
SVs <- svatwostep$sv %>% as.data.frame()
colnames(SVs) <- paste0("SV", 1:svatwostep$n.sv)
colData(rse_gene_smo) <- cbind(colData(rse_gene_smo) , SVs)

colnames(rse_tx_smo) <- rse_tx_smo$SAMPLE_ID
rse_tx_smo <- rse_tx_smo[, rse_gene_smo$SAMPLE_ID]
colnames(rse_exon_smo) <- rse_exon_smo$SAMPLE_ID
rse_exon_smo <- rse_exon_smo[, rse_gene_smo$SAMPLE_ID]
colnames(rse_jx_smo) <- rse_jx_smo$SAMPLE_ID
rse_jx_smo <- rse_jx_smo[, rse_gene_smo$SAMPLE_ID]

identical(rse_gene_smo$SAMPLE_ID, rse_tx_smo$SAMPLE_ID)
identical(rse_gene_smo$SAMPLE_ID, rse_exon_smo$SAMPLE_ID)
identical(rse_gene_smo$SAMPLE_ID, rse_jx_smo$SAMPLE_ID)

rse_tx_smo$cluster <- rse_exon_smo$cluster <- rse_jx_smo$cluster <- rse_gene_smo$cluster
colData(rse_tx_smo)[, paste0("SV", 1:10)] <- colData(rse_exon_smo)[, paste0("SV", 1:10)] <- colData(rse_jx_smo)[, paste0("SV", 1:10)] <- colData(rse_gene_smo)[, paste0("SV", 1:10)]

save(rse_gene_smo, file = paste0("processed-data/04_DEA/sensitivity_analyses/", 
                                 "rse_gene_brain_pups_smoking_Kmeans_cluster_SV.Rdata"))
save(rse_tx_smo, file = paste0("processed-data/04_DEA/sensitivity_analyses/", 
                               "rse_tx_brain_pups_smoking_Kmeans_cluster_SV.Rdata"))
save(rse_exon_smo, file = paste0("processed-data/04_DEA/sensitivity_analyses/", 
                                 "rse_exon_brain_pups_smoking_Kmeans_cluster_SV.Rdata"))
save(rse_jx_smo, file = paste0("processed-data/04_DEA/sensitivity_analyses/", 
                               "rse_jx_brain_pups_smoking_Kmeans_cluster_SV.Rdata"))


## Compare SVs to clusters
plot_SV_vs_clusters <- function(expt, SV, color_var){
  
  rse <- get(paste0("rse_gene_", expt))
  data <- as.data.frame(colData(rse))
  
  shape_var = "Sex"
  
  plot <- ggplot(data = data, 
                   aes(x = cluster, y = get(SV))) + 
    geom_point(aes(color = get(color_var), shape = get(shape_var)), 
               size = 1.3, stroke = 1) + 
    geom_boxplot(fill = NA, outliers = F, box.colour = "gray20", width = 0.5, 
                 box.linewidth = 0.2, median.linewidth = 0.2, whisker.linewidth = 0.2) + 
    scale_color_manual(values = colors[[color_var]]) +
    scale_shape_manual(values = shapes[[shape_var]]) +
    theme_classic() + 
    labs(x = "K-mean cluster",
         y = SV,
         color = color_var, shape = shape_var) +
    theme(legend.position = "right",
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 10),
          legend.key.height = unit(4, units = "mm"),
          axis.title = element_text(size = 10)) 
  
  return(plot)
}

for(expt in c("nic", "smo")){
  i = 1
  plots <- list()
  for(color_var in c("Group", "Sex", "plate", "flowcell")){
    for(SV in paste0("SV", 1:4)){
      plots[[i]] <- plot_SV_vs_clusters(expt, SV, color_var)
      i = i + 1
    }
  }
  plot_grid(plotlist = plots, ncol = 4, align = "vh")
  ggsave(filename = paste0("plots/04_DEA/sensitivity_analyses/SV_vs_cluster_", 
         expt, ".pdf"), width = 17, height = 8)
}

# ------------------------------------------------------------------------------
## Run DGE adjusting for SVs 
dea_SV_adjusted <- function(expt, features){
  
  rse <- get(paste0("rse_", features, "_", expt))
  rse_norm <- calcNormFactors(rse, method = "TMM")
  f <- paste0("~ Group + Sex + plate + flowcell + rRNA_rate + overallMapRate +
            totalAssignedGene + ERCCsumLogErr + mitoRate + ", 
              paste(grep("SV", colnames(colData(rse)), value = T), collapse = " + "))
  
  model = model.matrix(as.formula(f), data = colData(rse))
  v = voom(rse_norm, design = model)
  fit = lmFit(v)
  eB = eBayes(fit)
  top = topTable(eB, coef = "GroupExperimental", 
                     p.value = 1, number = nrow(rse), sort.by="none")

  return(top)
}

top_genes_nic_SV_adjusted <- dea_SV_adjusted("nic", "gene")
nom_genes_nic_SV_adjusted <- subset(top_genes_nic_SV_adjusted, P.Value<0.05)
de_genes_nic_SV_adjusted <- subset(top_genes_nic_SV_adjusted, adj.P.Val<0.05)
de_genes_nic_SV_adjusted <- de_genes_nic_SV_adjusted[order(de_genes_nic_SV_adjusted$adj.P.Val, decreasing = F), ]

top_genes_smo_SV_adjusted <- dea_SV_adjusted("smo")
nom_genes_smo_SV_adjusted <- subset(top_genes_smo_SV_adjusted, P.Value<0.05)
de_genes_smo_SV_adjusted <- subset(top_genes_smo_SV_adjusted, adj.P.Val<0.05)
de_genes_smo_SV_adjusted <- de_genes_smo_SV_adjusted[order(de_genes_smo_SV_adjusted$adj.P.Val, decreasing = F), ]

## Check overlap with unadjusted and cluster-adjusted DEGs
DEG_nic <- list(
  "PNE SV adjusted DEGs" = de_genes_nic_SV_adjusted$ensemblID,
  "PNE cluster adjusted DEGs" = de_genes_nic_cluster_adjusted$ensemblID,
  "PNE unadjusted DEGs"= de_genes_nic_unadjusted$ensemblID)
venn.diagram(DEG_nic, disable.logging = T, filename = NULL)

DEG_smo <- list(
  "MSDP SV adjusted DEGs" = de_genes_smo_SV_adjusted$ensemblID,
  "MSDP cluster adjusted DEGs" = de_genes_smo_cluster_adjusted$ensemblID,
  "MSDP unadjusted DEGs"= de_genes_smo_unadjusted$ensemblID)
venn.diagram(DEG_smo, disable.logging = T, filename = NULL)

## % of DEGs found after SV adjustment at FDR<0.05
# PNE
length(intersect(de_genes_nic_SV_adjusted$ensemblID, de_genes_nic_unadjusted$ensemblID))/dim(de_genes_nic_unadjusted)[1] *100
# [1] 52.67327
# MSDP
length(intersect(de_genes_smo_SV_adjusted$ensemblID, de_genes_smo_unadjusted$ensemblID))/dim(de_genes_smo_unadjusted)[1] *100
# [1] 23.62545

## % of DEGs found after SV adjustment at p<0.05
# PNE
length(intersect(nom_genes_nic_SV_adjusted$ensemblID, de_genes_nic_unadjusted$ensemblID))/dim(de_genes_nic_unadjusted)[1] *100
# [1] 99.10891
# MSDP
length(intersect(nom_genes_smo_SV_adjusted$ensemblID, de_genes_smo_unadjusted$ensemblID))/dim(de_genes_smo_unadjusted)[1] *100
# [1] 52.26891


## Compare t-stats
p1 <- t_stat_plot("nic", "genes", "unadjusted", "SV_adjusted")
p2 <- t_stat_plot("smo", "genes", "unadjusted", "SV_adjusted")
plot_grid(p1, p2, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_genes_SV_adjusted_vs_unadjusted.pdf", width = 9.5, height = 3.2)

p3 <- t_stat_plot("nic", "genes", "cluster_adjusted", "SV_adjusted")
p4 <- t_stat_plot("smo", "genes", "cluster_adjusted", "SV_adjusted")
plot_grid(p3, p4, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_genes_SV_vs_cluster_adjusted.pdf", width = 9.5, height = 3.2)

## Compare gene rankings for PNE DGE
rank_nic_cluster_adjusted <- rank(top_genes_nic_cluster_adjusted$P.Value)
rank_nic_SV_adjusted <- rank(top_genes_nic_SV_adjusted$P.Value)
rank_nic_unadjusted <- rank(top_genes_nic_unadjusted$P.Value)
ranks_nic <- cbind("cluster_adjusted" = rank_nic_cluster_adjusted, 
                   "unadjusted" = rank_nic_unadjusted,
                   "SV_adjusted" = rank_nic_SV_adjusted)

p1 <- ggplot(ranks_nic, aes(x = unadjusted, y = cluster_adjusted)) +
  geom_point(size = 1, alpha = 0.2, color = "gray20") +
  theme_classic() +
  labs(x = "Unadjusted ranking",
       y = "Cluster-adjusted ranking") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  coord_cartesian(xlim = c(1, nrow(ranks_nic)), expand = F)

p2 <- ggplot(ranks_nic, aes(x = unadjusted, y = SV_adjusted)) +
  geom_point(size = 1, alpha = 0.2, color = "gray20") +
  theme_classic() +
  labs(x = "Unadjusted ranking",
       y = "SV-adjusted ranking") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  coord_cartesian(xlim = c(1, nrow(ranks_nic)), expand = F)

p3 <- ggplot(ranks_nic, aes(x = cluster_adjusted, y = SV_adjusted)) +
  geom_point(size = 1, alpha = 0.2, color = "gray20") +
  theme_classic() +
  labs(x = "Cluster-adjusted ranking",
       y = "SV-adjusted ranking") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  coord_cartesian(xlim = c(1, nrow(ranks_nic)), expand = F)

plot_grid(p1, p2 ,p3, ncol = 3)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/Rankings_DGE_PNE.pdf", width = 9.5, height = 3)


## Compare gene rankings for MSDP DGE
rank_smo_cluster_adjusted <- rank(top_genes_smo_cluster_adjusted$P.Value)
rank_smo_SV_adjusted <- rank(top_genes_smo_SV_adjusted$P.Value)
rank_smo_unadjusted <- rank(top_genes_smo_unadjusted$P.Value)
ranks_smo <- cbind("cluster_adjusted" = rank_smo_cluster_adjusted, 
                   "unadjusted" = rank_smo_unadjusted,
                   "SV_adjusted" = rank_smo_SV_adjusted)

p1 <- ggplot(ranks_smo, aes(x = unadjusted, y = cluster_adjusted)) +
  geom_point(size = 1, alpha = 0.2, color = "gray20") +
  theme_classic() +
  labs(x = "Unadjusted ranking",
       y = "Cluster-adjusted ranking") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  coord_cartesian(xlim = c(1, nrow(ranks_nic)), expand = F)

p2 <- ggplot(ranks_smo, aes(x = unadjusted, y = SV_adjusted)) +
  geom_point(size = 1, alpha = 0.2, color = "gray20") +
  theme_classic() +
  labs(x = "Unadjusted ranking",
       y = "SV-adjusted ranking") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  coord_cartesian(xlim = c(1, nrow(ranks_nic)), expand = F)

p3 <- ggplot(ranks_smo, aes(x = cluster_adjusted, y = SV_adjusted)) +
  geom_point(size = 1, alpha = 0.2, color = "gray20") +
  theme_classic() +
  labs(x = "Cluster-adjusted ranking",
       y = "SV-adjusted ranking") + 
  geom_abline(slope = 1, intercept = 0, color = "red") +
  coord_cartesian(xlim = c(1, nrow(ranks_nic)), expand = F)

plot_grid(p1, p2 ,p3, ncol = 3)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/Rankings_DGE_MSDP.pdf", width = 9.5, height = 3)


## Supp tables with results 
colnames(top_genes_nic_SV_adjusted)[1] <- "chr"
write.table(top_genes_nic_SV_adjusted, file = "processed-data/04_DEA/Gene_analysis/top_genes_brain_pup_nicotine_SV_adjusted.csv", row.names = FALSE, col.names = TRUE, sep = '\t')

colnames(top_genes_smo_SV_adjusted)[1] <- "chr"
write.table(top_genes_smo_SV_adjusted, file = "processed-data/04_DEA/Gene_analysis/top_genes_brain_pup_smoking_SV_adjusted.csv", row.names = FALSE, col.names = TRUE, sep = '\t')


# ------------------------------------------------------------------------------
## Run DTE adjusting for SVs
dte_SV_adjusted <- function(expt){
  
  rse <- get(paste0("rse_tx_", expt))
  
  ## Model matrix 
  f <- paste0("~ Group + Sex + plate + flowcell + rRNA_rate + overallMapRate +
            totalAssignedGene + ERCCsumLogErr + mitoRate + ", 
              paste(grep("SV", colnames(colData(rse)), value = T), collapse = " + "))
  
  model = model.matrix(as.formula(f), data = colData(rse))
  
  ## Fit linear model for each transcript
  fitTx = lmFit(assays(rse)$logcounts, design = model)
  
  ## Compute moderated F and t-statistics, and log-odds of DE
  eBTx = eBayes(fitTx)
  
  ## Select top-ranked transcripts for Group 
  top_tx = topTable(eBTx, coef = "GroupExperimental", p.value = 1, 
                    number = nrow(rse), sort.by = "none")
  
  ## Add relevant info 
  top_tx$Symbol <- rowData(rse)$gene_name
  top_tx$ensembl_id <- rowData(rse)$gene_id
  top_tx$transcript_id <- rowData(rse)$transcript_id
  top_tx$transcript_name <- rowData(rse)$transcript_name

  return(top_tx)
}

top_tx_nic_SV_adjusted <- dte_SV_adjusted("nic")
nom_tx_nic_SV_adjusted <- subset(top_tx_nic_SV_adjusted, P.Value<0.05)
de_tx_nic_SV_adjusted <- subset(top_tx_nic_SV_adjusted, adj.P.Val<0.05)
de_tx_nic_SV_adjusted <- de_tx_nic_SV_adjusted[order(de_tx_nic_SV_adjusted$adj.P.Val, decreasing = F), ]

top_tx_smo_SV_adjusted <- dte_SV_adjusted("smo")
nom_tx_smo_SV_adjusted <- subset(top_tx_smo_SV_adjusted, P.Value<0.05)
de_tx_smo_SV_adjusted <- subset(top_tx_smo_SV_adjusted, adj.P.Val<0.05)
de_tx_smo_SV_adjusted <- de_tx_smo_SV_adjusted[order(de_tx_smo_SV_adjusted$adj.P.Val, decreasing = F), ]

## Check overlap with unadjusted DETs
top_tx_nic_unadjusted <- get(load("processed-data/04_DEA/Tx_analysis/top_tx_nic.Rdata"))
top_tx_smo_unadjusted <- get(load("processed-data/04_DEA/Tx_analysis/top_tx_smo.Rdata"))
de_tx_nic_unadjusted <- subset(top_tx_nic_unadjusted, adj.P.Val<0.05)
de_tx_smo_unadjusted <- subset(top_tx_smo_unadjusted, adj.P.Val<0.05)
nom_tx_nic_unadjusted <- subset(top_tx_nic_unadjusted, P.Value<0.05)
nom_tx_smo_unadjusted <- subset(top_tx_smo_unadjusted, P.Value<0.05)

DET_nic <- list(
  "PNE SV adjusted DETs" = de_tx_nic_SV_adjusted$transcript_id,
  "PNE unadjusted DETs"= de_tx_nic_unadjusted$transcript_id)
venn.diagram(DET_nic, disable.logging = T, filename = NULL)

DET_smo <- list(
  "MSDP SV adjusted DETs" = de_tx_smo_SV_adjusted$transcript_id,
  "MSDP unadjusted DETs"= de_tx_smo_unadjusted$transcript_id)
venn.diagram(DET_smo, disable.logging = T, filename = NULL)

## % of DETs found after SV adjustment at FDR<0.05
# PNE
length(intersect(de_tx_nic_SV_adjusted$transcript_id, de_tx_nic_unadjusted$transcript_id))/dim(de_tx_nic_unadjusted)[1] *100
# [1] 59.05172
# MSDP
length(intersect(de_tx_smo_SV_adjusted$transcript_id, de_tx_smo_unadjusted$transcript_id))/dim(de_tx_smo_unadjusted)[1] *100
# [1] 9.164819

## % of DETs found after SV adjustment at p<0.05
# PNE
length(intersect(nom_tx_nic_SV_adjusted$transcript_id, de_tx_nic_unadjusted$transcript_id))/dim(de_tx_nic_unadjusted)[1] *100
# [1] 100
# MSDP
length(intersect(nom_tx_smo_SV_adjusted$transcript_id, de_tx_smo_unadjusted$transcript_id))/dim(de_tx_smo_unadjusted)[1] *100
# [1] 56.56566

## Compare t-stats
p1 <- t_stat_plot("nic", "tx", "unadjusted", "SV_adjusted")
p2 <- t_stat_plot("smo", "tx", "unadjusted", "SV_adjusted")
plot_grid(p1, p2, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_txs_SV_adjusted_vs_unadjusted.pdf", width = 9.5, height = 3.2)


# ------------------------------------------------------------------------------
## Run DEE adjusting for SVs
top_exons_nic_SV_adjusted <- dea_SV_adjusted("nic", "exon")
nom_exon_nic_SV_adjusted <- subset(top_exons_nic_SV_adjusted, P.Value<0.05 & abs(logFC)>0.25)
de_exon_nic_SV_adjusted <- subset(top_exons_nic_SV_adjusted, adj.P.Val<0.05 & abs(logFC)>0.25)
de_exon_nic_SV_adjusted <- de_exon_nic_SV_adjusted[order(de_exon_nic_SV_adjusted$adj.P.Val, decreasing = F), ]

top_exons_smo_SV_adjusted <- dea_SV_adjusted("smo", "exon")
nom_exon_smo_SV_adjusted <- subset(top_exons_smo_SV_adjusted, P.Value<0.05 & abs(logFC)>0.25)
de_exon_smo_SV_adjusted <- subset(top_exons_smo_SV_adjusted, adj.P.Val<0.05 & abs(logFC)>0.25)
de_exon_smo_SV_adjusted <- de_exon_smo_SV_adjusted[order(de_exon_smo_SV_adjusted$adj.P.Val, decreasing = F), ]

## Check overlap with unadjusted DEEs
top_exons_nic_unadjusted <- get(load("processed-data/04_DEA/Exon_analysis/top_exons_nic.Rdata"))
top_exons_smo_unadjusted <- get(load("processed-data/04_DEA/Exon_analysis/top_exons_smo.Rdata"))
de_exon_nic_unadjusted <- subset(top_exons_nic_unadjusted, adj.P.Val<0.05 & abs(logFC)>0.25)
de_exon_smo_unadjusted <- subset(top_exons_smo_unadjusted, adj.P.Val<0.05 & abs(logFC)>0.25)
nom_exon_nic_unadjusted <- subset(top_exons_nic_unadjusted, P.Value<0.05 & abs(logFC)>0.25)
nom_exon_smo_unadjusted <- subset(top_exons_smo_unadjusted, P.Value<0.05 & abs(logFC)>0.25)

DEE_nic <- list(
  "PNE SV adjusted DEEs" = de_exon_nic_SV_adjusted$exon_gencodeID,
  "PNE unadjusted DEEs"= de_exon_nic_unadjusted$exon_gencodeID)
venn.diagram(DEE_nic, disable.logging = T, filename = NULL)

DEE_smo <- list(
  "MSDP SV adjusted DEEs" = de_exon_smo_SV_adjusted$exon_gencodeID,
  "MSDP unadjusted DEEs"= de_exon_smo_unadjusted$exon_gencodeID)
venn.diagram(DEE_smo, disable.logging = T, filename = NULL)

## % of DEEs found after SV adjustment at FDR<0.05
# PNE
length(intersect(de_exon_nic_SV_adjusted$exon_gencodeID, de_exon_nic_unadjusted$exon_gencodeID))/dim(de_exon_nic_unadjusted)[1] *100
# [1] 38.02691
# MSDP
length(intersect(de_exon_smo_SV_adjusted$exon_gencodeID, de_exon_smo_unadjusted$exon_gencodeID))/dim(de_exon_smo_unadjusted)[1] *100
# [1] 5.448772

## % of DEEs found after SV adjustment at p<0.05
# PNE
length(intersect(nom_exon_nic_SV_adjusted$exon_gencodeID, de_exon_nic_unadjusted$exon_gencodeID))/dim(de_exon_nic_unadjusted)[1] *100
# [1] 76.23318
# MSDP
length(intersect(nom_exon_smo_SV_adjusted$exon_gencodeID, de_exon_smo_unadjusted$exon_gencodeID))/dim(de_exon_smo_unadjusted)[1] *100
# [1] 28.26341

## Compare t-stats
p1 <- t_stat_plot("nic", "exons", "unadjusted", "SV_adjusted")
p2 <- t_stat_plot("smo", "exons", "unadjusted", "SV_adjusted")
plot_grid(p1, p2, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_exons_SV_adjusted_vs_unadjusted.pdf", width = 9.5, height = 3.2)


# ------------------------------------------------------------------------------
## Run DJE adjusting for SVs
top_jxs_nic_SV_adjusted <- dea_SV_adjusted("nic", "jx")
top_jxs_nic_SV_adjusted$jxn_ID <- rownames(top_jxs_nic_SV_adjusted)
nom_jx_nic_SV_adjusted <- subset(top_jxs_nic_SV_adjusted, P.Value<0.05 & abs(logFC)>0.25)
de_jx_nic_SV_adjusted <- subset(top_jxs_nic_SV_adjusted, adj.P.Val<0.05 & abs(logFC)>0.25)
de_jx_nic_SV_adjusted <- de_jx_nic_SV_adjusted[order(de_jx_nic_SV_adjusted$adj.P.Val, decreasing = F), ]

top_jxs_smo_SV_adjusted <- dea_SV_adjusted("smo", "jx")
top_jxs_smo_SV_adjusted$jxn_ID <- rownames(top_jxs_smo_SV_adjusted)
nom_jx_smo_SV_adjusted <- subset(top_jxs_smo_SV_adjusted, P.Value<0.05 & abs(logFC)>0.25)
de_jx_smo_SV_adjusted <- subset(top_jxs_smo_SV_adjusted, adj.P.Val<0.05 & abs(logFC)>0.25)
de_jx_smo_SV_adjusted <- de_jx_smo_SV_adjusted[order(de_jx_smo_SV_adjusted$adj.P.Val, decreasing = F), ]

## Check overlap with unadjusted DEJs
top_jxs_nic_unadjusted <- get(load("processed-data/04_DEA/Jx_analysis/top_jxns_nic.Rdata"))
top_jxs_smo_unadjusted <- get(load("processed-data/04_DEA/Jx_analysis/top_jxns_smo.Rdata"))
top_jxs_nic_unadjusted$jxn_ID <- rownames(top_jxs_nic_unadjusted)
top_jxs_smo_unadjusted$jxn_ID <- rownames(top_jxs_smo_unadjusted)
de_jx_nic_unadjusted <- subset(top_jxs_nic_unadjusted, adj.P.Val<0.05)
de_jx_smo_unadjusted <- subset(top_jxs_smo_unadjusted, adj.P.Val<0.05)
nom_jx_nic_unadjusted <- subset(top_jxs_nic_unadjusted, P.Value<0.05)
nom_jx_smo_unadjusted <- subset(top_jxs_smo_unadjusted, P.Value<0.05)

DEJ_nic <- list(
  "PNE SV adjusted DEJs" = de_jx_nic_SV_adjusted$jxn_ID,
  "PNE unadjusted DEJs"= de_jx_nic_unadjusted$jxn_ID)
venn.diagram(DEJ_nic, disable.logging = T, filename = NULL)

DEJ_smo <- list(
  "MSDP SV adjusted DEJs" = de_jx_smo_SV_adjusted$jxn_ID,
  "MSDP unadjusted DEJs"= de_jx_smo_unadjusted$jxn_ID)
venn.diagram(DEJ_smo, disable.logging = T, filename = NULL)


## % of DEJs found after SV adjustment at FDR<0.05
# PNE
length(intersect(de_jx_nic_SV_adjusted$jxn_ID, de_jx_nic_unadjusted$jxn_ID))/dim(de_jx_nic_unadjusted)[1] *100
# [1] 34.03141
# MSDP
length(intersect(de_jx_smo_SV_adjusted$jxn_ID, de_jx_smo_unadjusted$jxn_ID))/dim(de_jx_smo_unadjusted)[1] *100
# [1] 0.7020117

## % of DEJs found after SV adjustment at p<0.05
# PNE
length(intersect(nom_jx_nic_SV_adjusted$jxn_ID, de_jx_nic_unadjusted$jxn_ID))/dim(de_jx_nic_unadjusted)[1] *100
# [1] 61.7801
# MSDP
length(intersect(nom_jx_smo_SV_adjusted$jxn_ID, de_jx_smo_unadjusted$jxn_ID))/dim(de_jx_smo_unadjusted)[1] *100
# [1] jxn_ID

## Compare t-stats
p1 <- t_stat_plot("nic", "jxs", "unadjusted", "SV_adjusted")
p2 <- t_stat_plot("smo", "jxs", "unadjusted", "SV_adjusted")
plot_grid(p1, p2, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_jxs_SV_adjusted_vs_unadjusted.pdf", width = 9.5, height = 3.2)







session_info()
# setting  value
# version  R version 4.5.2 (2025-10-31)
# os       macOS Sequoia 15.3.2
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       Europe/London
# date     2026-06-21
# rstudio  2026.01.0+392 Apple Blossom (desktop)
# pandoc   NA
# quarto   1.8.25 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-8     2024-09-12 [1] CRAN (R 4.5.0)
# annotate               1.88.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# AnnotationDbi          1.72.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# Biobase              * 2.70.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# BiocFileCache          3.0.0     2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# BiocGenerics         * 0.56.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# BiocManager            1.30.27   2025-11-14 [1] CRAN (R 4.5.2)
# BiocParallel         * 1.44.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# biomaRt                2.66.1    2026-02-12 [1] https://bioc-release.r-universe.dev (R 4.5.2)
# biomartr             * 1.0.7     2023-12-02 [1] CRAN (R 4.5.0)
# Biostrings             2.78.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# bit                    4.6.0     2025-03-06 [1] CRAN (R 4.5.0)
# bit64                  4.8.0     2026-04-21 [1] CRAN (R 4.5.2)
# bitops                 1.0-9     2024-10-03 [1] CRAN (R 4.5.0)
# blob                   1.3.0     2026-01-14 [1] CRAN (R 4.5.2)
# cachem                 1.1.0     2024-05-16 [1] CRAN (R 4.5.0)
# cli                    3.6.6     2026-04-09 [1] CRAN (R 4.5.2)
# codetools              0.2-20    2024-03-31 [1] CRAN (R 4.5.2)
# cowplot              * 1.2.0     2025-07-07 [1] CRAN (R 4.5.0)
# crayon                 1.5.3     2024-06-20 [1] CRAN (R 4.5.0)
# curl                   7.1.0     2026-04-22 [1] CRAN (R 4.5.2)
# data.table             1.18.4    2026-05-06 [1] CRAN (R 4.5.2)
# DBI                    1.3.0     2026-02-25 [1] CRAN (R 4.5.2)
# dbplyr                 2.5.2     2026-02-13 [1] CRAN (R 4.5.2)
# DelayedArray           0.36.1    2026-03-31 [1] https://bioc-release.r-universe.dev (R 4.5.3)
# dplyr                * 1.2.1     2026-04-03 [1] CRAN (R 4.5.2)
# edgeR                * 4.8.2     2025-12-23 [1] https://bioc-release.r-universe.dev (R 4.5.2)
# farver                 2.1.2     2024-05-13 [1] CRAN (R 4.5.0)
# fastmap                1.2.0     2024-05-15 [1] CRAN (R 4.5.0)
# filelock               1.0.3     2023-12-11 [1] CRAN (R 4.5.0)
# formatR                1.14      2023-01-17 [1] CRAN (R 4.5.0)
# futile.logger        * 1.4.9     2025-12-29 [1] CRAN (R 4.5.2)
# futile.options         1.0.1     2018-04-20 [1] CRAN (R 4.5.0)
# genefilter           * 1.92.0    2025-10-29 [1] https://bioc-release.r-universe.dev (R 4.5.2)
# generics             * 0.1.4     2025-05-09 [1] CRAN (R 4.5.0)
# GenomicRanges        * 1.62.1    2025-12-08 [1] Bioconductor 3.22 (R 4.5.2)
# ggplot2              * 4.0.3     2026-04-22 [1] CRAN (R 4.5.2)
# ggrepel              * 0.9.8     2026-03-17 [1] CRAN (R 4.5.2)
# glue                   1.8.1     2026-04-17 [1] CRAN (R 4.5.2)
# gtable                 0.3.6     2024-10-25 [1] CRAN (R 4.5.0)
# here                 * 1.0.2     2025-09-15 [1] CRAN (R 4.5.0)
# hms                    1.1.4     2025-10-17 [1] CRAN (R 4.5.0)
# httr                   1.4.8     2026-02-13 [1] CRAN (R 4.5.2)
# httr2                  1.2.2     2025-12-08 [1] CRAN (R 4.5.2)
# IRanges              * 2.44.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# KEGGREST               1.50.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.5.0)
# lambda.r               1.2.4     2019-09-18 [1] CRAN (R 4.5.0)
# lattice                0.22-9    2026-02-09 [1] CRAN (R 4.5.2)
# lifecycle              1.0.5     2026-01-08 [1] CRAN (R 4.5.2)
# limma                * 3.66.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# locfit                 1.5-9.12  2025-03-05 [1] CRAN (R 4.5.0)
# magrittr               2.0.5     2026-04-04 [1] CRAN (R 4.5.2)
# Matrix                 1.7-5     2026-03-21 [1] CRAN (R 4.5.2)
# MatrixGenerics       * 1.22.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# matrixStats          * 1.5.0     2025-01-07 [1] CRAN (R 4.5.0)
# memoise                2.0.1     2021-11-26 [1] CRAN (R 4.5.0)
# mgcv                 * 1.9-4     2025-11-07 [1] CRAN (R 4.5.0)
# nlme                 * 3.1-169   2026-03-27 [1] CRAN (R 4.5.2)
# otel                   0.2.0     2025-08-29 [1] CRAN (R 4.5.0)
# pheatmap             * 1.0.13    2025-06-05 [1] CRAN (R 4.5.0)
# pillar                 1.11.1    2025-09-17 [1] CRAN (R 4.5.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.5.0)
# png                    0.1-9     2026-03-15 [1] CRAN (R 4.5.2)
# prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.5.0)
# progress               1.2.3     2023-12-06 [1] CRAN (R 4.5.0)
# purrr                  1.2.2     2026-04-10 [1] CRAN (R 4.5.2)
# R6                     2.6.1     2025-02-15 [1] CRAN (R 4.5.0)
# ragg                   1.5.2     2026-03-23 [1] CRAN (R 4.5.2)
# rappdirs               0.3.4     2026-01-17 [1] CRAN (R 4.5.2)
# RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.5.0)
# Rcpp                   1.1.1-1.1 2026-04-24 [1] CRAN (R 4.5.2)
# RCurl                  1.98-1.18 2026-03-21 [1] CRAN (R 4.5.2)
# rlang                * 1.2.0     2026-04-06 [1] CRAN (R 4.5.2)
# rprojroot              2.1.1     2025-08-26 [1] CRAN (R 4.5.0)
# rsconnect              1.8.0     2026-04-10 [1] CRAN (R 4.5.2)
# RSQLite                3.52.0    2026-05-10 [1] CRAN (R 4.5.2)
# rstudioapi             0.18.0    2026-01-16 [1] CRAN (R 4.5.2)
# S4Arrays               1.10.1    2025-12-08 [1] Bioconductor 3.22 (R 4.5.2)
# S4Vectors            * 0.48.1    2026-04-04 [1] https://bioc-release.r-universe.dev (R 4.5.3)
# S7                     0.2.2     2026-04-22 [1] CRAN (R 4.5.2)
# scales                 1.4.0     2025-04-24 [1] CRAN (R 4.5.0)
# Seqinfo              * 1.0.0     2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# sessioninfo          * 1.2.3     2025-02-05 [1] CRAN (R 4.5.0)
# SparseArray            1.10.10   2026-03-30 [1] https://bioc-release.r-universe.dev (R 4.5.3)
# statmod                1.5.1     2025-10-09 [1] CRAN (R 4.5.0)
# stringi                1.8.7     2025-03-27 [1] CRAN (R 4.5.0)
# stringr                1.6.0     2025-11-04 [1] CRAN (R 4.5.0)
# SummarizedExperiment * 1.40.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# survival               3.8-6     2026-01-16 [1] CRAN (R 4.5.2)
# sva                  * 3.58.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# systemfonts            1.3.2     2026-03-05 [1] CRAN (R 4.5.2)
# textshaping            1.0.5     2026-03-06 [1] CRAN (R 4.5.2)
# tibble                 3.3.1     2026-01-11 [1] CRAN (R 4.5.2)
# tidyr                * 1.3.2     2025-12-19 [1] CRAN (R 4.5.2)
# tidyselect             1.2.1     2024-03-11 [1] CRAN (R 4.5.0)
# vctrs                  0.7.3     2026-04-11 [1] CRAN (R 4.5.2)
# VennDiagram          * 1.8.2     2026-01-11 [1] CRAN (R 4.5.2)
# withr                  3.0.2     2024-10-28 [1] CRAN (R 4.5.0)
# XML                    3.99-0.23 2026-03-20 [1] CRAN (R 4.5.2)
# xml2                   1.5.2     2026-01-17 [1] CRAN (R 4.5.2)
# xtable                 1.8-8     2026-02-22 [1] CRAN (R 4.5.2)
# XVector                0.50.0    2025-10-29 [1] Bioconductor 3.22 (R 4.5.1)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library
# * ── Packages attached to the search path.
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
