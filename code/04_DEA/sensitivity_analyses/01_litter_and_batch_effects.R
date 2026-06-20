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
t_stat_plot <- function(expt, model1, model2){
  
  top_genes1 <- get(paste0("top_genes_", expt, "_", model1))
  top_genes2 <- get(paste0("top_genes_", expt, "_", model2))

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
                          t_stats$de1 == TRUE & t_stats$de2 == FALSE ~ "sig. unadj only",
                          t_stats$de1 == FALSE & t_stats$de2 == TRUE ~ "sig. adj only",
                          t_stats$de1 == FALSE & t_stats$de2 == FALSE ~ "n.s.")
  
  ## Colors and transparency
  cols <- c("#8B8B00", "#CDCD00","#EED8AE", "gray80") 
  alphas <- c( 1, 1, 1,0.5)  
  names(cols) <- names(alphas) <- c("sig. both", "sig. unadj only", "sig. adj only", "n.s.")
  
  plot <- ggplot(t_stats, aes(x = t1, y = t2, color=de, alpha=de)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = cols, drop = T) + 
    scale_alpha_manual(values = alphas, drop=T) +
    labs(x = paste("t-stats", name[model1]), 
         y = paste("t-stats", name[model2]),
         subtitle = as.expression(bquote(~ rho  == .(signif(res$rho, 2)) ~ ", " ~ italic(.("p")) == .(signif(res$rho_p, 2)))), 
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

p1 <- t_stat_plot("nic", "unadjusted", "cluster_adjusted")
p2 <- t_stat_plot("smo", "unadjusted", "cluster_adjusted")

plot_grid(p1, p2, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_cluster_adjusted_vs_unadjusted.pdf", width = 9.5, height = 3.2)


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

## For PNE
## Add to tx, exon and jx rse objects
## Make sure samples are in same order as in gene rse
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

## Add cluster and SVs
rse_tx_nic$cluster <- rse_exon_nic$cluster <- rse_jx_nic$cluster <- rse_gene_nic$cluster
colData(rse_tx_nic)[, paste0("SV", 1:4)] <- colData(rse_exon_nic)[, paste0("SV", 1:4)] <- colData(rse_jx_nic)[, paste0("SV", 1:4)] <- colData(rse_gene_nic)[, paste0("SV", 1:4)]

## SVs in MSDP:
model = model.matrix(f, data = colData(rse_gene_smo))
n.sv = num.sv(assays(rse_gene_smo)$logcounts, model, method = "be", B = 100, seed = 123)
svatwostep <- twostepsva.build(assays(rse_gene_smo)$logcounts, model, n.sv)
SVs <- svatwostep$sv %>% as.data.frame()
colnames(SVs) <- paste0("SV", 1:svatwostep$n.sv)
colData(rse_gene_smo) <- cbind(colData(rse_gene_smo) , SVs)

## For MSDP
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
## Run DGE adjusting for SVs HEREEEEEEEE
dge_SV_adjusted <- function(expt){
  
  rse <- get(paste0("rse_gene_", expt))
  rse_norm <- calcNormFactors(rse, method = "TMM")
  f <- paste0("~ Group + Sex + plate + flowcell + rRNA_rate + overallMapRate +
            totalAssignedGene + ERCCsumLogErr + mitoRate + ", 
              paste(grep("SV", colnames(colData(rse)), value = T), collapse = " + "))
  
  model = model.matrix(as.formula(f), data = colData(rse))
  vGene = voom(rse_norm, design = model)
  fitGene = lmFit(vGene)
  eBGene = eBayes(fitGene)
  top_genes = topTable(eBGene, coef = "GroupExperimental", 
                       p.value = 1, number = nrow(rse), sort.by="none")

  return(top_genes)
}

top_genes_nic_SV_adjusted <- dge_SV_adjusted("nic")
nom_genes_nic_SV_adjusted <- subset(top_genes_nic_SV_adjusted, P.Value<0.05)
de_genes_nic_SV_adjusted <- subset(top_genes_nic_SV_adjusted, adj.P.Val<0.05)
de_genes_nic_SV_adjusted <- de_genes_nic_SV_adjusted[order(de_genes_nic_SV_adjusted$adj.P.Val, decreasing = F), ]

top_genes_smo_SV_adjusted <- dge_SV_adjusted("smo")
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
p1 <- t_stat_plot("nic", "unadjusted", "SV_adjusted")
p2 <- t_stat_plot("smo", "unadjusted", "SV_adjusted")
plot_grid(p1, p2, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_SV_adjusted_vs_unadjusted.pdf", width = 9.5, height = 3.2)

p3 <- t_stat_plot("nic", "cluster_adjusted", "SV_adjusted")
p4 <- t_stat_plot("smo", "cluster_adjusted", "SV_adjusted")
plot_grid(p3, p4, ncol = 2)
ggsave(filename = "plots/04_DEA/sensitivity_analyses/t_stats_SV_vs_cluster_adjusted.pdf", width = 9.5, height = 3.2)

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


# ------------------------------------------------------------------------------
## Run DTE adjusting for SVs
dge_SV_adjusted <- function(expt){
  
  rse <- get(paste0("rse_tx_", expt))
  
  ## Model matrix using formula for the fitted model
  f <- ~ Group + Sex + plate + flowcell + rRNA_rate + totalAssignedGene + ERCCsumLogErr + overallMapRate + mitoRate
  model = model.matrix(f, data = colData(rse))
  
  ## Fit linear model for each transcript
  fitTx = lmFit(assays(rse)$logcounts, design = model)
  
  ## Compute moderated F and t-statistics, and log-odds of DE
  eBTx = eBayes(fitTx)
  
  ## Plot average log expression vs logFC
  limma::plotMA(eBTx, coef = "GroupExperimental", xlab = "Mean of normalized counts", 
                ylab="logFC")
  
  ## Plot -log(p-value) vs logFC
  volcanoplot(eBTx, coef = "GroupExperimental")
  
  ## Select top-ranked transcripts for Group 
  top_tx = topTable(eBTx, coef="GroupExperimental", p.value = 1, number=nrow(RSE), sort.by="none")
  ## Histogram of adjusted p values
  hist(top_tx$adj.P.Val, xlab="FDR", main="")
  
  
  
  
  
  rse_norm <- calcNormFactors(rse, method = "TMM")
  f <- paste0("~ Group + Sex + plate + flowcell + rRNA_rate + overallMapRate +
            totalAssignedGene + ERCCsumLogErr + mitoRate + ", 
              paste(grep("SV", colnames(colData(rse)), value = T), collapse = " + "))
  
  model = model.matrix(as.formula(f), data = colData(rse))
  vGene = voom(rse_norm, design = model)
  fitGene = lmFit(vGene)
  eBGene = eBayes(fitGene)
  top_genes = topTable(eBGene, coef = "GroupExperimental", 
                       p.value = 1, number = nrow(rse), sort.by="none")
  
  return(top_genes)
}