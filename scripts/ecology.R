######################## ECOLOGICAL ANALYSES ##########################
#=======================================================================

# The population structure was determined as the covariance among populations, an (\(\Omega\)) 
# matrix derived from putatively neutral loci (four-fold degenerate sites, introns, and intergenic regions),
# similar to \protect\citep{gautier2015}, where genetic distance between two populations is assumed to
# co-vary based on their shared history, fillowed by a PCA. PCA axes are the single value decomposition (SVD)
# of pairwise covariance.

#================================================================================
### Ecology vs CNV per family [with MCA] ------------------
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(data.table)
library(vegan)
## genetic diversity from long reads --------
svd_df<-read.csv("data/DefenseNstress/AT/genDist/AT_gen_div_covar_svd.csv")
# ** Run multiple regression on cnv vs env on the def & stress gene families ***
cnv_out<-read.csv("data/DefenseNstress/AT/cnv_within/cnv_percentage_per_family2.csv")
lr_met<-fread("data/ncbi_LR_136dataset_cleaned_w_sampleInfo_with_qgraph_genetic_group2.1.csv",h=T)
env_comb<-read.csv("data/DefenseNstress/AT/eco/AT_chelsa_bioclim_for_all_pop_Apr22_25.csv")
# ** drop N.American and mixed-up Japanese populations
drop_samples<-c("GCA_036926965.1","GCA_036937515.1","GCA_036940625.1","GCA_036926925.1","GCA_020911765.2","GCA_946499705.1","GCA_036926975.1","GCA_946407795.1","GCA_949796485.1")# last ones are Ip-Lor-16 & IP-Tri-0 with assembly issues
drop_index<-match(drop_samples,lr_met$gc_name)
lr_met<-lr_met[-drop_index,]
og_vs_ass<-read.table("data/DefenseNstress/AT/cnv_within/geneFamily_vs_assembly_cnv2.2.txt",h=T)
as_nams<-substr(colnames(og_vs_ass)[-1],1,15)
colnames(og_vs_ass)<-c("family_id",as_nams)
rownames(og_vs_ass)<-og_vs_ass$family_id
og_vs_ass<-t(og_vs_ass[,-1])
ass_index<-match(lr_met$gc_name,as_nams)
og_vs_ass<-og_vs_ass[ass_index,]

ev_index<-match(lr_met$gc_name,env_comb$accession)
pca_env<-env_comb[ev_index,]

pca_input <- as.data.frame(og_vs_ass)
nm<-lr_met$sample_name.1
rownames(pca_input)<-nm

## pure MCA first 
pca_matrix <- as.data.frame(pca_input)
pca_matrix[] <- lapply(pca_matrix, as.factor)

# MCA 
mca_result <- MCA(pca_matrix, graph = F)
# Extract variance explained
explained_variance <- mca_result$eig[, 2]  # Percentage of variance per component
pc1_var <- round(explained_variance[1], 1)
pc2_var <- round(explained_variance[2], 1)

# Extract MCA coordinates
mca_df <- as.data.frame(mca_result$ind$coord)
mca_df$Assembly <- rownames(mca_df)
colnames(mca_df)<-c(paste0("Dim",1:5),"Assembly")

# add groupings
ggroup<-lr_met[match(mca_df$Assembly,lr_met$sample_name.1),]
mca_df$group<-ggroup$geo_group
mca_df$relict<-ggroup$relict
sams<-stringr::str_split_fixed(mca_df$Assembly,"\\ ",n=2)[,1]
lr_clusters<-svd_df[match(sams,svd_df$Sample),]
mca_df$LR_cluster<-factor(lr_clusters$LR_cluster)


# Visualize MCA with variance in axis labels
p_uncorrected<-
  ggplot(mca_df, aes(x = Dim1, y = Dim2, label = Assembly)) +
  geom_point(
    size = 3,
    aes(
      color = group,      # Color by group
      shape = LR_cluster,           # Shape by cluster
      fill = relict            # Fill by relict status 
    ),
    stroke = 1                 # Border thickness
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24),
    name = "Genetic Cluster"
  ) +
  # Color scale for clusters
  scale_color_brewer(
    palette = "Set1",
    name = "Geographic Group"
  ) +
  # Fill scale for relict status
  scale_fill_manual(
    values = c("non-relict" = "white", "relict" = "gray50"),
    name = "Relict Status"
  ) +
  # Text labels for assemblies
  geom_text_repel(
    aes(label = Assembly),
    size = 3,
    color = "black",
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  theme_minimal() +
  labs(title = "A.thaliana: MCA of LR-Assemblies Based on CNV in Def.Stress Gene Families",
       x = paste0("Dimension 1 (", pc1_var, "% variance)"),
       y = paste0("Dimension 2 (", pc2_var, "% variance)"))+
  guides(
    color = guide_legend(override.aes = list(shape = 15)),  # Squares for color legend
    shape = guide_legend(override.aes = list(fill = "white")),  # White fill for shape legend
    fill = guide_legend(override.aes = list(shape = 21))  # Circles for fill legend
  )

print(p_uncorrected)


## correct for IBD ---
Variables<-pca_env
ev1<-Variables[-c(1:8)] # no need to use PCs here
ccr<-cor(ev1,use="pairwise.complete.obs")
# corrplot::corrplot(ccr,type="lower")
diag(ccr)<-0
thr<-max(abs(ccr))
while(thr>0.9){
  ccr<-cor(ev1,use="pairwise.complete.obs")
  diag(ccr)<-0
  thr<-max(abs(ccr))
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)
evars<-lyrs
eq<-paste0(evars,collapse = "+")

## Remove IBD ----
mca_res<-mca_result
mca_coord <- data.frame(mca_res$ind$coord)
# Fit environmental variables
envfit_res <- envfit(mca_coord, Variables[,evars], permutations = 999)
arrow_df <- as.data.frame(envfit_res$vectors$arrows)
arrow_df$label <- rownames(arrow_df)
arrow_df$p_value <- envfit_res$vectors$pvals
write.csv(arrow_df,"data/DefenseNstress/AT/eco/AT_MCA_model_sifinificance_before.csv")

pop_structure <- as.data.frame(lr_clusters[,c("PC1","PC2")])
# Regress out geography from MCA dimensions
mca_corrected <- apply(mca_res$ind$coord, 2, function(y) {
  residuals(lm(y ~ ., data = pop_structure))
})
mca_corrected<-data.frame(mca_corrected)
# fit envfit on geography-corrected coordinates
envfit_res <- envfit(mca_corrected, Variables[,evars], permutations = 999)

mca_corrected$Assembly <- rownames(mca_coord)
ggroup<-lr_met[match(mca_corrected$Assembly,lr_met$sample_name.1),]
mca_corrected$group<-ggroup$geo_group
mca_corrected$relict<-ggroup$relict
mca_corrected$LR_cluster<-factor(lr_clusters$LR_cluster)
# Get arrow data as a dataframe
arrow_df <- as.data.frame(envfit_res$vectors$arrows)
arrow_df$label <- rownames(arrow_df)
arrow_df$p_value <- envfit_res$vectors$pvals
write.csv(arrow_df,"data/DefenseNstress/AT/eco/AT_MCA_model_sifinificance_after.csv")
arrow_df<-arrow_df[arrow_df$p_value<0.05,]

# If pc1_var and pc2_var are in decimal form (e.g., 0.12 for 12%), convert:
pc1_label <- paste0("Dimension 1 (", pc1_var , "% variance)")
pc2_label <- paste0("Dimension 2 (", pc2_var , "% variance)")


mca_corrected$LR_cluster <- as.factor(mca_corrected$LR_cluster)
mca_corrected$group <- as.factor(mca_corrected$group)
mca_corrected$relict <- as.factor(mca_corrected$relict)

# Plot
p_corrected<-ggplot(mca_corrected, aes(x = Dim.1, y = Dim.2, label = Assembly)) +
  geom_point(
    size = 3,
    aes(
      color = group,      # Color by group
      shape = LR_cluster,           # Shape by cluster
      fill = relict            # Fill by relict status 
    ),
    stroke = 1                 # Border thickness
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 24),
    name = "Genetic Cluster"
  ) +
  # Color scale for clusters
  scale_color_brewer(
    palette = "Set1",
    name = "Geographic Group"
  ) +
  # Fill scale for relict status
  scale_fill_manual(
    values = c("non-relict" = "white", "relict" = "gray50"),
    name = "Relict Status"
  ) +
  # Text labels for assemblies
  geom_text_repel(
    aes(label = Assembly),
    size = 3,
    color = "black",
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  # Environment arrows
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black",
    inherit.aes = FALSE
  ) +
  # Arrow labels
  geom_text(
    data = arrow_df,
    aes(x = Dim.1, y = Dim.2, label = label),
    color = "red",
    hjust = 1,
    vjust = 1,
    inherit.aes = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "A. thaliana: MCA EAA Based on CNV in Def./Stress Gene Families",
    x = pc1_label,
    y = pc2_label
  ) +
  
  # Legend adjustments
  guides(
    color = guide_legend(override.aes = list(shape = 15)),  # Squares for color legend
    shape = guide_legend(override.aes = list(fill = "white")),  # White fill for shape legend
    fill = guide_legend(override.aes = list(shape = 21))  # Circles for fill legend
  )

print(p_corrected)


# Identify most influential gene families using geography-corrected MCA ----
mca_var <- mca_result$var
var_coord <- as.data.frame(mca_var$coord)
colnames(var_coord)<-paste0("Dim.",1:5)
var_coord$category <- rownames(var_coord)
# approximate influence: distance from origin
var_coord$contrib <- rowSums(var_coord[, 1:2]^2)
# gene family from CNV state category
var_coord$family <- stringr::str_split_fixed(var_coord$category, "_", 2)[,1]
# contribution by gene family
fam_contrib <- aggregate(contrib ~ family, data = var_coord, sum)
top_families <- fam_contrib[order(-fam_contrib$contrib), ]
top_families$fam_size<-cnv_out[match(top_families$family,cnv_out$family_id),2]
top_cat <- var_coord[order(-var_coord$contrib), ][1:10, ]
write.csv(top_families,"data/DefenseNstress/AT/eco/AT_MCA_top_influential_families.csv")

arrow_df <- as.data.frame(envfit_res$vectors$arrows)
arrow_df$label <- rownames(arrow_df)

fam_contrib<-ggplot(var_coord, aes(x = Dim.1, y = Dim.2)) +
  geom_point(alpha = 0.2, color = "gray20") +
  geom_point(data = top_cat, aes(x = Dim.1, y = Dim.2), color = 2, size = 2) +
  geom_text_repel(data = top_cat,
                  aes(label = family),
                  size = 3,
                  color = "black",
                  max.overlaps = 20) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = 4,
               inherit.aes = FALSE) +
  geom_text(data = arrow_df,
            aes(x = Dim.1, y = Dim.2, label = label),
            color = 4,
            hjust = 1, vjust = 1,
            inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Top Influential Gene Families with Environmental Fit (MCA)",
       x = "Dimension 1", y = "Dimension 2")


# plot all together in one documents
pdf("plots/AT_MCA_defNstres_families_136_assm_pop_corrct_and_uncorrct_v5.pdf",h=5,w=6)
print(p_uncorrected)
print(p_corrected)
print(fam_contrib)
dev.off()



# LYRATA ----------------------------------

## 1. Covariance matrix
# ***** calculation is in the script pop_structure.R

### Ecology vs CNV per family [with MCA] ------------------
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(data.table)
library(vegan)
# ** Run multiple regression on cnv vs env on the def & stress gene families ***
env_comb<-read.csv("data/DefenseNstress/AL/eco/AL_chelsa_bioclim_for_all_pop_Apr28_25.csv")
svd_df<-read.csv("data/DefenseNstress/AL/genDist/AL_gen_div_covar_svd.csv")
og_vs_ass<-read.table("data/DefenseNstress/AL/cnv_within/geneFamily_vs_assembly_cnv_v5.txt",h=T)
pca_input <- as.data.frame(og_vs_ass[,-c(20,24)]) # remove TE11 and MN47_v4; samples with assembly issues
colnames(pca_input)<-c("family_id","73_3a","AL08","al1","AL27","ALR","BAM12","BOR","F1_14","KN003","KN005","KN006","MN47","NLA","NT1","NT12","NT8","NT9","PU6","TE4","TE8","WS1")
rownames(pca_input)<-pca_input$family_id


## pure MCA first 
pca_matrix <- as.data.frame(t(pca_input[,-1]))
pca_matrix[] <- lapply(pca_matrix, as.factor)
# MCA 
mca_result <- MCA(pca_matrix, graph = F)
# Extract variance explained
explained_variance <- mca_result$eig[, 2]  # Percentage of variance per component
pc1_var <- round(explained_variance[1], 1)
pc2_var <- round(explained_variance[2], 1)

# Extract MCA coordinates
mca_df <- as.data.frame(mca_result$ind$coord)
mca_df$Assembly <- rownames(mca_df)
colnames(mca_df)<-c(paste0("Dim",1:5),"Assembly")

# add groupings
ggroup<-svd_df[match(mca_df$Assembly,svd_df$Sample),]
mca_df$group<-ggroup$Group
mca_df$LR_cluster<-ggroup$LR_cluster

# Visualize MCA with variance in axis labels
pdf("plots/AL_MCAdim34_of_defNstres_v5.pdf",h=5,w=8)
p_uncorrected<-
  ggplot(mca_df, aes(x = Dim1, y = Dim2, label = Assembly)) +
  geom_point(
    size = 3,
    aes(
      color = group,     
      shape = factor(LR_cluster)
    ),
    stroke = 1 # Border thickness
  ) +
  scale_shape_manual(
    values = c(21, 22, 23),
    name = "Genetic Cluster"
  ) +
  # Color scale for clusters
  scale_color_brewer(
    palette = "Set1",
    name = "Geographic Group"
  )+
  # Text labels for assemblies
  geom_text_repel(
    aes(label = Assembly),
    size = 3,
    color = "black",
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  theme_minimal() +
  labs(title = "A.lyrata: MCA of LR-Assemblies Based on CNV in Def.Stress Gene Families",
       x = paste0("Dimension 1 (", pc1_var, "% variance)"),
       y = paste0("Dimension 2 (", pc2_var, "% variance)"))+
  guides(
    color = guide_legend(override.aes = list(shape = 15)),  
    shape = guide_legend(override.aes = list(fill = "white"))
  )

print(p_uncorrected)
dev.off()


## correct for Structure ----
pca_env<-env_comb[match(rownames(mca_df),env_comb$Population_name),]
svd_df_match<-svd_df[match(pca_env$Population_name,svd_df$Sample),]

Variables<-pca_env
ev1<-Variables[-c(1:7)]
ccr<-cor(ev1,use="pairwise.complete.obs")
corrplot::corrplot(ccr,type="lower")
diag(ccr)<-0
thr<-max(abs(ccr))
while(thr>0.99){
  ccr<-cor(ev1,use="pairwise.complete.obs")
  diag(ccr)<-0
  thr<-max(abs(ccr))
  ev1<-ev1[,-(which.max(colSums(abs(ccr))))]
}
lyrs<-colnames(ev1)
evars<-lyrs
eq<-paste0(evars,collapse = "+")
ccr<-cor(ev1,use="pairwise.complete.obs")
corrplot::corrplot(ccr,type="lower")

## Remove structure ----
mca_res<-mca_result
mca_coord <- data.frame(mca_res$ind$coord)
# Fit environmental variables
envfit_res <- envfit(mca_coord, Variables[,evars], permutations = 999)
arrow_df <- as.data.frame(envfit_res$vectors$arrows)
arrow_df$label <- rownames(arrow_df)
arrow_df$p_value <- envfit_res$vectors$pvals
write.csv(arrow_df,"data/DefenseNstress/AL/eco/AL_MCA_significance_before.csv")

arrow_df0 <- as.data.frame(envfit_res$vectors$arrows)
arrow_df0$p_value <- envfit_res$vectors$pvals
pop_structure <- as.data.frame(svd_df_match[,c("PC2")])/3
pop_structure<-data.frame(apply(pop_structure,2,scale))
# Regress out geography from MCA dimensions
mca_corrected <- apply(mca_res$ind$coord, 2, function(y) {
  residuals(lm(y ~ ., data = pop_structure))
})
mca_corrected<-data.frame(mca_corrected)
# fit envfit on geography-corrected coordinates
envfit_res <- envfit(mca_corrected, Variables[,evars], permutations = 999)
mca_corrected$Assembly <- rownames(mca_coord)
ggroup<-svd_df[match(mca_df$Assembly,svd_df$Sample),]
mca_corrected$group<-ggroup$Group
mca_corrected$LR_cluster<-factor(ggroup$LR_cluster)
# Get arrow data as a dataframe
arrow_df <- as.data.frame(envfit_res$vectors$arrows)
arrow_df$label <- rownames(arrow_df)
arrow_df$p_value <- envfit_res$vectors$pvals
arrow_df<-arrow_df[arrow_df0$p_value<0.1,]
write.csv(arrow_df,"data/DefenseNstress/AL/eco/AL_MCA_significance_after.csv")


# If pc1_var and pc2_var are in decimal form (e.g., 0.12 for 12%), convert:
pc1_label <- paste0("Dimension 1 (", pc1_var , "% variance)")
pc2_label <- paste0("Dimension 2 (", pc2_var , "% variance)")

mca_corrected$LR_cluster <- as.factor(mca_corrected$LR_cluster)
mca_corrected$group <- as.factor(mca_corrected$group)


# Plot
p_corrected<-ggplot(mca_corrected, aes(x = Dim.1, y = Dim.2, label = Assembly)) +
  geom_point(
    size = 3,
    aes(
      color = group,
      shape = LR_cluster 
    ),
    stroke = 1 # Border thickness
  ) +
  scale_shape_manual(
    values = c(21, 22, 23),
    name = "Genetic Cluster"
  ) +
  # Color scale for clusters
  scale_color_brewer(
    palette = "Set1",
    name = "Geographic Group"
  ) +
  # Text labels for assemblies
  geom_text_repel(
    aes(label = Assembly),
    size = 3,
    color = "black",
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  # Environment arrows
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black",
    inherit.aes = FALSE
  ) +
  # Arrow labels
  geom_text(
    data = arrow_df,
    aes(x = Dim.1, y = Dim.2, label = label),
    color = "red",
    hjust = 1,
    vjust = 1,
    inherit.aes = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "A. lyrata: MCA EAA Based on CNV in Def./Stress Gene Families",
    x = pc1_label,
    y = pc2_label
  ) +
  
  # Legend adjustments
  guides(
    color = guide_legend(override.aes = list(shape = 15)),  # Squares for color legend
    shape = guide_legend(override.aes = list(fill = "white"))  # White fill for shape legend
  )


# Identify most influential gene families using geography-corrected MCA ----
mca_var <- mca_result$var
var_coord <- as.data.frame(mca_var$coord)
colnames(var_coord)<-paste0("Dim.",1:5)
var_coord$category <- rownames(var_coord)
# approximate influence: distance from origin
var_coord$contrib <- rowSums(var_coord[, 1:2]^2)
# gene family from CNV state category
var_coord$family <- stringr::str_split_fixed(var_coord$category, "_", 2)[,1]
# contribution by gene family
fam_contrib <- aggregate(contrib ~ family, data = var_coord, sum)
top_families <- fam_contrib[order(-fam_contrib$contrib), ]
top_families$fam_size<-cnv_out[match(top_families$family,cnv_out$family_id),2]
top_cat <- var_coord[order(-var_coord$contrib), ][1:10, ]
write.csv(top_families,"data/DefenseNstress/AL/eco/AL_MCA_top_influential_families.csv")

arrow_df <- as.data.frame(envfit_res$vectors$arrows)
arrow_df$label <- rownames(arrow_df)

fam_contrib<-ggplot(var_coord, aes(x = Dim.1, y = Dim.2)) +
  geom_point(alpha = 0.2, color = "gray20") +
  geom_point(data = top_cat, aes(x = Dim.1, y = Dim.2), color = 2, size = 2) +
  geom_text_repel(data = top_cat,
                  aes(label = family),
                  size = 3,
                  color = "black",
                  max.overlaps = 20) +
  geom_segment(data = arrow_df,
               aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = 4,
               inherit.aes = FALSE) +
  geom_text(data = arrow_df,
            aes(x = Dim.1, y = Dim.2, label = label),
            color = 4,
            hjust = 1, vjust = 1,
            inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Top Influential Gene Families with Environmental Fit (MCA)",
       x = "Dimension 1", y = "Dimension 2")


pdf("plots/AL_MCA_defNstres_families_corrct_and_uncorrct_v5.pdf",h=5,w=6)
print(p_uncorrected)
print(p_corrected)
print(fam_contrib)
dev.off()

