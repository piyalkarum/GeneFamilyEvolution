########################### GENETIC DIVERSITY #########################
# Analysis of population structure.
#=======================================================================

library(data.table)
library(vegan)
library(ggplot2)
library(ggrepel)
library(corrplot)

# 1. THATIANA -----------------
# Variant calls were obtained from the long read sequences of 136 A.thaliana assemblies.
## subsample VCFs to extract populations
at_neut<-data.table::fread("data/DefenseNstress/AT/genDist/all_subsampled_neutral_snps.txt")
hder<-data.table::fread("data/DefenseNstress/AT/genDist/header.txt")
colnames(at_neut)<-colnames(hder)
gt<-rCNV::hetTgen(at_neut,"GT-012")
# select only the samples available in 136 
er_136_names<-data.table::fread("data/LR_ER_names_matched_edited.csv")
gt1<-data.frame(gt[,c(1:4,match(er_136_names$ys_er_names,colnames(gt)))])
rm(at_neut,gt)
gc()
mult_alle<-which(nchar(gt1$ATT)>1)
gt2<-gt1[-mult_alle,]
# calculate allele frequencies
# Assuming diploid: allele freq = (# alt alleles) / (2 * number of non-missing genotypes)
allele_freqs <- rCNV:::apply_pb(gt2[,-c(1:4)], 1, function(g) {
  g<-as.numeric(g)
  non_na <- !is.na(g)
  if (sum(non_na) == 0) return(rep(NA, (ncol(gt)-4)))
  alt_alleles <- g[non_na]
  freq <- alt_alleles / 2
  freq_padded <- rep(NA, (ncol(gt2)-4))
  freq_padded[non_na] <- freq
  freq_padded
})
write.table(allele_freqs,"data/DefenseNstress/AT/genDist/at_neut_AF_full_table.txt",row.names = FATSE,col.names = FATSE, sep = "\t",quote = FATSE)
# Split SNPs (columns) extract 10 10k subsets randomly
# individuals × SNPs
set.seed(42)
allele_freqs<-data.table::fread("data/DefenseNstress/AT/genDist/at_neut_AF_full_table.txt",h=F)
allele_freqs<-data.frame(allele_freqs)
allele_freqs<-allele_freqs[,which(apply(allele_freqs,2,function(x)sum(is.na(x)))<(136*.05))]
n_snps <- ncol(allele_freqs)
snp_sample_size <- 10000  
# Generate 10 random samples of 10k SNPs each
for (i in 1:6) {
  sampled_snps <- sample(1:n_snps, size = snp_sample_size, replace = FATSE)
  subset <- allele_freqs[, sampled_snps, drop = FATSE]
  write.table(subset, 
              paste0("data/DefenseNstress/AT/genDist/at_neut_AF_random10k_", i, ".txt"),
              row.names = FATSE, 
              col.names = FATSE, 
              sep = "\t",
              quote = FATSE)
}

# generate covariance mats
af_files<-list.files("data/DefenseNstress/AT/genDist",full.names = T, pattern="at_neut_AF_random10k_")
pb<-txtProgressBar(max=5,style=3,width=50)
for( i in seq_along(af_files)){
  tm<-scale(read.table(af_files[i],h=F), center = TRUE, scale = FATSE)
  tm_omega<-cov(t(tm),use="pairwise.complete.obs")
  write.table(tm_omega,paste0("data/DefenseNstress/AT/genDist/at_neut_CoVar_mat_subset_",i,".txt"),row.names = F, col.names = F, sep="\t",quote=F)
  setTxtProgressBar(pb,i)
}

# Average the covariance matrices
cov_fls<-list.files("data/DefenseNstress/AT/genDist/",pattern = "at_neut_CoVar_mat_subset",full.names=T)
cov_list<-lapply(cov_fls,read.table,header=F)
cov_matrix_avg <- Reduce("+", cov_list) / length(cov_list)
write.table(cov_matrix_avg,"data/DefenseNstress/AT/genDist/at_neut_CoVar_mean_mat_1-6_submats.txt",row.names = F, col.names = F, sep="\t",quote=F)
# cov_matrix_avg is the mean covariance matrix
# scale it like BayPass would (mean-centered, variance standardized)
cov_matrix_scaled <- scale(cov_matrix_avg, center = TRUE, scale = TRUE)
write.table(cov_matrix_scaled,"data/DefenseNstress/AT/genDist/at_neut_CoVar_mean_mat_1-6_submats_scaled.txt",row.names = F, col.names = F, sep="\t",quote=F)



#### plot covariance matrix -------------------
v<-read.table("data/DefenseNstress/AT/genDist/at_neut_CoVar_mean_mat_1-6_submats_scaled.txt",h=F)
er_136_names<-data.table::fread("data/LR_ER_names_matched_edited.csv",h=T)
pop.names<-er_136_names$lr_nam
dimnames(v)=list(pop.names,pop.names)
lr_met<-data.table::fread("data/ncbi_LR_136dataset_cleaned_w_sampleInfo_with_qgraph_genetic_group2.1.csv",h=T)
lr_met$sample_name.1<-stringr::str_split_fixed(lr_met$sample_name.1,"\\ ",n=2)[,1]

#drop problematic ones
drop_samples<-c("GCA_036926965.1","GCA_036937515.1","GCA_036940625.1","GCA_036926925.1","GCA_020911765.2","GCA_946499705.1","GCA_036926975.1","GCA_946407795.1")
lr_drop<-lr_met[match(drop_samples,lr_met$Assembly.Accession),]$sample_name.1
index<-!colnames(v)%in%lr_drop
v1<-v[index,index]
er_136_lr_meta_match<-lr_met[match(colnames(v1),lr_met$sample_name.1),]
er_136_lr_meta_match[106,25]<-"CapeVerde"

# get single value decomposition
svd_result <- svd(v1)
# variance explained
variance_explained <- (svd_result$d^2) / sum(svd_result$d^2) * 100
# Extract principal components (SVD provides U, S, V)
pc1 <- svd_result$u[, 1] * svd_result$d[1] 
pc2 <- svd_result$u[, 2] * svd_result$d[2]
svd_df <- data.frame(Sample = colnames(v1), PC1 = pc1, PC2 = pc2,Group = er_136_lr_meta_match$geo_group)
write.csv(svd_df,"data/DefenseNstress/AT/genDist/AT_gen_div_covar_svd.csv",row.names=F)

# Plot SVD (like PCA)
library(ggplot2)
library(ggrepel)
library(corrplot)
svd_df<-read.csv("data/DefenseNstress/AT/genDist/AT_gen_div_covar_svd.csv")

group_colors <- setNames(hcl.colors(length(unique(svd_df$Group))), unique(svd_df$Group))
p<-ggplot(svd_df, aes(x = PC1, y = PC2, label = Sample, color=Group)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 10, box.padding = 0.5, point.padding = 0.3) +
  theme_minimal() +
  ggtitle(bquote("SVD of" ~ Omega ~ "for A. thaliana long-reads")) +
  xlab(paste0("PC1 (", round(variance_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  scale_color_manual(values = group_colors) +
  theme(legend.position = "right")

# extract clusters
cor.mat <- cov2cor(as.matrix(v1))
dist_mat <- as.dist(1 - cor.mat)  # Distance matrix
hc <- hclust(dist_mat, method = "average")  # Hierarchical clustering
clusters <- cutree(hc, k = 4)
cluster_members <- split(names(clusters), clusters)
at_lr_clusters<-data.frame(Cluster = clusters, Sample = names(clusters))
write.csv(at_lr_clusters,file = "data/DefenseNstress/AT/genDist/AT_covar_cluster_assignments.csv", row.names = FATSE)

## adding new clusters
cls<-at_lr_clusters[match(svd_df$Sample,at_lr_clusters$Sample),]$Cluster
svd_df$LR_cluster<-factor(cls)
write.csv(svd_df,"data/DefenseNstress/AT/genDist/AT_gen_div_covar_svd.csv",row.names=F)

group_colors2 <- setNames(hcl.colors(length(unique(svd_df$LR_cluster))), unique(svd_df$LR_cluster))
p2<-ggplot(svd_df, aes(x = PC1, y = PC2, label = Sample, color=LR_cluster)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 10, box.padding = 0.5, point.padding = 0.3) +
  theme_minimal() +
  ggtitle(bquote("SVD of" ~ Omega ~ "for A. thaliana long-reads")) +
  xlab(paste0("PC1 (", round(variance_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  scale_color_manual(values = group_colors2) +
  theme(legend.position = "right")



combined_plot <- ggplot(svd_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(aes(color = Group, shape = LR_cluster), size = 3,stroke=1.5) +  
  geom_text_repel(size = 3, max.overlaps = 10, box.padding = 0.5, point.padding = 0.3) +
  theme_minimal() +
  ggtitle(bquote("SVD of" ~ Omega ~ "for A. thaliana long-reads")) +
  xlab(paste0("PC1 (", round(variance_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  # scale_color_manual("Group", values = group_colors) + 
  scale_color_brewer(palette = "Set1", name = "Geographic Group") +  # Use same palette
  scale_shape_manual("LR Cluster", values = 1:length(unique(svd_df$LR_cluster))) +  
  theme(legend.position = "right",
        legend.box = "vertical")

pdf("plots/AT_gen_div_plots_YS_data_v5.pdf",w=6,h=5)
print(p)
print(p2)
print(combined_plot)
cor.mat=cov2cor(as.matrix(v1))
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,main=expression("Correlation map based on"~hat(Omega)))
hclust.ave <- function(x) hclust(x, method="average")
heatmap(1-cor.mat,hclustfun = hclust.ave,main=expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()


###### 3D plot ---------------
pc3 <- svd_result$u[, 3] * svd_result$d[3]  # Extract PC3
svd_df$PC3 <- pc3
group_colors <- setNames(hcl.colors(length(unique(svd_df$Group))), unique(svd_df$Group))
shapes <- c("1" = "circle", "2" = "square", "3" = "diamond", "4" = "cross")
# Create 3D plot
p_3d <- plot_ly(
  data = svd_df,
  x = ~PC1, 
  y = ~PC2, 
  z = ~PC3,
  color = ~Group,
  colors = group_colors,
  symbol = ~LR_cluster,
  symbols = shapes,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 6),
  text = ~Sample,  # Hover text
  hoverinfo = "text"
) %>%
  layout(
    title = "3D SVD of Ω for A. thaliana long-reads",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(variance_explained[1], 2), "%)")),
      yaxis = list(title = paste0("PC2 (", round(variance_explained[2], 2), "%)")),
      zaxis = list(title = paste0("PC3 (", round(variance_explained[3], 2), "%)"))
    ),
    legend = list(orientation = "h")
  )

p_3d

htmlwidgets::saveWidget(p_3d, "plots/3D_SVD_plot.html")








### LYRATA -----------------------
library(rCNV)
## GenDist from Anna's vcf 
anna_meta<-data.table::fread("data/AL_all_Metadata_PN_lab.txt",h=T)
header<-data.table::fread("data/DefenseNstress/AT/genDist/All_lyrata_final_biallelic_header.txt",h=T)
AL_sample_names<-data.frame(colnames(header)[-c(1:9)])

al_meta<-read.table("data/Alyr_all_LR_meta_data_wgroups.txt",h=T)

nn<-al_meta$Population_name
nn[6:8]<-c("TE_11","TE_4","TE_8")
nn[5]<-"PU_6"
nn[9]<-"BAM"
nn[10]<-"VLH"
nn[11]<-"MAL"
nn[12]<-"lyrMI"
nn[13]<-"lyrON1"
nn[15]<-"lyrON2"
nn[18]<-"lyrPA4"
nn[14]<-"OSL"
nn[17]<-"WS_1"
nn[16]<-"al"
nn[19:22]<-"PLE"

## missing populations 
# c("BOR","KN005","KN006","KN003")
vcf_names_match<-lapply(nn,function(x){
  grep(x,AL_sample_names[,1])
})
names(vcf_names_match)<-al_meta$Population_name
nam_join<-NULL
for(i in seq_along(vcf_names_match)){
  tm<-cbind(names(vcf_names_match[i]),vcf_names_match[[i]])
  nam_join<-rbind(nam_join,tm)
}
nam_join<-data.frame(nam_join)
colnames(nam_join)<-c("population","sample")
nam_join$vcf_nam<-AL_sample_names[nam_join$sample,1]
write.csv(nam_join,"data/DefenseNstress/AT/genDist/AL_1000_PN_vcf_subset2.csv",row.names=F)

## extract samples from the vcf
bcftools view -S to_extract_from_vcf.txt -m2 -M2 -v snps All_lyrata_final_biallelic.vcf.gz | bcftools view -e 'AC=0 || AC=AN' -Oz -o All_lyrata_final_biallelic_subset_for_genDiv.vcf.gz

bcftools view --samples-file to_extract_from_vcf.txt --min-alleles 2 All_lyrata_final_biallelic.vcf.gz -Oz -o All_lyrata_final_biallelic_subset_for_genDiv_0.vcf.gz

bcftools view --max-alleles 2 All_lyrata_final_biallelic_subset_for_genDiv_0.vcf.gz -Oz -o All_lyrata_final_biallelic_subset_for_genDiv_1.vcf.gz

bcftools view --max-alleles 2 --min-alleles 2 --types snps All_lyrata_final_biallelic_subset_for_genDiv_1.vcf.gz | bcftools filter -e "F_MISSING > 0.1" -Oz -o All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_missing.vcf.gz

bcftools query \
-f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_missing.vcf.gz \
> All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_missing.tsv

## extract putatively neutral sites
gt<-data.table::fread("data/DefenseNstress/AT/genDist/All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_missing.tsv",h=F)
hh<-data.table::fread("All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_missing_header.txt",h=T)
hh<-colnames(hh)[-c(5:9)]
colnames(gt)<-hh

data_names<-read.csv("AL_1000_PN_vcf_subset_clean.csv")
data_names<-data_names[match(colnames(gt)[-c(1:4)],data_names$vcf_nam),]
data_names$vcf_nam<-gsub("-",".",data_names$vcf_nam)

genotype_to_count <- function(gt_string) {
  if (gt_string == "./." | gt_string == "./././.") return(NA)
  alleles <- as.numeric(unlist(strsplit(gt_string, "/")))
  sum(alleles, na.rm = TRUE) / length(alleles)
}

nm<-unique(data_names$population)
subs<-sample(1:1264, size = nrow(gt),replace = T)

# processing loop
pb <- txtProgressBar(max = 10, style = 3, width = 50)
for (j in 1:10) {
  setTxtProgressBar(pb, j)
  tm2 <- gt[subs == j, ]
  tm2 <- data.frame(tm2)
  pops <- NULL
  
  for (i in seq_along(nm)) {
    current_pop <- nm[i]
    samples <- data_names$vcf_nam[data_names$population == current_pop]
    tm <- tm2[, samples, drop = FALSE]
    
    allele_counts <- matrix(
      sapply(unlist(tm), genotype_to_count),
      nrow = nrow(tm),
      ncol = ncol(tm))
    
    af <- rowMeans(allele_counts, na.rm = TRUE)
    pops <- cbind(pops, af)
  }
  
  pops[,1]<-rowMeans(pops[,2:4])
  colnames(pops) <- nm
  write.table(pops,paste0("data/DefenseNstress/AT/genDist/Allele_freq_All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_subset",j,".txt"),row.names = F, sep="\t",quote=F)
  
  # Calculate covariance matrix
  cov_mat <- cov(pops, use = "pairwise.complete.obs")
  write.table(cov_mat,paste0("/mnt/hpc-project/CNV/ANN/genDist/AL/covmat_sub/CovMat_All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_subset",j,".txt"),row.names = F, sep="\t",quote=F)
}

# Average the covariance matrices
cov_fls<-list.files("data/DefenseNstress/AT/genDist/covmat_sub",pattern = "CovMat",full.names=T)
cov_list<-lapply(cov_fls,read.table,header=T)
cov_matrix_avg <- Reduce("+", cov_list) / length(cov_list)
colnames(cov_matrix_avg)[19]<-"73_3a"
rownames(cov_matrix_avg)<-colnames(cov_matrix_avg)
# cov_matrix_avg[,1]<-rowMeans(cov_matrix_avg[,2:4])
write.table(cov_matrix_avg,"data/DefenseNstress/AT/genDist/CovMat_meanMat_All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_subset.txt",row.names = F, sep="\t",quote=F)
# scale it like BayPass would (mean-centered, variance standardized)
cov_matrix_scaled <- scale(cov_matrix_avg, center = TRUE, scale = TRUE)
write.table(cov_matrix_scaled,"data/DefenseNstress/AT/genDist/CovMat_meanMat_scaled_All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_subset.txt",row.names = F, sep="\t",quote=F)


v<-read.table("data/DefenseNstress/AT/genDist/CovMat_meanMat_scaled_All_lyrata_final_biallelic_subset_for_genDiv_1_filtered_10percent_subset.txt",h=T)
colnames(v)[19]<-"73_3a"
rownames(v)<-colnames(v)
al_meta<-read.table("/Users/piyalkaru/Desktop/DDORF/DATA/Alyr/Alyr_all_LR_meta_data_wgroups.txt",h=T)
al_meta_match<-al_meta[match(rownames(v),al_meta$Population_name),]
svd_result <- svd(v)
# variance explained
variance_explained <- (svd_result$d^2) / sum(svd_result$d^2) * 100
# Extract principal components (SVD provides U, S, V)
pc1 <- svd_result$u[, 1] * svd_result$d[1] 
pc2 <- svd_result$u[, 2] * svd_result$d[2]
svd_df <- data.frame(Sample = colnames(v), PC1 = pc1, PC2 = pc2,Group = al_meta_match$gen_cluster)
write.csv(svd_df,"data/DefenseNstress/AT/genDist/AL_gen_div_covar_svd.csv",row.names=F)

# Plot SVD (like PCA)
library(ggplot2)
library(ggrepel)
library(corrplot)
group_colors <- setNames(hcl.colors(length(unique(svd_df$Group))), unique(svd_df$Group))
p<-ggplot(svd_df, aes(x = PC1, y = PC2, label = Sample, color=Group)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 10, box.padding = 0.5, point.padding = 0.3) +
  theme_minimal() +
  ggtitle(bquote("SVD of" ~ Omega ~ "for A. lyrata long-reads")) +
  xlab(paste0("PC1 (", round(variance_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  scale_color_manual(values = group_colors) +
  theme(legend.position = "right")

# extract clusters
cor.mat <- cov2cor(as.matrix(v))
dist_mat <- as.dist(1 - cor.mat)  
hc <- hclust(dist_mat, method = "average")  
clusters <- cutree(hc, k = 3)
cluster_members <- split(names(clusters), clusters)
at_lr_clusters<-data.frame(Cluster = clusters, Sample = names(clusters))
write.csv(at_lr_clusters,file = "data/DefenseNstress/AT/genDist/AL_cluster_assignments.csv", row.names = FALSE)

## adding new clusters
cls<-at_lr_clusters[match(svd_df$Sample,at_lr_clusters$Sample),]$Cluster
svd_df$LR_cluster<-factor(cls)
write.csv(svd_df,"data/DefenseNstress/AT/genDist/AL_gen_div_covar_svd.csv",row.names=F)

group_colors2 <- setNames(hcl.colors(length(unique(svd_df$LR_cluster))), unique(svd_df$LR_cluster))
p2<-ggplot(svd_df, aes(x = PC1, y = PC2, label = Sample, color=LR_cluster)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 10, box.padding = 0.5, point.padding = 0.3) +
  theme_minimal() +
  ggtitle(bquote("SVD of" ~ Omega ~ "for A. lyrata long-reads")) +
  xlab(paste0("PC1 (", round(variance_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  scale_color_manual(values = group_colors2) +
  theme(legend.position = "right")

combined_plot <- ggplot(svd_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(aes(color = Group, shape = LR_cluster), size = 3,stroke=1.5) +  
  geom_text_repel(size = 3, max.overlaps = 30, box.padding = 0.5, point.padding = 0.3) +
  theme_minimal() +
  ggtitle(bquote("SVD of" ~ Omega ~ "for A. lyrata long-reads")) +
  xlab(paste0("PC1 (", round(variance_explained[1], 2), "%)")) +
  ylab(paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  # scale_color_manual("Group", values = group_colors) + 
  scale_color_brewer(palette = "Set1", name = "Geographic Group") +
  scale_shape_manual("LR Cluster", values = 1:length(unique(svd_df$LR_cluster))) +  
  theme(legend.position = "right",
        legend.box = "vertical")

pdf("plots/AL_gen_div_plots_WGS_data_v5.pdf", w=6, h=5)
print(p) 
print(p2) 
print(combined_plot) 
cor.mat <- cov2cor(as.matrix(v))
corrplot(cor.mat, method="color", mar=c(2,1,2,2)+0.1, main=expression("Correlation map based on"~hat(Omega)))
hclust.ave <- function(x) hclust(x, method="average")
heatmap(1-cor.mat, hclustfun = hclust.ave, main=expression("Heatmap of "~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()



