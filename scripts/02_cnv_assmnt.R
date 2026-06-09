#########################################################
####                 CNV ASSESSMENT                 #####
#########################################################

setwd("/Users/piyalkaru/Desktop/DDORF/Ann/231ogs/cnv_geneFam")

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(vegan)
  library(patchwork)
  library(ggrepel)
  library(pheatmap)
  library(data.table)
})

dir.create("plots", showWarnings = FALSE)

#########################################################
####                    INPUT                       #####
#########################################################

species1_name <- "AT"
species2_name <- "AL"

sp1 <- read.table("data/231/AT/cnv_within/Athaliana231_geneFamily_vs_assembly_cnv_v1.1.txt",h = TRUE)
sp2 <- read.table("data/231/AL/cnv_within/Alyrata231_geneFamily_vs_assembly_cnv_v1.1.txt",h = TRUE)

# Remove problematic / duplicated assemblies
sp1 <- sp1[, -c(2, 69, 122)]
sp2 <- sp2[, -15]

at_meta<-read.csv("data/231/AT/genDist/AT_gen_div_covar_svd_v2.csv")
al_meta<-read.csv("data/231/AL/genDist/AL_gen_div_covar_svd.csv")

# fix sample and population names
at_names<-at_meta[match(colnames(sp1)[-1],at_meta$Sample),2]
colnames(sp1)[2:ncol(sp1)]<-at_names

at_meta<-data.frame(Sample=at_meta$Population,at_meta[,-c(1,2)])
al_meta$Sample[19]<-"X73_3a"


#########################################################
####                HELPER FUNCTIONS                #####
#########################################################
source("scripts/cnv_helpers.R")

#########################################################
####          1. COMBINED DATA PREPARATION          #####
#########################################################

dat <- make_full_long(sp1, sp2, species1_name, species2_name)
df <- dat$df
df_full <- dat$df_full

wm <- make_wide_matrix(df)
wide <- wm$wide
mat <- wm$mat

#########################################################
####      2. COMBINED WITHIN / BETWEEN SPECIES      #####
#########################################################

family_var <- df %>%
  group_by(species, family_id) %>%
  summarise(
    mean_copy = mean(copy),
    var_copy = var(copy),
    cv = sd(copy) / mean(copy),
    .groups = "drop"
  )

m1 <- glmmTMB(
  copy ~ 1 + (1 | family_id) + (1 | assembly),
  data = df,
  family = nbinom2
)

m2 <- glmmTMB(
  copy ~ species + (1 | family_id) + (1 | assembly),
  data = df,
  family = nbinom2
)

summary(m1)
summary(m2)
anova(m1, m2)

family_tests <- df %>%
  group_by(family_id) %>%
  summarise(
    p_value = tryCatch(wilcox.test(copy ~ species)$p.value, error = function(e) NA_real_),
    mean_sp1 = mean(copy[species == species1_name]),
    mean_sp2 = mean(copy[species == species2_name]),
    logFC = log2((mean_sp2 + 1) / (mean_sp1 + 1)),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

write.csv(family_tests,"data/231/between_sp_cnv_family_test.csv")

#########################################################
####         3. COMBINED ASSEMBLY-BIAS CHECK        #####
#########################################################

total_copy_df <- data.frame(
  assembly = wide$assembly,
  species = wide$species,
  total_copy = rowSums(mat)
)


#########################################################
####            4. COMBINED PCA ANALYSES            #####
#########################################################

# Raw combined PCA
pca_raw_combined <- prcomp(log1p(mat), scale. = TRUE)
p_raw_combined <- make_pca_plot(pca_raw_combined, wide,clusters = NULL, max_overlap = 10,"Combined PCA - raw")

# Within-species z-normalized combined PCA
mat_sp1 <- mat[wide$species == species1_name, , drop = FALSE]
mat_sp2 <- mat[wide$species == species2_name, , drop = FALSE]

mat_sp1_z <- t(scale(t(log1p(mat_sp1))))
mat_sp2_z <- t(scale(t(log1p(mat_sp2))))
mat_z_combined <- rbind(mat_sp1_z, mat_sp2_z)

pca_z_combined <- prcomp(mat_z_combined, scale. = FALSE)
p_z_combined <- make_pca_plot(pca_z_combined, wide,clusters = NULL,  max_overlap = 10,"Combined PCA - within-species z-normalized")

# Family-centered combined PCA
mat_fc_combined <- scale(log1p(mat), center = TRUE, scale = FALSE)
pca_fc_combined <- prcomp(mat_fc_combined, scale. = FALSE)
p_fc_combined <- make_pca_plot(pca_fc_combined, wide,clusters = NULL,  max_overlap = 10,"Combined PCA - family-centered")

ggsave("plots/PCA_combined_raw.pdf", p_raw_combined, width = 9, height = 7)
ggsave("plots/PCA_combined_zscore.pdf", p_z_combined, width = 9, height = 7)
ggsave("plots/PCA_combined_family_centered.pdf", p_fc_combined, width = 9, height = 7)
ggsave("plots/PCA_combined_all_three.pdf", p_raw_combined + p_z_combined + p_fc_combined, width = 20, height = 7)

#########################################################
####           5. SEPARATE SPECIES PCA              #####
#########################################################

res_AT <- run_species_pca(sp1, species1_name, clusters = at_meta)
res_AL <- run_species_pca(sp2, species2_name,clusters = al_meta)

ggsave("plots/PCA_AT_raw.pdf", res_AT$p_raw, width = 6, height = 4)
ggsave("plots/PCA_AT_zscore.pdf", res_AT$p_z, width = 8, height = 6)
ggsave("plots/PCA_AT_family_centered.pdf", res_AT$p_fc, width = 8, height = 6)

ggsave("plots/PCA_AL_raw.pdf", res_AL$p_raw, width = 8, height = 6)
ggsave("plots/PCA_AL_zscore.pdf", res_AL$p_z, width = 5.5, height = 4)
ggsave("plots/PCA_AL_family_centered.pdf", res_AL$p_fc, width = 8, height = 6)

ggsave("plots/PCA_species_raw_side_by_side.pdf", res_AT$p_raw + res_AL$p_raw, width = 16, height = 6)
ggsave("plots/PCA_species_zscore_side_by_side.pdf", res_AT$p_z + res_AL$p_z, width = 16, height = 6)
ggsave("plots/PCA_species_family_centered_side_by_side.pdf", res_AT$p_fc + res_AL$p_fc, width = 16, height = 6)

#########################################################
####         6. SEPARATE SPECIES BIAS CHECK         #####
#########################################################

total_copy_AT <- data.frame(
  assembly = res_AT$wide$assembly,
  species = res_AT$wide$species,
  total_copy = rowSums(res_AT$mat)
)

total_copy_AL <- data.frame(
  assembly = res_AL$wide$assembly,
  species = res_AL$wide$species,
  total_copy = rowSums(res_AL$mat)
)

p_total_AT <- ggplot(total_copy_AT, aes(x = assembly, y = total_copy)) +
  geom_col(fill = "#1f78b4") +
  theme_bw() +
  labs(title = paste0(species1_name, ": total copy count per assembly")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_total_AL <- ggplot(total_copy_AL, aes(x = assembly, y = total_copy)) +
  geom_col(fill = "#e31a1c") +
  theme_bw() +
  labs(title = paste0(species2_name, ": total copy count per assembly")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("plots/total_copy_per_assembly_by_species.pdf", p_total_AT + p_total_AL, width = 14, height = 5)

#########################################################
####          7. COMBINED PERMANOVA                 #####
#########################################################

dist_mat <- vegdist(mat_z_combined, method = "euclidean")
permanova_res <- adonis2(dist_mat ~ species, data = wide)
betadisper_res <- betadisper(dist_mat, wide$species)

print(permanova_res)
print(anova(betadisper_res))

#########################################################
####             8. SPECIES HEATMAPS                #####
#########################################################

mat_heat_sp1 <- sp1 %>%
  column_to_rownames("family_id") %>%
  as.matrix()

mat_heat_sp2 <- sp2 %>%
  column_to_rownames("family_id") %>%
  as.matrix()

mat_heat_sp1_z <- t(scale(t(log1p(mat_heat_sp1))))
mat_heat_sp2_z <- t(scale(t(log1p(mat_heat_sp2))))

pdf("plots/heatmap_AT_row_zscore.pdf", width = 8, height = 10)
pheatmap(
  mat_heat_sp1_z,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = paste0(species1_name, ": row-scaled CNV heatmap")
)
dev.off()

pdf("plots/heatmap_AL_row_zscore.pdf", width = 8, height = 10)
pheatmap(
  mat_heat_sp2_z,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = paste0(species2_name, ": row-scaled CNV heatmap")
)
dev.off()

#########################################################
####      9. COMBINED FAMILY SUMMARY PLOTTING       #####
#########################################################

family_stats <- df_full %>%
  group_by(family_id, species) %>%
  summarise(
    n_assemblies = n(),
    max_copy = max(copy),
    min_copy = min(copy),
    median_copy = median(copy),
    mean_copy = mean(copy),
    sd_copy = sd(copy),
    var_copy = var(copy),
    cv = ifelse(mean_copy > 0, sd_copy / mean_copy, NA_real_),
    dispersion = ifelse(mean_copy > 0, var_copy / mean_copy, NA_real_),
    .groups = "drop"
  )

variable_families <- family_stats %>%
  group_by(family_id) %>%
  summarise(is_variable = any(var_copy > 0, na.rm = TRUE), .groups = "drop") %>%
  filter(is_variable) %>%
  pull(family_id)

family_wide <- family_stats %>%
  filter(family_id %in% variable_families) %>%
  select(family_id, species, mean_copy, cv, dispersion) %>%
  pivot_wider(
    names_from = species,
    values_from = c(mean_copy, cv, dispersion),
    names_sep = "__"
  )

family_order <- family_wide %>%
  mutate(order_metric = rowSums(across(starts_with("mean_copy__")), na.rm = TRUE)) %>%
  arrange(order_metric) %>%
  pull(family_id)

family_wide <- family_wide %>%
  mutate(family_id = factor(family_id, levels = family_order))

p_mean <- make_mirror_plot(
  family_wide,
  col_sp1 = paste0("mean_copy__", species1_name),
  col_sp2 = paste0("mean_copy__", species2_name),
  xlab = "Mean copy number",
  title_text = "Mean copy number per variable gene family",
  species1_label = species1_name,
  species2_label = species2_name
)

p_cv <- make_mirror_plot(
  family_wide,
  col_sp1 = paste0("cv__", species1_name),
  col_sp2 = paste0("cv__", species2_name),
  xlab = "Coefficient of variation (CV)",
  title_text = "Copy-number variability per variable gene family",
  species1_label = species1_name,
  species2_label = species2_name
)

p_disp <- make_mirror_plot(
  family_wide,
  col_sp1 = paste0("dispersion__", species1_name),
  col_sp2 = paste0("dispersion__", species2_name),
  xlab = "Dispersion (variance / mean)",
  title_text = "Dispersion per variable gene family",
  species1_label = species1_name,
  species2_label = species2_name
)

ggsave("plots/mirror_plot_mean.pdf", p_mean, width = 9, height = 12)
ggsave("plots/mirror_plot_cv.pdf", p_cv, width = 9, height = 12)
ggsave("plots/mirror_plot_dispersion.pdf", p_disp, width = 9, height = 12)
ggsave("plots/mirror_plot_all_panels.pdf", p_mean / p_cv / p_disp, width = 10, height = 18)

#########################################################
####                  SAVE TABLES                    #####
#########################################################

write.table(family_var, "outputs/family_variation_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(family_tests, "outputs/family_species_tests.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(total_copy_df, "outputs/assembly_total_copy_counts_combined.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(family_stats, "outputs/family_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(as.data.frame(permanova_res), "outputs/permanova_results.tsv", sep = "\t", row.names = TRUE, quote = FALSE)




####### summary for the manuscript #########
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# =========================================================
# Input
# =========================================================
var_file  <- "outputs/family_variation_summary.tsv"
test_file <- "outputs/family_species_tests.tsv"

var_df  <- fread(var_file)
test_df <- fread(test_file)

# =========================================================
# 1. Within-species summary
# =========================================================
# A family is counted as showing within-species CNV if variance > 0
within_species_summary <- var_df %>%
  group_by(species) %>%
  summarise(
    n_families = n(),
    n_variable_families = sum(var_copy > 0, na.rm = TRUE),
    prop_variable_families = n_variable_families / n_families,
    median_variance = median(var_copy, na.rm = TRUE),
    median_cv = median(cv, na.rm = TRUE),
    mean_cv = mean(cv, na.rm = TRUE),
    .groups = "drop"
  )

# =========================================================
# 2. Between-species summary
# =========================================================
between_species_summary <- test_df %>%
  summarise(
    n_families_tested = n(),
    n_sig_p05 = sum(p_value < 0.05, na.rm = TRUE),
    n_sig_fdr05 = sum(FDR < 0.05, na.rm = TRUE),
    prop_sig_p05 = n_sig_p05 / n_families_tested,
    prop_sig_fdr05 = n_sig_fdr05 / n_families_tested,
    median_abs_logFC = median(abs(logFC), na.rm = TRUE),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE)
  )

# =========================================================
# 3. Compact publication-style table
# =========================================================
# reshape within-species table to one row
within_wide <- within_species_summary %>%
  select(species, n_variable_families, prop_variable_families, median_variance, median_cv) %>%
  pivot_wider(
    names_from = species,
    values_from = c(n_variable_families, prop_variable_families, median_variance, median_cv)
  )

small_summary_table <- bind_cols(within_wide, between_species_summary)

# =========================================================
# 4. Optional: per-family classification table
# =========================================================
# This table lets you say which families vary within species and/or between species
var_wide <- var_df %>%
  select(species, family_id, var_copy, cv) %>%
  pivot_wider(
    names_from = species,
    values_from = c(var_copy, cv)
  )

family_level_summary <- test_df %>%
  left_join(var_wide, by = "family_id") %>%
  mutate(
    between_sig_fdr05 = FDR < 0.05,
    between_sig_p05   = p_value < 0.05,
    within_var_sp1 = ifelse(!is.na(var_copy_AT), var_copy_AT > 0,
                            ifelse(!is.na(var_copy_AL), var_copy_AL > 0, NA)),
    within_var_sp2 = case_when(
      "var_copy_AL" %in% names(.) ~ var_copy_AL > 0,
      TRUE ~ NA
    )
  )

# =========================================================
# 5. Save outputs
# =========================================================
write.table(within_species_summary,
            file = "outputs/within_species_cnv_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(between_species_summary,
            file = "outputs/between_species_cnv_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(small_summary_table,
            file = "outputs/small_cnv_summary_table.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(family_level_summary,
            file = "outputs/family_level_cnv_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)



range(var_df$cv[var_df$species == "AL"], na.rm = TRUE)
range(var_df$cv[var_df$species == "AT"], na.rm = TRUE)

cnv_long %>%
  group_by(species, family_id) %>%
  summarise(mean_copy = mean(copy_number)) %>%
  group_by(species) %>%
  summarise(range = paste0(min(mean_copy), "--", max(mean_copy)))

     