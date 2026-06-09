
# THALIANA ------------------

# =========================================================
### 7. ENRICHMENT / OVERLAP ANALYSES ------------------------
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(purrr)
  library(broom)
  library(ggrepel)
  library(patchwork)
  library(tidyverse)
  library(scales)
})
source("scripts/rda_helpers.R")

in_file <- "outputs/AT/eco/main_df_final_CNV_ENV_ready_imputed.csv"
out_dir <- "outputs/AT/eco/cnv_env_univariate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

main_df <- data.table::fread(in_file)
main_df<-main_df[!main_df$accession%in%c("GCA_020911765.2", "GCA_036942435.1", "GCA_946499705.1"),]


structure_cols <- intersect(c("PC1", "PC2"), names(main_df))
if (length(structure_cols) < 2) stop("PC1 and/or PC2 missing")

og_cols <- grep("^OG", names(main_df), value = TRUE)
if (length(og_cols) == 0) stop("No OG columns found")

bio_cols  <- grep("^bio\\d+$", names(main_df), value = TRUE)
soil_cols <- grep("250m", names(main_df), value = TRUE)

# leave out TE
env_cols <- unique(c(bio_cols, soil_cols))
env_cols <- intersect(env_cols, names(main_df))
if (length(env_cols) == 0) stop("No climate/soil variables found")

at_meta<-read.csv("data/231/AT/genDist/AT_gen_div_covar_svd_v2.csv")
ogList<-data.frame(fread("data/ogList_with_functions_v2.csv"))
function_short<-data.frame(fread("data/og_fam_function_short.tsv",h=T))
function_short<-function_short[match(ogList$Family.ID,function_short$Family.ID),]
ogList$func_short<-function_short$short_title

main_df$clean_name <- clean_sample_name(main_df$pop_name)

# exact match between assembly and cluster sample name
geo_group <- at_meta[match(main_df$clean_name, at_meta$Population), , drop = FALSE]
main_df<-data.frame(main_df,geo_group[,c("Group","LR_cluster")])

family_env_z_results<-read.csv(file.path(out_dir, "familywise_zscore_environment_associations.csv"))

# add family functinal information
og_matched<-ogList[match(family_env_z_results$family_id,ogList$Family.ID),c("func_short","FunctionalCategory")]
family_env_z_results<-data.frame(short_func=og_matched[,1],func_cat=og_matched[,2],family_env_z_results)

# -----------------------------
# A. define significant families
# -----------------------------
# use FDR hits if present, otherwise nominal p < 0.05
if (sum(family_env_z_results$q.value < 0.05, na.rm = TRUE) > 0) {
  sig_fams <- family_env_z_results %>%
    filter(q.value < 0.05) %>%
    distinct(family_id) %>%
    pull(family_id)
  
  sig_label <- "FDR < 0.05"
} else {
  sig_fams <- family_env_z_results %>%
    filter(p.value < 0.05) %>%
    distinct(family_id) %>%
    pull(family_id)
  
  sig_label <- "p < 0.05"
}

tested_fams <- family_env_z_results %>%
  distinct(family_id) %>%
  pull(family_id)

cat("Using significant-family definition:", sig_label, "\n")
cat("Number of significant families:", length(sig_fams), "\n")
cat("Number of tested families:", length(tested_fams), "\n")

# =========================================================
### 7A. Functional-category enrichment ------------------------
# =========================================================

fam_annot <- ogList %>%
  transmute(
    family_id = Family.ID,
    FunctionalCategory = FunctionalCategory
  ) %>%
  filter(family_id %in% tested_fams) %>%
  distinct()

fam_annot$FunctionalCategory[is.na(fam_annot$FunctionalCategory)] <- "Unknown"

sig_tbl <- tibble(family_id = tested_fams) %>%
  mutate(is_sig = family_id %in% sig_fams) %>%
  left_join(fam_annot, by = "family_id")

functional_enrichment <- sig_tbl %>%
  group_by(FunctionalCategory) %>%
  summarise(
    sig_in_cat    = sum(is_sig, na.rm = TRUE),
    nonsig_in_cat = sum(!is_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sig_not_cat    = sum(sig_tbl$is_sig) - sig_in_cat,
    nonsig_not_cat = sum(!sig_tbl$is_sig) - nonsig_in_cat
  ) %>%
  rowwise() %>%
  mutate(
    fisher_p = fisher.test(
      matrix(c(sig_in_cat, nonsig_in_cat, sig_not_cat, nonsig_not_cat),
             nrow = 2, byrow = TRUE),
      alternative = "greater"
    )$p.value,
    odds_ratio = {
      m <- matrix(c(sig_in_cat, nonsig_in_cat, sig_not_cat, nonsig_not_cat),
                  nrow = 2, byrow = TRUE)
      ft <- fisher.test(m, alternative = "greater")
      unname(ft$estimate)
    },
    prop_sig_in_cat = sig_in_cat / (sig_in_cat + nonsig_in_cat),
    prop_sig_global = sum(sig_tbl$is_sig) / nrow(sig_tbl)
  ) %>%
  ungroup() %>%
  mutate(
    q.value = p.adjust(fisher_p, method = "BH")
  ) %>%
  arrange(q.value, fisher_p, desc(odds_ratio))

write.csv(
  functional_enrichment,
  file.path(out_dir, "functional_category_enrichment_in_significant_families.csv"),
  row.names = FALSE
)

# Visualize it 
enrich_df<-functional_enrichment

# clean up
plot_df <- enrich_df %>%
  mutate(
    FunctionalCategory = as.character(FunctionalCategory),
    odds_ratio = as.numeric(odds_ratio),
    fisher_p = as.numeric(fisher_p),
    q.value = as.numeric(q.value),
    sig_in_cat = as.numeric(sig_in_cat),
    prop_sig_in_cat = as.numeric(prop_sig_in_cat),
    neglog10_q = -log10(pmax(q.value, 1e-300)),
    neglog10_p = -log10(pmax(fisher_p, 1e-300)),
    is_fdr_sig = q.value < 0.05
  ) %>%
  arrange(desc(odds_ratio), desc(sig_in_cat))

# optional: remove categories with zero significant families
plot_df2 <- plot_df %>%
  filter(sig_in_cat > 0)

# reorder categories for plotting
plot_df2 <- plot_df2 %>%
  mutate(
    FunctionalCategory = forcats::fct_reorder(FunctionalCategory, odds_ratio)
  )

# dot plot
p_enrich_dot <- ggplot(plot_df2, aes(x = odds_ratio, y = FunctionalCategory)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
  geom_point(aes(size = sig_in_cat, color = neglog10_q)) +
  scale_color_gradient(
    low = "grey75",
    high = "#b2182b",
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name = "No. significant\nfamilies") +
  theme_classic(base_size = 11) +
  labs(
    x = "Odds ratio (enrichment)",
    y = NULL,
    title = "Functional category enrichment among significant CNV–environment families"
  )

ggsave(
  file.path(out_dir, "Fig_functional_category_enrichment_dotplot.pdf"),
  p_enrich_dot, width = 6.5, height = 4
)

p_enrich_dot


# =========================================================
### 7B. Sample-level group/cluster association per family ----
# =========================================================
# Here we test whether each family shows CNV differences by Group or LR_cluster.
# Then we ask whether env-significant families are overrepresented among those.

run_family_group_lm <- function(df, group_var, covars = c("PC1", "PC2")) {
  # require at least 2 levels
  if (length(unique(na.omit(df[[group_var]]))) < 2) return(tibble())
  
  df[[group_var]] <- as.factor(df[[group_var]])
  
  rhs <- paste(c(group_var, covars), collapse = " + ")
  form <- as.formula(paste0("cnv_z ~ ", rhs))
  
  fit <- lm(form, data = df)
  a <- anova(fit)
  
  if (!(group_var %in% rownames(a))) return(tibble())
  
  tibble(
    p.value = a[group_var, "Pr(>F)"],
    statistic = a[group_var, "F value"]
  )
}

# rebuild long table including Group and LR_cluster
cnv_long_groups <- main_df %>%
  select(all_of(c("assembly", "accession", "Group", "LR_cluster", structure_cols, og_cols))) %>%
  pivot_longer(
    cols = all_of(og_cols),
    names_to = "family_id",
    values_to = "copy_number"
  ) %>%
  mutate(copy_number = as.numeric(copy_number)) %>%
  filter(family_id %in% tested_fams)

# family-wise z-score again
fam_stats2 <- cnv_long_groups %>%
  group_by(family_id) %>%
  summarise(
    mean_copy = mean(copy_number, na.rm = TRUE),
    sd_copy = sd(copy_number, na.rm = TRUE),
    .groups = "drop"
  )

cnv_long_groups <- cnv_long_groups %>%
  left_join(fam_stats2, by = "family_id") %>%
  mutate(
    cnv_z = ifelse(sd_copy > 0, (copy_number - mean_copy) / sd_copy, NA_real_)
  ) %>%
  select(-mean_copy, -sd_copy) %>%
  mutate(
    across(all_of(structure_cols), ~ as.numeric(scale(.x))),
    Group = as.factor(Group),
    LR_cluster = as.factor(LR_cluster)
  )

# family-wise tests for Group
family_group_results <- cnv_long_groups %>%
  filter(!is.na(cnv_z), !is.na(Group), !is.na(PC1), !is.na(PC2)) %>%
  group_by(family_id) %>%
  group_modify(~ {
    tryCatch(
      run_family_group_lm(.x, group_var = "Group", covars = structure_cols),
      error = function(e) tibble()
    )
  }) %>%
  ungroup() %>%
  mutate(
    group_var = "Group",
    q.value = p.adjust(p.value, method = "BH")
  )

# family-wise tests for LR_cluster
family_cluster_results <- cnv_long_groups %>%
  filter(!is.na(cnv_z), !is.na(LR_cluster), !is.na(PC1), !is.na(PC2)) %>%
  group_by(family_id) %>%
  group_modify(~ {
    tryCatch(
      run_family_group_lm(.x, group_var = "LR_cluster", covars = structure_cols),
      error = function(e) tibble()
    )
  }) %>%
  ungroup() %>%
  mutate(
    group_var = "LR_cluster",
    q.value = p.adjust(p.value, method = "BH")
  )

group_cluster_family_tests <- bind_rows(family_group_results, family_cluster_results)

write.csv(
  group_cluster_family_tests,
  file.path(out_dir, "familywise_group_cluster_associations.csv"),
  row.names = FALSE
)

# =========================================================
### 7C. Are env-significant families enriched for Group/LR_cluster effects?
# =========================================================

test_overlap_enrichment <- function(test_df, label_name) {
  # choose FDR if present, otherwise nominal
  if (sum(test_df$q.value < 0.05, na.rm = TRUE) > 0) {
    assoc_fams <- test_df %>%
      filter(q.value < 0.05) %>%
      distinct(family_id) %>%
      pull(family_id)
    assoc_label <- "FDR < 0.05"
  } else {
    assoc_fams <- test_df %>%
      filter(p.value < 0.05) %>%
      distinct(family_id) %>%
      pull(family_id)
    assoc_label <- "p < 0.05"
  }
  
  a <- sum(tested_fams %in% sig_fams & tested_fams %in% assoc_fams)
  b <- sum(tested_fams %in% sig_fams & !(tested_fams %in% assoc_fams))
  c <- sum(!(tested_fams %in% sig_fams) & tested_fams %in% assoc_fams)
  d <- sum(!(tested_fams %in% sig_fams) & !(tested_fams %in% assoc_fams))
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
  
  tibble(
    variable = label_name,
    assoc_definition = assoc_label,
    n_tested = length(tested_fams),
    n_env_sig = length(intersect(sig_fams, tested_fams)),
    n_assoc = length(intersect(assoc_fams, tested_fams)),
    n_overlap = a,
    odds_ratio = unname(ft$estimate),
    fisher_p = ft$p.value
  )
}

overlap_enrichment <- bind_rows(
  test_overlap_enrichment(family_group_results, "Group"),
  test_overlap_enrichment(family_cluster_results, "LR_cluster")
) %>%
  mutate(q.value = p.adjust(fisher_p, method = "BH"))

write.csv(
  overlap_enrichment,
  file.path(out_dir, "overlap_enrichment_env_significant_vs_group_cluster_associated.csv"),
  row.names = FALSE
)

# =========================================================
### 7D. Optional: which Group / cluster tends to show highest CNV
# =========================================================
# This is descriptive, not a formal enrichment test.

family_group_direction <- cnv_long_groups %>%
  group_by(family_id, Group) %>%
  summarise(mean_cnv_z = mean(cnv_z, na.rm = TRUE), .groups = "drop_last") %>%
  slice_max(order_by = mean_cnv_z, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(is_env_sig = family_id %in% sig_fams)

family_cluster_direction <- cnv_long_groups %>%
  group_by(family_id, LR_cluster) %>%
  summarise(mean_cnv_z = mean(cnv_z, na.rm = TRUE), .groups = "drop_last") %>%
  slice_max(order_by = mean_cnv_z, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(is_env_sig = family_id %in% sig_fams)

write.csv(
  family_group_direction,
  file.path(out_dir, "descriptive_top_group_per_family.csv"),
  row.names = FALSE
)

write.csv(
  family_cluster_direction,
  file.path(out_dir, "descriptive_top_cluster_per_family.csv"),
  row.names = FALSE
)

## visualize
df <- family_group_direction

# -------------------------
# counts
# -------------------------
count_df <- df %>%
  group_by(Group, is_env_sig) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = is_env_sig, values_from = n, values_fill = 0) %>%
  rename(
    n_non_sig = `FALSE`,
    n_sig = `TRUE`
  ) %>%
  mutate(
    total = n_sig + n_non_sig,
    prop_sig = n_sig / total
  )

# -------------------------
# plot
# -------------------------
p_group_prop <- ggplot(count_df, aes(x = Group, y = prop_sig)) +
  geom_col(fill = "#4575b4") +
  geom_hline(
    yintercept = mean(df$is_env_sig),
    linetype = 2,
    color = "grey40"
  ) +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 11) +
  labs(
    x = "Geographic group",
    y = "Proportion of significant families",
    title = "Relative enrichment of CNV–environment associations by geographic group",
    subtitle = "Dashed line = all-family-wide baseline"
  )

ggsave(
  file.path(out_dir, "Fig_group_enrichment_proportion.pdf"),
  p_group_prop, width = 5, height = 3.5
)





# LYRATA ------------------

# =========================================================
### 7. ENRICHMENT / OVERLAP ANALYSES ------------------------
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(purrr)
  library(broom)
  library(ggrepel)
  library(patchwork)
  library(tidyverse)
  library(scales)
})
source("scripts/rda_helpers.R")

in_file <- "outputs/AL/eco/main_df_final_CNV_ENV_ready_imputed.csv"
out_dir <- "outputs/AL/eco/cnv_env_univariate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

main_df <- data.table::fread(in_file)
# tetraploids
tet<-c("WS1","al1","AL27","PU6")
## add ploidy
main_df$ploidy <- ifelse(main_df$Population_name %in% tet, 4, 2)

structure_cols <- intersect(c("PC1", "PC2"), names(main_df))
if (length(structure_cols) < 2) stop("PC1 and/or PC2 missing")

og_cols <- grep("^OG", names(main_df), value = TRUE)
if (length(og_cols) == 0) stop("No OG columns found")

bio_cols  <- grep("^bio\\d+$", names(main_df), value = TRUE)
soil_cols <- grep("250m", names(main_df), value = TRUE)

# leave out TE
env_cols <- unique(c(bio_cols, soil_cols))
env_cols <- intersect(env_cols, names(main_df))
if (length(env_cols) == 0) stop("No climate/soil variables found")



al_meta<-read.csv("data/231/AL/genDist/AL_gen_div_covar_svd.csv")
ogList<-data.frame(fread("data/ogList_with_functions_v2.csv"))
function_short<-data.frame(fread("data/og_fam_function_short.tsv",h=T))
function_short<-function_short[match(ogList$Family.ID,function_short$Family.ID),]
ogList$func_short<-function_short$short_title

# exact match between assembly and cluster sample name
geo_group <- al_meta[match(main_df$Population_name, al_meta$Sample), , drop = FALSE]
main_df<-data.frame(main_df,geo_group[,c("Group","LR_cluster")])

family_env_z_results<-read.csv(file.path(out_dir, "familywise_zscore_environment_associations.csv"))

# add family functinal information
og_matched<-ogList[match(family_env_z_results$family_id,ogList$Family.ID),c("func_short","FunctionalCategory")]
family_env_z_results<-data.frame(short_func=og_matched[,1],func_cat=og_matched[,2],family_env_z_results)

# -----------------------------
# A. define significant families
# -----------------------------
# use FDR hits if present, otherwise nominal p < 0.05
if (sum(family_env_z_results$q.value < 0.05, na.rm = TRUE) > 0) {
  sig_fams <- family_env_z_results %>%
    filter(q.value < 0.05) %>%
    distinct(family_id) %>%
    pull(family_id)
  
  sig_label <- "FDR < 0.05"
} else {
  sig_fams <- family_env_z_results %>%
    filter(p.value < 0.05) %>%
    # filter(abs(estimate) > 0.5) %>%
    distinct(family_id) %>%
    pull(family_id)
  
  sig_label <- "p < 0.05"
}

tested_fams <- family_env_z_results %>%
  distinct(family_id) %>%
  pull(family_id)

cat("Using significant-family definition:", sig_label, "\n")
cat("Number of significant families:", length(sig_fams), "\n")
cat("Number of tested families:", length(tested_fams), "\n")

# =========================================================
### 7A. Functional-category enrichment ------------------------
# =========================================================

fam_annot <- ogList %>%
  transmute(
    family_id = Family.ID,
    FunctionalCategory = FunctionalCategory
  ) %>%
  filter(family_id %in% tested_fams) %>%
  distinct()

fam_annot$FunctionalCategory[is.na(fam_annot$FunctionalCategory)] <- "Unknown"

sig_tbl <- tibble(family_id = tested_fams) %>%
  mutate(is_sig = family_id %in% sig_fams) %>%
  left_join(fam_annot, by = "family_id")

functional_enrichment <- sig_tbl %>%
  group_by(FunctionalCategory) %>%
  summarise(
    sig_in_cat    = sum(is_sig, na.rm = TRUE),
    nonsig_in_cat = sum(!is_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sig_not_cat    = sum(sig_tbl$is_sig) - sig_in_cat,
    nonsig_not_cat = sum(!sig_tbl$is_sig) - nonsig_in_cat
  ) %>%
  rowwise() %>%
  mutate(
    fisher_p = fisher.test(
      matrix(c(sig_in_cat, nonsig_in_cat, sig_not_cat, nonsig_not_cat),
             nrow = 2, byrow = TRUE),
      alternative = "greater"
    )$p.value,
    odds_ratio = {
      m <- matrix(c(sig_in_cat, nonsig_in_cat, sig_not_cat, nonsig_not_cat),
                  nrow = 2, byrow = TRUE)
      ft <- fisher.test(m, alternative = "greater")
      unname(ft$estimate)
    },
    prop_sig_in_cat = sig_in_cat / (sig_in_cat + nonsig_in_cat),
    prop_sig_global = sum(sig_tbl$is_sig) / nrow(sig_tbl)
  ) %>%
  ungroup() %>%
  mutate(
    q.value = p.adjust(fisher_p, method = "BH")
  ) %>%
  arrange(q.value, fisher_p, desc(odds_ratio))

write.csv(
  functional_enrichment,
  file.path(out_dir, "functional_category_enrichment_in_significant_families.csv"),
  row.names = FALSE
)

# Visualize it 
enrich_df<-functional_enrichment

# clean up
plot_df <- enrich_df %>%
  mutate(
    FunctionalCategory = as.character(FunctionalCategory),
    odds_ratio = as.numeric(odds_ratio),
    fisher_p = as.numeric(fisher_p),
    q.value = as.numeric(q.value),
    sig_in_cat = as.numeric(sig_in_cat),
    prop_sig_in_cat = as.numeric(prop_sig_in_cat),
    neglog10_q = -log10(pmax(q.value, 1e-300)),
    neglog10_p = -log10(pmax(fisher_p, 1e-300)),
    is_fdr_sig = q.value < 0.05
  ) %>%
  arrange(desc(odds_ratio), desc(sig_in_cat))

# optional: remove categories with zero significant families
plot_df2 <- plot_df %>%
  filter(sig_in_cat > 0)

# reorder categories for plotting
plot_df2 <- plot_df2 %>%
  mutate(
    FunctionalCategory = forcats::fct_reorder(FunctionalCategory, odds_ratio)
  )

# dot plot
p_enrich_dot <- ggplot(plot_df2, aes(x = odds_ratio, y = FunctionalCategory)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
  geom_point(aes(size = sig_in_cat, color = neglog10_q)) +
  scale_color_gradient(
    low = "grey75",
    high = "#b2182b",
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name = "No. significant\nfamilies") +
  theme_classic(base_size = 11) +
  labs(
    x = "Odds ratio (enrichment)",
    y = NULL,
    title = "Functional category enrichment among significant CNV–environment families"
  )

ggsave(
  file.path(out_dir, "Fig_functional_category_enrichment_dotplot.pdf"),
  p_enrich_dot, width = 6.5, height = 4
)

p_enrich_dot


# =========================================================
### 7B. Sample-level group/cluster association per family ----
# =========================================================
# Here we test whether each family shows CNV differences by Group or LR_cluster.
# Then we ask whether env-significant families are overrepresented among those.

run_family_group_lm <- function(df, group_var, covars = c("PC1", "PC2")) {
  # require at least 2 levels
  if (length(unique(na.omit(df[[group_var]]))) < 2) return(tibble())
  
  df[[group_var]] <- as.factor(df[[group_var]])
  
  rhs <- paste(c(group_var, covars), collapse = " + ")
  form <- as.formula(paste0("cnv_z ~ ", rhs))
  
  fit <- lm(form, data = df)
  a <- anova(fit)
  
  if (!(group_var %in% rownames(a))) return(tibble())
  
  tibble(
    p.value = a[group_var, "Pr(>F)"],
    statistic = a[group_var, "F value"]
  )
}

# rebuild long table including Group and LR_cluster
cnv_long_groups <- main_df %>%
  select(all_of(c("Population_name", "Group", "LR_cluster", structure_cols, og_cols))) %>%
  pivot_longer(
    cols = all_of(og_cols),
    names_to = "family_id",
    values_to = "copy_number"
  ) %>%
  mutate(copy_number = as.numeric(copy_number)) %>%
  filter(family_id %in% tested_fams)

# family-wise z-score again
fam_stats2 <- cnv_long_groups %>%
  group_by(family_id) %>%
  summarise(
    mean_copy = mean(copy_number, na.rm = TRUE),
    sd_copy = sd(copy_number, na.rm = TRUE),
    .groups = "drop"
  )

cnv_long_groups <- cnv_long_groups %>%
  left_join(fam_stats2, by = "family_id") %>%
  mutate(
    cnv_z = ifelse(sd_copy > 0, (copy_number - mean_copy) / sd_copy, NA_real_)
  ) %>%
  select(-mean_copy, -sd_copy) %>%
  mutate(
    across(all_of(structure_cols), ~ as.numeric(scale(.x))),
    Group = as.factor(Group),
    LR_cluster = as.factor(LR_cluster)
  )

# family-wise tests for Group
family_group_results <- cnv_long_groups %>%
  filter(!is.na(cnv_z), !is.na(Group), !is.na(PC1), !is.na(PC2)) %>%
  group_by(family_id) %>%
  group_modify(~ {
    tryCatch(
      run_family_group_lm(.x, group_var = "Group", covars = structure_cols),
      error = function(e) tibble()
    )
  }) %>%
  ungroup() %>%
  mutate(
    group_var = "Group",
    q.value = p.adjust(p.value, method = "BH")
  )

# family-wise tests for LR_cluster
family_cluster_results <- cnv_long_groups %>%
  filter(!is.na(cnv_z), !is.na(LR_cluster), !is.na(PC1), !is.na(PC2)) %>%
  group_by(family_id) %>%
  group_modify(~ {
    tryCatch(
      run_family_group_lm(.x, group_var = "LR_cluster", covars = structure_cols),
      error = function(e) tibble()
    )
  }) %>%
  ungroup() %>%
  mutate(
    group_var = "LR_cluster",
    q.value = p.adjust(p.value, method = "BH")
  )

group_cluster_family_tests <- bind_rows(family_group_results, family_cluster_results)

write.csv(
  group_cluster_family_tests,
  file.path(out_dir, "familywise_group_cluster_associations.csv"),
  row.names = FALSE
)

# =========================================================
### 7C. Are env-significant families enriched for Group/LR_cluster effects?
# =========================================================

test_overlap_enrichment <- function(test_df, label_name) {
  # choose FDR if present, otherwise nominal
  if (sum(test_df$q.value < 0.05, na.rm = TRUE) > 0) {
    assoc_fams <- test_df %>%
      filter(q.value < 0.05) %>%
      distinct(family_id) %>%
      pull(family_id)
    assoc_label <- "FDR < 0.05"
  } else {
    assoc_fams <- test_df %>%
      filter(p.value < 0.05) %>%
      distinct(family_id) %>%
      pull(family_id)
    assoc_label <- "p < 0.05"
  }
  
  a <- sum(tested_fams %in% sig_fams & tested_fams %in% assoc_fams)
  b <- sum(tested_fams %in% sig_fams & !(tested_fams %in% assoc_fams))
  c <- sum(!(tested_fams %in% sig_fams) & tested_fams %in% assoc_fams)
  d <- sum(!(tested_fams %in% sig_fams) & !(tested_fams %in% assoc_fams))
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
  
  tibble(
    variable = label_name,
    assoc_definition = assoc_label,
    n_tested = length(tested_fams),
    n_env_sig = length(intersect(sig_fams, tested_fams)),
    n_assoc = length(intersect(assoc_fams, tested_fams)),
    n_overlap = a,
    odds_ratio = unname(ft$estimate),
    fisher_p = ft$p.value
  )
}

overlap_enrichment <- bind_rows(
  test_overlap_enrichment(family_group_results, "Group"),
  test_overlap_enrichment(family_cluster_results, "LR_cluster")
) %>%
  mutate(q.value = p.adjust(fisher_p, method = "BH"))

write.csv(
  overlap_enrichment,
  file.path(out_dir, "overlap_enrichment_env_significant_vs_group_cluster_associated.csv"),
  row.names = FALSE
)

# =========================================================
### 7D. Optional: which Group / cluster tends to show highest CNV
# =========================================================
# This is descriptive, not a formal enrichment test.

family_group_direction <- cnv_long_groups %>%
  group_by(family_id, Group) %>%
  summarise(mean_cnv_z = mean(cnv_z, na.rm = TRUE), .groups = "drop_last") %>%
  slice_max(order_by = mean_cnv_z, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(is_env_sig = family_id %in% sig_fams)

family_cluster_direction <- cnv_long_groups %>%
  group_by(family_id, LR_cluster) %>%
  summarise(mean_cnv_z = mean(cnv_z, na.rm = TRUE), .groups = "drop_last") %>%
  slice_max(order_by = mean_cnv_z, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(is_env_sig = family_id %in% sig_fams)

write.csv(
  family_group_direction,
  file.path(out_dir, "descriptive_top_group_per_family.csv"),
  row.names = FALSE
)

write.csv(
  family_cluster_direction,
  file.path(out_dir, "descriptive_top_cluster_per_family.csv"),
  row.names = FALSE
)

## visualize
df <- family_group_direction

# -------------------------
# counts
# -------------------------
count_df <- df %>%
  group_by(Group, is_env_sig) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = is_env_sig, values_from = n, values_fill = 0) %>%
  rename(
    n_non_sig = `FALSE`,
    n_sig = `TRUE`
  ) %>%
  mutate(
    total = n_sig + n_non_sig,
    prop_sig = n_sig / total
  )

# -------------------------
# plot
# -------------------------
p_group_prop <- ggplot(count_df, aes(x = Group, y = prop_sig)) +
  geom_col(fill = "#4575b4") +
  geom_hline(
    yintercept = mean(df$is_env_sig),
    linetype = 2,
    color = "grey40"
  ) +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 11) +
  labs(
    x = "Geographic group",
    y = "Proportion of significant families",
    title = "Relative enrichment of CNV–environment associations by geographic group",
    subtitle = "Dashed line = all-family-wide baseline"
  )

ggsave(
  file.path(out_dir, "Fig_group_enrichment_proportion.pdf"),
  p_group_prop, width = 5, height = 3.5
)

