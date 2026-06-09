################## ECOLOGY ASSESSMENT ###############

# THALIANA ------------------------------------
setwd("/Users/piyalkaru/Desktop/DDORF/Ann/231ogs/cnv_geneFam")
# All data joined earlier
main_df <- data.table::fread("data/231/AT/eco/AT_main_data_for_eco_v1.1.tsv", h = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(vegan)
  library(caret)
  library(purrr)
  library(broom)
  library(car)
})

out_dir <- "outputs/eco"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1. Define columns ------------------------
## =========================================================

id_cols <- c("accession", "assembly", "sample_name", "location",
             "latitude", "longitude", "pop_name")

# keep structure as available
structure_cols <- c("PC1", "PC2")

tech_cols <- c("busco_complete_pct",
               "contig_n50_mbp",
               "total_length_bp")

# TE as genomic/ecological covariate
te_cols <- c("te_int_frac")

og_cols   <- grep("^OG", names(main_df), value = TRUE)
bio_cols  <- grep("^bio\\d+$", names(main_df), value = TRUE)
soil_cols <- grep("250m", names(main_df), value = TRUE)

# full predictor pool for CNV-ENV analysis
env_cols <- c(bio_cols, soil_cols, te_cols)
env_cols <- unique(env_cols)

## optional: convert assembly size to Mbp
main_df <- main_df %>%
  mutate(total_length_mbp = total_length_bp / 1e6)

tech_cols2 <- c("busco_complete_pct",
                "contig_n50_mbp",
                "total_length_mbp")

## =========================================================
## 2. Technical covariate screening ------------------------
## =========================================================

tech_df <- main_df %>%
  select(all_of(c("assembly", tech_cols2, structure_cols))) %>%
  distinct()

cor_mat <- cor(tech_df[, ..tech_cols2], use = "pairwise.complete.obs", method = "spearman")

write.csv(cor_mat,
          file.path(out_dir, "technical_covariate_correlations.csv"),
          row.names = TRUE)

cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Var1", "Var2", "rho")

p_cor <- ggplot(cor_long, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", midpoint = 0) +
  theme_classic(base_size = 12) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Spearman rho",
       title = "Correlation among technical covariates")

ggsave(file.path(out_dir, "Fig_S1_technical_covariate_correlation_heatmap.pdf"),
       p_cor, width = 6.2, height = 5.2)

highly_cor <- caret::findCorrelation(cor_mat, cutoff = 0.8, names = TRUE)
tech_selected <- setdiff(tech_cols2, highly_cor)
if (length(tech_selected) == 0) tech_selected <- tech_cols2
message("Technical covariates retained: ",
        paste(tech_selected, collapse = ", "))

## =========================================================
## 3. Reshape CNV and keep variable families ------------------------
## =========================================================

cnv_long <- main_df %>%
  select(all_of(c("assembly", structure_cols, tech_cols2, og_cols))) %>%
  pivot_longer(cols = all_of(og_cols),
               names_to = "family_id",
               values_to = "copy_number") %>%
  mutate(copy_number = as.numeric(copy_number))

family_filter <- cnv_long %>%
  group_by(family_id) %>%
  summarise(
    n = n(),
    n_unique = n_distinct(copy_number),
    sd_copy = sd(copy_number, na.rm = TRUE),
    mean_copy = mean(copy_number, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_unique > 1, sd_copy > 0)

cnv_long_filt <- cnv_long %>%
  filter(family_id %in% family_filter$family_id)

## =========================================================
## 4. Standardize copy number within family ------------------------
## =========================================================

cnv_long_filt <- cnv_long_filt %>%
  group_by(family_id) %>%
  mutate(
    cnv_dev = copy_number - median(copy_number, na.rm = TRUE),
    cnv_z = {
      s <- sd(copy_number, na.rm = TRUE)
      if (is.na(s) || s == 0) rep(0, dplyr::n())
      else (copy_number - mean(copy_number, na.rm = TRUE)) / s
    }
  ) %>%
  ungroup()

## =========================================================
## 5. Multivariate screening: ------------------------
##    Does technical variation explain CNV matrix
##    after controlling for structure?
## =========================================================

cnv_mat <- cnv_long_filt %>%
  select(assembly, family_id, cnv_z) %>%
  pivot_wider(names_from = family_id, values_from = cnv_z) %>%
  as.data.frame()

rownames(cnv_mat) <- cnv_mat$assembly
cnv_mat$assembly <- NULL

cnv_mat <- cnv_mat[, colSums(is.na(cnv_mat)) == 0, drop = FALSE]
cnv_mat <- cnv_mat[, apply(cnv_mat, 2, sd) > 0, drop = FALSE]

meta_df <- main_df %>%
  select(all_of(c("assembly", structure_cols, tech_selected))) %>%
  distinct() %>%
  as.data.frame()

rownames(meta_df) <- meta_df$assembly
meta_df$assembly <- NULL

common_ids <- intersect(rownames(cnv_mat), rownames(meta_df))
cnv_mat <- cnv_mat[common_ids, , drop = FALSE]
meta_df <- meta_df[common_ids, , drop = FALSE]

tech_formula <- paste(tech_selected, collapse = " + ")
rda_formula <- as.formula(paste0("cnv_mat ~ ", tech_formula, " + Condition(PC1 + PC2)"))

rda_tech <- vegan::rda(rda_formula, data = meta_df)

rda_overall <- anova.cca(rda_tech, permutations = 999)
rda_terms   <- anova.cca(rda_tech, by = "term", permutations = 999)
rda_rsq     <- RsquareAdj(rda_tech)

write.csv(as.data.frame(rda_overall),
          file.path(out_dir, "technical_rda_overall.csv"),
          row.names = TRUE)
write.csv(as.data.frame(rda_terms),
          file.path(out_dir, "technical_rda_term_tests.csv"),
          row.names = TRUE)

rda_term_df <- as.data.frame(rda_terms)
rda_term_df$term <- rownames(rda_term_df)
rownames(rda_term_df) <- NULL
pcol_rda <- grep("Pr", names(rda_term_df), value = TRUE)[1]

vcol_rda <- if ("Variance" %in% names(rda_term_df)) {
  "Variance"
} else if ("ChiSquare" %in% names(rda_term_df)) {
  "ChiSquare"
} else {
  names(rda_term_df)[2]
}

p_rda <- ggplot(rda_term_df %>% filter(term != "Residual"),
                aes(x = reorder(term, .data[[vcol_rda]]),
                    y = .data[[vcol_rda]],
                    fill = .data[[pcol_rda]] < 0.05)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "gray70")) +
  theme_classic(base_size = 12) +
  labs(x = NULL,
       y = "Explained constrained variance",
       fill = "p < 0.05",
       title = "Multivariate association between CNV and technical covariates",
       subtitle = paste0("Partial RDA controlling for PC1 and PC2; adjusted R² = ",
                         round(rda_rsq$adj.r.squared, 3)))

ggsave(file.path(out_dir, "Fig_S2_multivariate_CNV_vs_technical_covariates.pdf"),
       p_rda, width = 7, height = 4.8)

## =========================================================
## 6. Family-wise technical screening ------------------------
## =========================================================

run_family_lm <- function(df, tech_var) {
  fit <- lm(as.formula(paste0("cnv_z ~ ", tech_var, " + PC1 + PC2")), data = df)
  td <- broom::tidy(fit)
  gl <- broom::glance(fit)
  term_row <- td %>% filter(term == tech_var)
  
  tibble(
    estimate = term_row$estimate,
    std.error = term_row$std.error,
    statistic = term_row$statistic,
    p.value = term_row$p.value,
    r.squared = gl$r.squared,
    adj.r.squared = gl$adj.r.squared
  )
}

family_results <- map_dfr(tech_selected, function(tv) {
  cnv_long_filt %>%
    group_by(family_id) %>%
    group_modify(~ {
      out <- tryCatch(run_family_lm(.x, tv), error = function(e) NULL)
      if (is.null(out)) {
        return(tibble(
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          r.squared = NA_real_,
          adj.r.squared = NA_real_
        ))
      }
      out
    }) %>%
    ungroup() %>%
    mutate(technical_var = tv)
}) %>%
  group_by(technical_var) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write.csv(family_results,
          file.path(out_dir, "familywise_technical_associations.csv"),
          row.names = FALSE)

family_summary <- family_results %>%
  group_by(technical_var) %>%
  summarise(
    n_tested = sum(!is.na(p.value)),
    n_sig_p05 = sum(p.value < 0.05, na.rm = TRUE),
    n_sig_fdr05 = sum(q.value < 0.05, na.rm = TRUE),
    median_abs_effect = median(abs(estimate), na.rm = TRUE),
    median_adj_r2 = median(adj.r.squared, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_fdr05))

write.csv(family_summary,
          file.path(out_dir, "familywise_technical_summary.csv"),
          row.names = FALSE)

p_bar <- family_summary %>%
  mutate(technical_var = fct_reorder(technical_var, n_sig_fdr05)) %>%
  ggplot(aes(x = technical_var, y = n_sig_fdr05)) +
  geom_col(fill = "#377eb8") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL,
       y = "Number of families (FDR < 0.05)",
       title = "Family-wise associations with assembly-quality covariates")

ggsave(file.path(out_dir, "Fig_S3_number_of_families_associated_with_technical_covariates.pdf"),
       p_bar, width = 6.5, height = 4.5)

heat_df <- family_results %>%
  filter(q.value < 0.05) %>%
  mutate(
    technical_var = factor(technical_var, levels = tech_selected),
    family_id = fct_reorder(family_id, abs(estimate), .fun = max, .desc = TRUE)
  )

if (nrow(heat_df) > 0) {
  p_heat <- ggplot(heat_df,
                   aes(x = technical_var, y = family_id, fill = estimate)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
    theme_classic(base_size = 11) +
    labs(x = NULL, y = "Gene family",
         fill = "Effect size",
         title = "Effect sizes of significant technical covariates on family CNV")
  
  # ggsave(file.path(out_dir, "Fig_S4_familywise_effect_heatmap_technical_covariates.pdf"),
  #        p_heat, width = 6.8, height = 8.5)
}

flagged_families <- family_results %>%
  group_by(family_id) %>%
  summarise(
    any_sig_fdr05 = any(q.value < 0.05, na.rm = TRUE),
    n_sig_covars = sum(q.value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(any_sig_fdr05)

write.csv(flagged_families,
          file.path(out_dir, "flagged_families_technical_confounds.csv"),
          row.names = FALSE)

## =========================================================
## 6b. TE content screening: ------------------------
##     Does CNV associate with TE content after controlling
##     for population structure?
## =========================================================

if (!"te_int_frac" %in% names(main_df)) {
  warning("te_int_frac not found in main_df; skipping TE screening")
} else {
  
  ## -------------------------------------------------------
  ## A. Multivariate TE ~ CNV matrix (partial RDA)
  ## -------------------------------------------------------
  meta_te <- main_df %>%
    select(assembly, PC1, PC2, te_int_frac) %>%
    distinct() %>%
    as.data.frame()
  
  rownames(meta_te) <- meta_te$assembly
  meta_te$assembly <- NULL
  
  common_ids_te <- intersect(rownames(cnv_mat), rownames(meta_te))
  cnv_mat_te <- cnv_mat[common_ids_te, , drop = FALSE]
  meta_te <- meta_te[common_ids_te, , drop = FALSE]
  
  # remove rows with missing TE or PCs
  keep_te <- complete.cases(meta_te[, c("PC1", "PC2", "te_int_frac")])
  cnv_mat_te <- cnv_mat_te[keep_te, , drop = FALSE]
  meta_te <- meta_te[keep_te, , drop = FALSE]
  
  if (nrow(meta_te) > 5) {
    rda_te <- vegan::rda(cnv_mat_te ~ te_int_frac + Condition(PC1 + PC2), data = meta_te)
    
    rda_te_overall <- anova.cca(rda_te, permutations = 999)
    rda_te_term <- anova.cca(rda_te, by = "term", permutations = 999)
    rda_te_rsq <- RsquareAdj(rda_te)
    
    write.csv(as.data.frame(rda_te_overall),
              file.path(out_dir, "te_rda_overall.csv"),
              row.names = TRUE)
    
    write.csv(as.data.frame(rda_te_term),
              file.path(out_dir, "te_rda_term_test.csv"),
              row.names = TRUE)
    
    te_term_df <- as.data.frame(rda_te_term)
    te_term_df$term <- rownames(te_term_df)
    rownames(te_term_df) <- NULL
    
    pcol_te <- grep("Pr", names(te_term_df), value = TRUE)[1]
    vcol_te <- if ("Variance" %in% names(te_term_df)) {
      "Variance"
    } else if ("ChiSquare" %in% names(te_term_df)) {
      "ChiSquare"
    } else {
      names(te_term_df)[2]
    }
    
    p_te_rda <- ggplot(te_term_df %>% filter(term != "Residual"),
                       aes(x = reorder(term, .data[[vcol_te]]),
                           y = .data[[vcol_te]],
                           fill = .data[[pcol_te]] < 0.05)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "gray70")) +
      theme_classic(base_size = 12) +
      labs(x = NULL,
           y = "Explained constrained variance",
           fill = "p < 0.05",
           title = "Multivariate association between CNV and TE content",
           subtitle = paste0("Partial RDA controlling for PC1 and PC2; adjusted R² = ",
                             round(rda_te_rsq$adj.r.squared, 3)))
    
    ggsave(file.path(out_dir, "Fig_S4_multivariate_CNV_vs_TE_content.pdf"),
           p_te_rda, width = 6.8, height = 4.5)
  }
  
  ## -------------------------------------------------------
  ## B. Family-wise screening: ------------------------
  ##    cnv_z ~ te_int_frac + PC1 + PC2
  ## -------------------------------------------------------
  run_family_lm_te <- function(df) {
    fit <- lm(cnv_z ~ te_int_frac + PC1 + PC2, data = df)
    td <- broom::tidy(fit)
    gl <- broom::glance(fit)
    term_row <- td %>% filter(term == "te_int_frac")
    
    tibble(
      estimate = term_row$estimate,
      std.error = term_row$std.error,
      statistic = term_row$statistic,
      p.value = term_row$p.value,
      r.squared = gl$r.squared,
      adj.r.squared = gl$adj.r.squared
    )
  }
  
  cnv_long_te <- cnv_long_filt %>%
    left_join(
      main_df %>%
        select(assembly, te_int_frac) %>%
        distinct(),
      by = "assembly"
    ) %>%
    filter(!is.na(te_int_frac), !is.na(PC1), !is.na(PC2))
  
  family_results_te <- cnv_long_te %>%
    group_by(family_id) %>%
    group_modify(~ {
      out <- tryCatch(run_family_lm_te(.x), error = function(e) NULL)
      if (is.null(out)) {
        return(tibble(
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          r.squared = NA_real_,
          adj.r.squared = NA_real_
        ))
      }
      out
    }) %>%
    ungroup() %>%
    mutate(q.value = p.adjust(p.value, method = "BH"))
  
  write.csv(family_results_te,
            file.path(out_dir, "familywise_TE_associations.csv"),
            row.names = FALSE)
  
  family_summary_te <- family_results_te %>%
    summarise(
      n_tested = sum(!is.na(p.value)),
      n_sig_p05 = sum(p.value < 0.05, na.rm = TRUE),
      n_sig_fdr05 = sum(q.value < 0.05, na.rm = TRUE),
      median_abs_effect = median(abs(estimate), na.rm = TRUE),
      median_adj_r2 = median(adj.r.squared, na.rm = TRUE)
    )
  
  write.csv(family_summary_te,
            file.path(out_dir, "familywise_TE_summary.csv"),
            row.names = FALSE)
  
  p_te_bar <- family_summary_te %>%
    pivot_longer(cols = c(n_sig_p05, n_sig_fdr05),
                 names_to = "metric",
                 values_to = "n_families") %>%
    ggplot(aes(x = metric, y = n_families)) +
    geom_col(fill = "#984ea3") +
    theme_classic(base_size = 12) +
    labs(x = NULL,
         y = "Number of families",
         title = "Family-wise associations between CNV and TE content") +
    scale_x_discrete(labels = c("n_sig_p05" = "p < 0.05",
                                "n_sig_fdr05" = "FDR < 0.05"))
  
  ggsave(file.path(out_dir, "Fig_S5_number_of_families_associated_with_TE_content.pdf"),
         p_te_bar, width = 4.8, height = 4.2)
  
  heat_df_te <- family_results_te %>%
    filter(q.value < 0.05) %>%
    mutate(family_id = fct_reorder(family_id, abs(estimate), .desc = TRUE))
  
  if (nrow(heat_df_te) > 0) {
    p_heat_te <- ggplot(heat_df_te,
                        aes(x = "te_int_frac", y = family_id, fill = estimate)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
      theme_classic(base_size = 11) +
      labs(x = NULL, y = "Gene family",
           fill = "Effect size",
           title = "Effect sizes of TE content on family CNV")
    
    ggsave(file.path(out_dir, "Fig_S6_familywise_effect_heatmap_TE_content.pdf"),
           p_heat_te, width = 5.5, height = 8)
  }
  
  flagged_families_te <- family_results_te %>%
    filter(q.value < 0.05)
  
  write.csv(flagged_families_te,
            file.path(out_dir, "flagged_families_TE_confounds.csv"),
            row.names = FALSE)
}

## =========================================================
## 7. Filter gene families before CNV-ENV analysis ------------------------
## =========================================================

family_stats <- main_df %>%
  select(all_of(c("assembly", og_cols))) %>%
  pivot_longer(cols = all_of(og_cols),
               names_to = "family_id",
               values_to = "copy_number") %>%
  group_by(family_id) %>%
  summarise(
    n_individuals = n(),
    mean_copy = mean(copy_number, na.rm = TRUE),
    sd_copy = sd(copy_number, na.rm = TRUE),
    n_states = n_distinct(copy_number),
    n_nonzero = sum(copy_number > 0, na.rm = TRUE),
    zero_fraction = mean(copy_number == 0, na.rm = TRUE),
    min_copy = min(copy_number, na.rm = TRUE),
    max_copy = max(copy_number, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(family_stats,
          file.path(out_dir, "family_filtering_statistics.csv"),
          row.names = FALSE)

min_nonzero <- 5
min_states <- 2
max_zero_fraction <- 0.95

family_stats <- family_stats %>%
  mutate(
    keep_sd = sd_copy > 0,
    keep_nonzero = n_nonzero >= min_nonzero,
    keep_states = n_states >= min_states,
    keep_zero = zero_fraction <= max_zero_fraction,
    keep_final = keep_sd & keep_nonzero & keep_states & keep_zero
  )

families_kept <- family_stats %>% filter(keep_final) %>% pull(family_id)
families_removed <- setdiff(og_cols, families_kept)

write.csv(family_stats,
          file.path(out_dir, "family_filtering_decisions.csv"),
          row.names = FALSE)

write.table(families_kept,
            file.path(out_dir, "families_retained.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(families_removed,
            file.path(out_dir, "families_removed.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

filter_summary <- tibble(
  criterion = c("Invariant (sd = 0)",
                paste0("< ", min_nonzero, " nonzero individuals"),
                paste0("< ", min_states, " copy-number states"),
                paste0("Zero fraction > ", max_zero_fraction),
                "Retained"),
  n_families = c(
    sum(!family_stats$keep_sd),
    sum(!family_stats$keep_nonzero),
    sum(!family_stats$keep_states),
    sum(!family_stats$keep_zero),
    sum(family_stats$keep_final)
  )
)

p_filter <- ggplot(filter_summary,
                   aes(x = fct_reorder(criterion, n_families), y = n_families)) +
  geom_col(fill = "#377eb8") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Number of gene families",
       title = "Filtering of gene families prior to CNV–environment analysis")

ggsave(file.path(out_dir, "Fig_S5_family_filtering_summary.pdf"),
       p_filter, width = 7, height = 4.5)

## =========================================================
## 8. Environmental / soil / TE variable assessment ------------------------
## =========================================================

env_raw <- main_df %>%
  select(all_of(c("assembly", env_cols))) %>%
  distinct()

# 8.1 Missingness assessment
env_missing <- env_raw %>%
  summarise(across(all_of(env_cols), ~ mean(is.na(.x)))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "missing_fraction") %>%
  arrange(desc(missing_fraction))

write.csv(env_missing,
          file.path(out_dir, "environment_missingness_summary.csv"),
          row.names = FALSE)

p_missing <- env_missing %>%
  mutate(variable = fct_reorder(variable, missing_fraction)) %>%
  ggplot(aes(x = variable, y = missing_fraction)) +
  geom_col(fill = "#8c6bb1") +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Fraction missing",
       title = "Missingness across climate, soil, and TE variables")

ggsave(file.path(out_dir, "Fig_S6_environment_missingness.pdf"),
       p_missing, width = 7, height = 6)

# threshold for environmental missingness
max_env_missing <- 0.25

env_keep_missing <- env_missing %>%
  filter(missing_fraction <= max_env_missing) %>%
  pull(variable)

env_drop_missing <- setdiff(env_cols, env_keep_missing)

# 8.2 Build reduced env table
env_df <- env_raw %>%
  select(all_of(c("assembly", env_keep_missing)))

# 8.3 Median imputation for remaining NAs
env_imputed <- env_df %>%
  mutate(across(-assembly, ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))

# 8.4 Remove near-zero variance variables
nzv_idx <- nearZeroVar(env_imputed %>% select(-assembly))
env_nzv_removed <- names(env_imputed %>% select(-assembly))[nzv_idx]

if (length(nzv_idx) > 0) {
  env_df2 <- env_imputed %>% select(-assembly) %>% .[, -nzv_idx, drop = FALSE]
} else {
  env_df2 <- env_imputed %>% select(-assembly)
}

# 8.5 Correlation matrix before collinearity filtering
cor_before <- cor(env_df2, use = "pairwise.complete.obs", method = "spearman")
write.csv(cor_before,
          file.path(out_dir, "environment_correlations_before_filtering.csv"),
          row.names = TRUE)

cor_before_long <- as.data.frame(as.table(cor_before))
colnames(cor_before_long) <- c("Var1", "Var2", "rho")

p_cor_before <- ggplot(cor_before_long, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", midpoint = 0) +
  theme_classic(base_size = 10) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Spearman rho",
       title = "Environmental correlation matrix before filtering") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(out_dir, "Fig_S7_environment_correlation_before.pdf"),
       p_cor_before, width = 8, height = 7)

# 8.6 Remove highly correlated variables
cor_cutoff <- 0.8
high_cor_idx <- findCorrelation(cor_before, cutoff = cor_cutoff, names = FALSE)

if (length(high_cor_idx) > 0) {
  cor_removed <- names(env_df2)[high_cor_idx]
  keep_cols <- setdiff(names(env_df2), cor_removed)
  env_df3 <- env_df2[, ..keep_cols]
} else {
  cor_removed <- character(0)
  env_df3 <- env_df2
}

# 8.7 Iterative VIF filtering
vif_cutoff <- 6

calc_vif_manual <- function(df) {
  vars <- names(df)
  sapply(vars, function(v) {
    others <- setdiff(vars, v)
    if (length(others) == 0) return(1)
    f <- as.formula(paste(v, "~", paste(others, collapse = " + ")))
    r2 <- summary(lm(f, data = df))$r.squared
    1 / (1 - r2)
  })
}

env_df4 <- data.frame(env_df3)
vif_removed <- character(0)

repeat {
  vif_vals <- calc_vif_manual(env_df4)
  max_vif <- max(vif_vals, na.rm = TRUE)
  if (max_vif <= vif_cutoff || ncol(env_df4) <= 1) break
  drop_var <- names(which.max(vif_vals))
  vif_removed <- c(vif_removed, drop_var)
  env_df4 <- env_df4[, setdiff(names(env_df4), drop_var), drop = FALSE]
}

final_env_vars <- names(env_df4)

env_filter_summary <- tibble(
  step = c("Initial variables",
           paste0("Removed missingness > ", max_env_missing),
           "Removed near-zero variance",
           paste0("Removed high correlation (|rho| > ", cor_cutoff, ")"),
           paste0("Removed by VIF > ", vif_cutoff),
           "Final retained variables"),
  n_variables = c(length(env_cols),
                  length(env_drop_missing),
                  length(env_nzv_removed),
                  length(cor_removed),
                  length(vif_removed),
                  length(final_env_vars)),
  variables = c(
    paste(env_cols, collapse = ", "),
    paste(env_drop_missing, collapse = ", "),
    paste(env_nzv_removed, collapse = ", "),
    paste(cor_removed, collapse = ", "),
    paste(vif_removed, collapse = ", "),
    paste(final_env_vars, collapse = ", ")
  )
)

write.csv(env_filter_summary,
          file.path(out_dir, "environment_filtering_summary.csv"),
          row.names = FALSE)

write.table(final_env_vars,
            file.path(out_dir, "environment_variables_retained.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cor_after <- cor(env_df4, use = "pairwise.complete.obs", method = "spearman")
write.csv(cor_after,
          file.path(out_dir, "environment_correlations_after_filtering.csv"),
          row.names = TRUE)

cor_after_long <- as.data.frame(as.table(cor_after))
colnames(cor_after_long) <- c("Var1", "Var2", "rho")

p_cor_after <- ggplot(cor_after_long, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", midpoint = 0) +
  theme_classic(base_size = 10) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Spearman rho",
       title = "Environmental correlation matrix after filtering") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(out_dir, "Fig_S8_environment_correlation_after.pdf"),
       p_cor_after, width = 7, height = 6)

## =========================================================
## 9. Final datasets for CNV-ENV analysis ------------------------
## =========================================================

# raw table (retains NAs)
main_df_final_raw <- main_df %>%
  select(all_of(c(id_cols, structure_cols, tech_cols2, final_env_vars, families_kept)))

write.csv(main_df_final_raw,
          file.path(out_dir, "main_df_final_CNV_ENV_ready_raw.csv"),
          row.names = FALSE)

# imputed table for multivariate / regression workflows
env_imputed_final <- env_imputed %>%
  select(all_of(c("assembly", final_env_vars)))

main_df_final_imputed <- main_df %>%
  select(all_of(c(id_cols, structure_cols, tech_cols2, families_kept))) %>%
  left_join(env_imputed_final, by = "assembly")

write.csv(main_df_final_imputed,
          file.path(out_dir, "main_df_final_CNV_ENV_ready_imputed.csv"),
          row.names = FALSE)

## also save a compact variable-class table
variable_classes <- tibble(
  variable = c(tech_selected, final_env_vars),
  class = c(rep("technical", length(tech_selected)),
            ifelse(final_env_vars %in% bio_cols, "climate",
                   ifelse(final_env_vars %in% soil_cols, "soil", "TE")))
)

write.csv(variable_classes,
          file.path(out_dir, "variable_classes_for_CNV_ENV.csv"),
          row.names = FALSE)


## CNV–ENV ANALYSIS ------------------------------------
#=======================================================

################## CNV–ENV PARTIAL RDA ##################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(vegan)
  library(glmmTMB)
  library(broom)
  library(purrr)
  library(stringr)
  library(scales)
  library(ggrepel)
  library(patchwork)
})

# ====================================================
### 0. Input ------------------------
# ====================================================

in_file <- "outputs/AT/eco/main_df_final_CNV_ENV_ready_imputed.csv"
varclass_file <- "outputs/AT/eco/variable_classes_for_CNV_ENV.csv"
out_dir <- "outputs/AT/eco/cnv_env_rda"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

main_df <- fread(in_file)
var_classes <- fread(varclass_file)

# drop problematic assemblies << assmelby issues
main_df<-main_df[!main_df$accession%in%c("GCA_020911765.2", "GCA_036942435.1", "GCA_946499705.1"),]

# ====================================================
### 1. Define variable sets ------------------------
# ====================================================

id_cols <- intersect(c("accession", "assembly", "sample_name", "location",
                       "latitude", "longitude", "pop_name", "Population"), names(main_df))

structure_cols <- intersect(c("PC1", "PC2"), names(main_df))
if (length(structure_cols) < 2) stop("PC1 and/or PC2 missing")

og_cols <- grep("^OG", names(main_df), value = TRUE)
if (length(og_cols) == 0) stop("No OG columns found")

# remove TE and technical variables from analysis
climate_vars <- var_classes %>%
  filter(class == "climate") %>%
  pull(variable) %>%
  intersect(names(main_df))

soil_vars <- var_classes %>%
  filter(class == "soil") %>%
  pull(variable) %>%
  intersect(names(main_df))

env_vars <- unique(c(climate_vars, soil_vars))

if (length(env_vars) == 0) stop("No climate/soil variables found")
if (length(climate_vars) == 0) warning("No climate variables found")
if (length(soil_vars) == 0) warning("No soil variables found")

message("Climate vars: ", paste(climate_vars, collapse = ", "))
message("Soil vars: ", paste(soil_vars, collapse = ", "))

# ====================================================
### 2. Build response matrices ------------------------
# ====================================================

cnv_wide <- main_df %>%
  select(accession, all_of(og_cols)) %>%
  distinct() %>%
  as.data.frame()

rownames(cnv_wide) <- cnv_wide$accession
cnv_wide$accession <- NULL

# keep only variable families
keep_fams <- apply(cnv_wide, 2, function(x) sd(x, na.rm = TRUE) > 0)
cnv_raw <- cnv_wide[, keep_fams, drop = FALSE]

# z-score within family
cnv_z <- scale(cnv_raw, center = TRUE, scale = TRUE)
cnv_z <- as.data.frame(cnv_z)

# remove any families that became all-NA after scaling
cnv_z <- cnv_z[, colSums(is.na(cnv_z)) == 0, drop = FALSE]

# align raw to z columns
shared_fams <- intersect(colnames(cnv_raw), colnames(cnv_z))
cnv_raw <- cnv_raw[, shared_fams, drop = FALSE]
cnv_z   <- cnv_z[, shared_fams, drop = FALSE]

# ====================================================
### 3. Build predictor matrices ------------------------
# ====================================================

meta_df <- main_df %>%
  select(all_of(c("accession", id_cols, structure_cols, climate_vars, soil_vars))) %>%
  distinct() %>%
  as.data.frame()

rownames(meta_df) <- meta_df$accession
meta_df$accession <- NULL

# scale predictors for RDA
scale_if_present <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (length(vars) > 0) {
    df[, vars] <- lapply(df[, vars, drop = FALSE], function(x) as.numeric(scale(x)))
  }
  df
}

meta_df <- scale_if_present(meta_df, c(structure_cols, climate_vars, soil_vars))

# align rows
common_ids <- Reduce(intersect, list(
  rownames(cnv_raw),
  rownames(cnv_z),
  rownames(meta_df)
))

cnv_raw <- cnv_raw[common_ids, , drop = FALSE]
cnv_z   <- cnv_z[common_ids, , drop = FALSE]
meta_df <- meta_df[common_ids, , drop = FALSE]

# ====================================================
### 4. Helper functions ------------------------
source("scripts/rda_helpers.R")
# ====================================================
### 5. Partial RDA analyses ------------------------
# ====================================================

# A. combined environment after conditioning on population structure
res_raw_env <- run_partial_rda(
  Y = cnv_raw,
  meta_df = meta_df,
  explanatory = env_vars,
  conditional = structure_cols,
  label = "Raw CNV ~ climate + soil | structure"
)

res_z_env <- run_partial_rda(
  Y = cnv_z,
  meta_df = meta_df,
  explanatory = env_vars,
  conditional = structure_cols,
  label = "Z-score CNV ~ climate + soil | structure"
)

# B. pure climate after conditioning on soil + structure
res_raw_clim <- if (length(climate_vars) > 0) {
  run_partial_rda(
    Y = cnv_raw,
    meta_df = meta_df,
    explanatory = climate_vars,
    conditional = unique(c(structure_cols, soil_vars)),
    label = "Raw CNV ~ climate | soil + structure"
  )
} else NULL

res_z_clim <- if (length(climate_vars) > 0) {
  run_partial_rda(
    Y = cnv_z,
    meta_df = meta_df,
    explanatory = climate_vars,
    conditional = unique(c(structure_cols, soil_vars)),
    label = "Z-score CNV ~ climate | soil + structure"
  )
} else NULL

# C. pure soil after conditioning on climate + structure
res_raw_soil <- if (length(soil_vars) > 0) {
  run_partial_rda(
    Y = cnv_raw,
    meta_df = meta_df,
    explanatory = soil_vars,
    conditional = unique(c(structure_cols, climate_vars)),
    label = "Raw CNV ~ soil | climate + structure"
  )
} else NULL

res_z_soil <- if (length(soil_vars) > 0) {
  run_partial_rda(
    Y = cnv_z,
    meta_df = meta_df,
    explanatory = soil_vars,
    conditional = unique(c(structure_cols, climate_vars)),
    label = "Z-score CNV ~ soil | climate + structure"
  )
} else NULL

# save core results
save_rda_outputs <- function(res_obj, prefix) {
  if (is.null(res_obj)) return(NULL)
  write.csv(res_obj$overall, file.path(out_dir, paste0(prefix, "_overall.csv")), row.names = TRUE)
  write.csv(res_obj$terms,   file.path(out_dir, paste0(prefix, "_terms.csv")), row.names = TRUE)
  write.csv(res_obj$r2,      file.path(out_dir, paste0(prefix, "_r2.csv")), row.names = FALSE)
}

save_rda_outputs(res_raw_env,  "rda_raw_env")
save_rda_outputs(res_z_env,    "rda_z_env")
save_rda_outputs(res_raw_clim, "rda_raw_climate")
save_rda_outputs(res_z_clim,   "rda_z_climate")
save_rda_outputs(res_raw_soil, "rda_raw_soil")
save_rda_outputs(res_z_soil,   "rda_z_soil")

# ====================================================
### 6. Variance partitioning ------------------------
# ====================================================

# among climate, soil, structure
vp_raw <- varpart(cnv_raw,
                  meta_df[, climate_vars, drop = FALSE],
                  meta_df[, soil_vars, drop = FALSE],
                  meta_df[, structure_cols, drop = FALSE])

vp_z <- varpart(cnv_z,
                meta_df[, climate_vars, drop = FALSE],
                meta_df[, soil_vars, drop = FALSE],
                meta_df[, structure_cols, drop = FALSE])

vp_raw_df <- extract_varpart_summary(vp_raw, "Raw CNV")
vp_z_df   <- extract_varpart_summary(vp_z, "Z-score CNV")

write.csv(vp_raw_df, file.path(out_dir, "variance_partition_raw.csv"), row.names = FALSE)
write.csv(vp_z_df,   file.path(out_dir, "variance_partition_zscore.csv"), row.names = FALSE)

# explicit model-level summary table
model_summary <- bind_rows(
  res_raw_env$r2 %>% mutate(model = "Raw: climate+soil | structure"),
  res_z_env$r2   %>% mutate(model = "Z-score: climate+soil | structure"),
  if (!is.null(res_raw_clim)) res_raw_clim$r2 %>% mutate(model = "Raw: climate | soil+structure"),
  if (!is.null(res_z_clim))   res_z_clim$r2   %>% mutate(model = "Z-score: climate | soil+structure"),
  if (!is.null(res_raw_soil)) res_raw_soil$r2 %>% mutate(model = "Raw: soil | climate+structure"),
  if (!is.null(res_z_soil))   res_z_soil$r2   %>% mutate(model = "Z-score: soil | climate+structure")
)

write.csv(model_summary, file.path(out_dir, "rda_model_summary.csv"), row.names = FALSE)

# ====================================================
### 7. Figures ------------------------
# ====================================================

# 7.1 compare raw vs z-score for key partial RDA models
p_model_comp <- model_summary %>%
  pivot_longer(c(r2, adj_r2), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = fct_reorder(model, value), y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c("adj_r2" = "#377eb8","r2" = "#4daf4a"),
                    labels = c("Adjusted R²","R²")) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Explained variance",
       fill = NULL,
       title = "Comparison of partial RDA models")

ggsave(file.path(out_dir, "Fig_1_partial_rda_model_comparison.pdf"),
       p_model_comp, width = 8.2, height = 5.6)

# 7.2 term-wise environmental contributions for combined-environment models
prep_term_plot <- function(tt, label) {
  tt <- as.data.frame(tt)
  tt$term <- rownames(tt)
  rownames(tt) <- NULL
  
  pcol <- grep("Pr", names(tt), value = TRUE)[1]
  vcol <- if ("Variance" %in% names(tt)) {
    "Variance"
  } else if ("ChiSquare" %in% names(tt)) {
    "ChiSquare"
  } else {
    names(tt)[2]
  }
  
  tt %>%
    filter(term != "Residual") %>%
    mutate(
      model = label,
      significant = .data[[pcol]] < 0.05,
      variance = .data[[vcol]]
    ) %>%
    select(term, model, significant, variance)
}

term_df <- bind_rows(
  prep_term_plot(res_raw_env$terms, "Raw CNV"),
  prep_term_plot(res_z_env$terms, "Z-score CNV")
)

write.csv(term_df, file.path(out_dir, "rda_env_term_plot_table.csv"), row.names = FALSE)

all_terms <- unique(term_df$term)
bio_terms <- sort(all_terms[str_detect(all_terms, "^bio")])
other_terms <- sort(all_terms[!str_detect(all_terms, "^bio")])

term_order <- c(bio_terms, other_terms)
term_df <- term_df %>%
  mutate(term = factor(term, levels = term_order))

p_terms <- ggplot(term_df,
                  aes(x = term, y = variance, fill = significant)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~model, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "gray70")) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Constrained variance",
       fill = "p < 0.05",
       title = "Environmental terms contributing to CNV structure")

ggsave(file.path(out_dir, "Fig_2_partial_rda_env_term_variance.pdf"),
       p_terms, width = 5.5, height = 3)

# 7.3 ordinations
at_meta<-read.csv("data/231/AT/genDist/AT_gen_div_covar_svd_v2.csv")
at_meta<-data.frame(Sample=at_meta$Population,at_meta[,-c(1,2)])

pop_col_use <- if ("pop_name" %in% names(res_raw_env$sites)) {
  "pop_name"
} else if ("Population" %in% names(res_raw_env$sites)) {
  "Population"
} else NULL

p_ord_raw <- plot_rda_ordination(res_raw_env, pop_col = pop_col_use, max_arrows = 8, clusters=at_meta)
p_ord_z   <- plot_rda_ordination(res_z_env,   pop_col = pop_col_use, max_arrows = 8, clusters=at_meta)

ggsave(file.path(out_dir, "Fig_3_partial_rda_ordination_raw.pdf"),
       p_ord_raw, width = 7.5, height = 6)
ggsave(file.path(out_dir, "Fig_4_partial_rda_ordination_zscore.pdf"),
       p_ord_z, width = 7.5, height = 4)
ggsave(file.path(out_dir, "Fig_5_partial_rda_ordination_comparison.pdf"),
       p_ord_raw + p_ord_z, width = 14, height = 6)

# 7.4 variance partitioning diagrams (base plotting)
pdf(file.path(out_dir, "Fig_6_variance_partitioning_raw.pdf"), width = 4, height = 4)
plot(vp_raw, bg = c("#8dd3c7", "#bebada", "#fb8072"),
     Xnames = c("Climate", "Soil", "Structure"),
     id.size = 1)
title("Variance partitioning: Raw CNV")
dev.off()

pdf(file.path(out_dir, "Fig_7_variance_partitioning_zscore.pdf"), width = 5, height = 5)
plot(vp_z, bg = c("#8dd3c7", "#bebada", "#fb8072"),
     Xnames = c("Climate", "Soil", "Structure"),
     id.size = 1)
title("Variance partitioning: Z-score CNV")
dev.off()

# ====================================================
### 8. Concise topline summary ------------------------
# ====================================================

analysis_summary <- tibble(
  analysis = c(
    "Raw: climate+soil | structure",
    "Z-score: climate+soil | structure",
    "Raw: climate | soil+structure",
    "Z-score: climate | soil+structure",
    "Raw: soil | climate+structure",
    "Z-score: soil | climate+structure"
  ),
  adj_r2 = c(
    res_raw_env$r2$adj_r2,
    res_z_env$r2$adj_r2,
    if (!is.null(res_raw_clim)) res_raw_clim$r2$adj_r2 else NA_real_,
    if (!is.null(res_z_clim))   res_z_clim$r2$adj_r2 else NA_real_,
    if (!is.null(res_raw_soil)) res_raw_soil$r2$adj_r2 else NA_real_,
    if (!is.null(res_z_soil))   res_z_soil$r2$adj_r2 else NA_real_
  ),
  r2 = c(
    res_raw_env$r2$r2,
    res_z_env$r2$r2,
    if (!is.null(res_raw_clim)) res_raw_clim$r2$r2 else NA_real_,
    if (!is.null(res_z_clim))   res_z_clim$r2$r2 else NA_real_,
    if (!is.null(res_raw_soil)) res_raw_soil$r2$r2 else NA_real_,
    if (!is.null(res_z_soil))   res_z_soil$r2$r2 else NA_real_
  )
)

write.csv(analysis_summary,
          file.path(out_dir, "analysis_summary_topline.csv"),
          row.names = FALSE)

# ====================================================
## 5. Family-wise association testing ------------------------
#    raw counts only, NB models
#    climate + soil predictors only
#    controlling only for population structure
# ====================================================
################## Z-SCORE CNV–ENV ANALYSIS ##################

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
})

# =========================================================
### 0. INPUT ------------------------
# =========================================================
in_file <- "outputs/AT/eco/main_df_final_CNV_ENV_ready_imputed.csv"
out_dir <- "outputs/AT/eco/cnv_env_univariate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

main_df <- data.table::fread(in_file)
main_df<-main_df[!main_df$accession%in%c("GCA_020911765.2", "GCA_036942435.1", "GCA_946499705.1"),]

# =========================================================
### 1. DEFINE COLUMNS ------------------------
# =========================================================
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

# =========================================================
### 2. BUILD LONG TABLE + FAMILY-WISE Z-SCORE ------------------------
# =========================================================
cnv_long_z <- main_df %>%
  select(all_of(c("assembly", "accession", structure_cols, env_cols, og_cols))) %>%
  pivot_longer(
    cols = all_of(og_cols),
    names_to = "family_id",
    values_to = "copy_number"
  ) %>%
  mutate(copy_number = as.numeric(copy_number))

# keep only variable families
family_stats <- cnv_long_z %>%
  group_by(family_id) %>%
  summarise(
    mean_copy = mean(copy_number, na.rm = TRUE),
    sd_copy = sd(copy_number, na.rm = TRUE),
    n_states = n_distinct(copy_number),
    .groups = "drop"
  )

families_keep <- family_stats %>%
  filter(sd_copy > 0, n_states > 1) %>%
  pull(family_id)

cnv_long_z <- cnv_long_z %>%
  filter(family_id %in% families_keep) %>%
  left_join(family_stats, by = "family_id") %>%
  mutate(
    cnv_z = ifelse(sd_copy > 0, (copy_number - mean_copy) / sd_copy, 0)
  ) %>%
  select(-mean_copy, -sd_copy, -n_states)

# standardize predictors for comparable coefficients
cnv_long_z <- cnv_long_z %>%
  mutate(across(all_of(c(structure_cols, env_cols)), ~ as.numeric(scale(.x))))

write.csv(
  family_stats,
  file.path(out_dir, "family_zscore_baseline_statistics.csv"),
  row.names = FALSE
)

# =========================================================
### 3. HELPER FUNCTIONS ------------------------
# =========================================================
run_family_lm_z <- function(df, env_var, covars = c("PC1", "PC2")) {
  rhs <- paste(c(env_var, covars), collapse = " + ")
  form <- as.formula(paste0("cnv_z ~ ", rhs))
  
  fit <- lm(form, data = df)
  td <- broom::tidy(fit)
  gl <- broom::glance(fit)
  
  tr <- td %>% filter(term == env_var)
  if (nrow(tr) == 0) return(tibble())
  
  tibble(
    estimate = tr$estimate,
    std.error = tr$std.error,
    statistic = tr$statistic,
    p.value = tr$p.value,
    r.squared = gl$r.squared,
    adj.r.squared = gl$adj.r.squared
  )
}

# =========================================================
### 4. FAMILY-WISE SCREENING WITH Z-SCORED CNV ------------------------
# =========================================================
ogList<-data.frame(fread("data/ogList_with_functions_v2.csv"))
function_short<-data.frame(fread("data/og_fam_function_short.tsv",h=T))
function_short<-function_short[match(ogList$Family.ID,function_short$Family.ID),]
ogList$func_short<-function_short$short_title

family_env_z_results <- purrr::map_dfr(env_cols, function(ev) {
  cnv_long_z %>%
    filter(!is.na(.data[[ev]]), !is.na(cnv_z), !is.na(PC1), !is.na(PC2)) %>%
    group_by(family_id) %>%
    group_modify(~ {
      out <- tryCatch(
        run_family_lm_z(.x, env_var = ev, covars = structure_cols),
        error = function(e) tibble()
      )
      out
    }) %>%
    ungroup() %>%
    mutate(env_var = ev)
}) %>%
  group_by(env_var) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write.csv(
  family_env_z_results,
  file.path(out_dir, "familywise_zscore_environment_associations.csv"),
  row.names = FALSE
)

# add family functinal information
og_matched<-ogList[match(family_env_z_results$family_id,ogList$Family.ID),c("func_short","FunctionalCategory")]
family_env_z_results<-data.frame(short_func=og_matched[,1],func_cat=og_matched[,2],family_env_z_results)


family_env_z_summary <- family_env_z_results %>%
  group_by(env_var) %>%
  summarise(
    n_tested = sum(!is.na(p.value)),
    n_sig_p05 = sum(p.value < 0.05, na.rm = TRUE),
    n_sig_fdr05 = sum(q.value < 0.05, na.rm = TRUE),
    median_abs_effect = median(abs(estimate), na.rm = TRUE),
    median_adj_r2 = median(adj.r.squared, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_fdr05), desc(n_sig_p05))

write.csv(
  family_env_z_summary,
  file.path(out_dir, "familywise_zscore_environment_summary.csv"),
  row.names = FALSE
)

# =========================================================
### 5. FIGURE 1: SUMMARY BARPLOT ------------------------
# =========================================================
all_env_vars <- unique(family_env_z_summary$env_var)

bio_vars <- sort(all_env_vars[str_detect(all_env_vars, "^bio")])
other_vars <- sort(all_env_vars[!str_detect(all_env_vars, "^bio")])
var_order <- c(bio_vars, other_vars)
family_env_z_summary_ordered <- family_env_z_summary %>%
  mutate(env_var = factor(env_var, levels = var_order))

p_env_bar_z <- family_env_z_summary_ordered %>%
  pivot_longer(
    cols = c(n_sig_p05, n_sig_fdr05),
    names_to = "metric",
    values_to = "n_families"
  ) %>%
  ggplot(aes(x = env_var, y = n_families, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  coord_flip() +
  scale_fill_manual(
    values = c("n_sig_p05" = "#984ea3", "n_sig_fdr05" = "#4daf4a"),
    labels = c("n_sig_p05" = "p < 0.05", "n_sig_fdr05" = "FDR < 0.05")
  ) +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "Number of families",
    fill = NULL,
    title = "Family-wise CNV–environment associations using z-scored CNV"
  )

p_env_bar_z

ggsave(
  file.path(out_dir, "Fig_1_zscore_familywise_hits.pdf"),
  p_env_bar_z, width = 4, height = 6
)

# =========================================================
### 6. FIGURE 2: HEATMAP OF SIGNIFICANT / STRONG FAMILIES ------------------------
# =========================================================
# use FDR hits if present, otherwise fall back to nominal hits

library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)

heat_df_z <- family_env_z_results %>% filter(!is.na(p.value))

if (sum(heat_df_z$q.value < 0.05, na.rm = TRUE) > 0) {
  heat_df_z <- heat_df_z %>% mutate(sig_class = ifelse(q.value < 0.05, "FDR", "NS"))
  fam_rank <- heat_df_z %>%
    filter(q.value < 0.05) %>%
    group_by(family_id) %>%
    summarise(
      best_sig = min(q.value, na.rm = TRUE),
      max_abs_effect = max(abs(estimate), na.rm = TRUE),
      .groups = "drop"
    )
  heat_subtitle <- "Stars indicate FDR < 0.05"
} else {
  heat_df_z <- heat_df_z %>% mutate(sig_class = ifelse(p.value < 0.05, "p05", "NS"))
  fam_rank <- heat_df_z %>%
    filter(p.value < 0.05) %>%
    group_by(family_id) %>%
    summarise(
      best_sig = min(p.value, na.rm = TRUE),
      max_abs_effect = max(abs(estimate), na.rm = TRUE),
      .groups = "drop"
    )
  heat_subtitle <- "Stars indicate nominal p < 0.05"
}

top_fams_heat <- fam_rank %>%
  arrange(best_sig, desc(max_abs_effect)) %>%
  slice(1:min(50, n())) %>%
  pull(family_id)

func_levels <- sort(unique(heat_df_z$func_cat))
func_cols <- setNames(hue_pal(l = 65, c = 100)(length(func_levels)), func_levels)

fam_meta <- heat_df_z %>%
  distinct(family_id, short_func, func_cat) %>%
  filter(family_id %in% top_fams_heat) %>%
  left_join(fam_rank, by = "family_id") %>%
  arrange(func_cat, best_sig, desc(max_abs_effect), short_func, family_id) %>%
  mutate(
    family_id = factor(family_id, levels = rev(family_id)),
    y_lab = paste0("<span style='color:", func_cols[func_cat], ";'>", short_func, "</span>")
  )

heat_df_plot <- heat_df_z %>%
  filter(family_id %in% top_fams_heat) %>%
  left_join(fam_meta %>% select(family_id, short_func, func_cat, y_lab), by = "family_id") %>%
  mutate(
    env_var = factor(env_var, levels = unique(heat_df_z$env_var)),
    family_id = factor(family_id, levels = levels(fam_meta$family_id))
  )

n_env <- nlevels(heat_df_plot$env_var)

right_lab_df <- fam_meta %>%
  mutate(xpos = n_env + 0.35)

legend_df <- data.frame(
  func_cat = factor(func_levels, levels = func_levels),
  x = Inf,
  y = Inf
)

p_heat_z <- ggplot(heat_df_plot, aes(x = as.numeric(env_var), y = family_id, fill = estimate)) +
  geom_tile(color = "white") +
  geom_point(
    data = subset(heat_df_plot, sig_class != "NS"),
    shape = 8, size = 1.7, color = "black"
  ) +
  geom_text(
    data = right_lab_df,
    aes(x = xpos, y = family_id, label = as.character(family_id)),
    inherit.aes = FALSE, hjust = 0, size = 2.8
  ) +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, color = func_cat),
    inherit.aes = FALSE, alpha = 0, size = 3, show.legend = TRUE
  ) +
  scale_x_continuous(
    breaks = seq_len(n_env),
    labels = levels(heat_df_plot$env_var),
    expand = expansion(mult = c(0, 0), add = c(0, 1.8))
  ) +
  scale_y_discrete(
    labels = setNames(as.character(fam_meta$y_lab), as.character(fam_meta$family_id))
  ) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    name = "Effect (Z-score)"
  ) +
  scale_color_manual(values = func_cols, name = "Functional category") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = ggtext::element_markdown(),
    plot.margin = margin(5.5, 140, 5.5, 5.5),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    color = guide_legend(order = 2, override.aes = list(alpha = 1, size = 4))
  ) +
  labs(
    x = NULL,
    y = "Gene family",
    title = "Family × environment associations using z-scored CNV",
    subtitle = heat_subtitle
  )

ggsave(
  file.path(out_dir, "Fig_2_zscore_family_environment_heatmap_v2.pdf"),
  p_heat_z, width = 8, height = 10
)


topfams<-heat_df_plot[heat_df_plot$sig_class=="FDR",]
topfams<-topfams[!duplicated(topfams$family_id),]
topfams<-topfams[order(topfams$func_cat.y),]


write.csv(
  topfams,
  file.path(out_dir,"AT_top_families_summary.csv"),
  row.names = F
)

# =========================================================
### 7. FIGURE 3: REACTION-NORM STYLE PLOTS ------------------------
# =========================================================
# choose strongest env variables
top_envs <- family_env_z_summary %>%
  slice(1:min(3, n())) %>%
  pull(env_var)

# choose strongest family × env combinations
top_hits_z <- family_env_z_results %>%
  filter(env_var %in% top_envs) %>%
  arrange(q.value, p.value) %>%
  slice(1:min(12, n())) %>%
  select(family_id, env_var) %>%
  distinct()

fam_lookup <- family_env_z_results %>%
  distinct(family_id, short_func, func_cat)

reaction_plot_list <- purrr::pmap(
  list(top_hits_z$family_id, top_hits_z$env_var),
  function(fid, ev) {
    
    meta <- fam_lookup %>% filter(family_id == fid) %>% slice(1)
    
    dat <- cnv_long_z %>%
      filter(family_id == fid) %>%
      filter(!is.na(.data[[ev]]), !is.na(cnv_z), !is.na(PC1), !is.na(PC2))
    
    ggplot(dat, aes(x = .data[[ev]], y = cnv_z)) +
      geom_point(size = 2, alpha = 0.75) +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      theme_classic(base_size = 10) +
      labs(
        title = paste0(meta$short_func, " (", fid, ")"),
        subtitle = meta$func_cat,
        x = ev,
        y = "CNV (z-score)"
      )
  }
)

if (length(reaction_plot_list) > 0) {
  pdf(file.path(out_dir, "Fig_3_zscore_reaction_norms_top_hits.pdf"), width = 10, height = 12)
  print(patchwork::wrap_plots(reaction_plot_list, ncol = 2))
  dev.off()
}

# =========================================================
### TOPLINE COMPARISON TABLE ------------------------
# =========================================================
analysis_summary_z <- family_env_z_summary %>%
  transmute(
    env_var,
    n_tested,
    n_sig_p05,
    n_sig_fdr05,
    median_abs_effect,
    median_adj_r2
  )

write.csv(
  analysis_summary_z,
  file.path(out_dir, "analysis_summary_zscore_familywise.csv"),
  row.names = FALSE
)



#################################################################################
#################################################################################
#################################################################################
#################################################################################


# LYRATA ------------------------------------------------
rm(list=ls())


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(vegan)
  library(caret)
  library(purrr)
  library(broom)
  library(car)
  library(glmmTMB)
  library(broom)
  library(scales)
  library(ggrepel)
  library(patchwork)
})

setwd("/Users/piyalkaru/Desktop/DDORF/Ann/231ogs/cnv_geneFam")

# All data joined earlier
main_df <- data.table::fread("data/231/AL/eco/AL_main_data_for_eco_v1.1.tsv", h = TRUE)
# drop ploblematic assembly NT1
main_df<-main_df[-which(main_df$Population_name=="NT1"),]

out_dir <- "outputs/AL/eco"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1. Define columns ------------------------
## =========================================================

id_cols <- c("Population_name","latitude", "longitude")

# keep structure as available
structure_cols <- c("PC1", "PC2")

tech_cols <- c("busco_complete_pct",
               "contig_n50_mbp",
               "total_length_bp")

# TE as genomic/ecological covariate
te_cols <- c("te_int_frac")

og_cols   <- grep("^OG", names(main_df), value = TRUE)
bio_cols  <- grep("^bio\\d+$", names(main_df), value = TRUE)
soil_cols <- grep("250m", names(main_df), value = TRUE)

# full predictor pool for CNV-ENV analysis
env_cols <- c(bio_cols, soil_cols, te_cols)
env_cols <- unique(env_cols)

## optional: convert assembly size to Mbp
main_df <- main_df %>%
  mutate(total_length_mbp = total_length_bp / 1e6)

tech_cols2 <- c("busco_complete_pct",
                "contig_n50_mbp",
                "total_length_mbp")

## =========================================================
## 2. Technical covariate screening ------------------------
## =========================================================

tech_df <- main_df %>%
  select(all_of(c("Population_name", tech_cols2, structure_cols))) %>%
  distinct()

cor_mat <- cor(tech_df[, ..tech_cols2], use = "pairwise.complete.obs", method = "spearman")

write.csv(cor_mat,
          file.path(out_dir, "technical_covariate_correlations.csv"),
          row.names = TRUE)

cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("Var1", "Var2", "rho")

p_cor <- ggplot(cor_long, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", midpoint = 0) +
  theme_classic(base_size = 12) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Spearman rho",
       title = "Correlation among technical covariates")

ggsave(file.path(out_dir, "Fig_S1_technical_covariate_correlation_heatmap.pdf"),
       p_cor, width = 6.2, height = 5.2)

highly_cor <- caret::findCorrelation(cor_mat, cutoff = 0.8, names = TRUE)
tech_selected <- setdiff(tech_cols2, highly_cor)
if (length(tech_selected) == 0) tech_selected <- tech_cols2
message("Technical covariates retained: ",
        paste(tech_selected, collapse = ", "))

## =========================================================
## 3. Reshape CNV and keep variable families ------------------------
## =========================================================

cnv_long <- main_df %>%
  select(all_of(c("Population_name", structure_cols, tech_cols2, og_cols))) %>%
  pivot_longer(cols = all_of(og_cols),
               names_to = "family_id",
               values_to = "copy_number") %>%
  mutate(copy_number = as.numeric(copy_number))

family_filter <- cnv_long %>%
  group_by(family_id) %>%
  summarise(
    n = n(),
    n_unique = n_distinct(copy_number),
    sd_copy = sd(copy_number, na.rm = TRUE),
    mean_copy = mean(copy_number, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_unique > 1, sd_copy > 0)

cnv_long_filt <- cnv_long %>%
  filter(family_id %in% family_filter$family_id)

## =========================================================
## 4. Standardize copy number within family ------------------------
## =========================================================

cnv_long_filt <- cnv_long_filt %>%
  group_by(family_id) %>%
  mutate(
    cnv_dev = copy_number - median(copy_number, na.rm = TRUE),
    cnv_z = {
      s <- sd(copy_number, na.rm = TRUE)
      if (is.na(s) || s == 0) rep(0, dplyr::n())
      else (copy_number - mean(copy_number, na.rm = TRUE)) / s
    }
  ) %>%
  ungroup()

## =========================================================
## 5. Multivariate screening: ------------------------
##    Does technical variation explain CNV matrix
##    after controlling for structure?
## =========================================================

cnv_mat <- cnv_long_filt %>%
  select(Population_name, family_id, cnv_z) %>%
  pivot_wider(names_from = family_id, values_from = cnv_z) %>%
  as.data.frame()

rownames(cnv_mat) <- cnv_mat$Population_name
cnv_mat$Population_name <- NULL

cnv_mat <- cnv_mat[, colSums(is.na(cnv_mat)) == 0, drop = FALSE]
cnv_mat <- cnv_mat[, apply(cnv_mat, 2, sd) > 0, drop = FALSE]

meta_df <- main_df %>%
  select(all_of(c("Population_name", structure_cols, tech_selected))) %>%
  distinct() %>%
  as.data.frame()

rownames(meta_df) <- meta_df$Population_name
meta_df$Population_name <- NULL

common_ids <- intersect(rownames(cnv_mat), rownames(meta_df))
cnv_mat <- cnv_mat[common_ids, , drop = FALSE]
meta_df <- meta_df[common_ids, , drop = FALSE]

tech_formula <- paste(tech_selected, collapse = " + ")
rda_formula <- as.formula(paste0("cnv_mat ~ ", tech_formula, " + Condition(PC1 + PC2)"))

rda_tech <- vegan::rda(rda_formula, data = meta_df)

rda_overall <- anova.cca(rda_tech, permutations = 999)
rda_terms   <- anova.cca(rda_tech, by = "term", permutations = 999)
rda_rsq     <- RsquareAdj(rda_tech)

write.csv(as.data.frame(rda_overall),
          file.path(out_dir, "technical_rda_overall.csv"),
          row.names = TRUE)
write.csv(as.data.frame(rda_terms),
          file.path(out_dir, "technical_rda_term_tests.csv"),
          row.names = TRUE)

rda_term_df <- as.data.frame(rda_terms)
rda_term_df$term <- rownames(rda_term_df)
rownames(rda_term_df) <- NULL
pcol_rda <- grep("Pr", names(rda_term_df), value = TRUE)[1]

vcol_rda <- if ("Variance" %in% names(rda_term_df)) {
  "Variance"
} else if ("ChiSquare" %in% names(rda_term_df)) {
  "ChiSquare"
} else {
  names(rda_term_df)[2]
}

p_rda <- ggplot(rda_term_df %>% filter(term != "Residual"),
                aes(x = reorder(term, .data[[vcol_rda]]),
                    y = .data[[vcol_rda]],
                    fill = .data[[pcol_rda]] < 0.05)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "gray70")) +
  theme_classic(base_size = 12) +
  labs(x = NULL,
       y = "Explained constrained variance",
       fill = "p < 0.05",
       title = "Multivariate association between CNV and technical covariates",
       subtitle = paste0("Partial RDA controlling for PC1 and PC2; adjusted R² = ",
                         round(rda_rsq$adj.r.squared, 3)))

ggsave(file.path(out_dir, "Fig_S2_multivariate_CNV_vs_technical_covariates.pdf"),
       p_rda, width = 7, height = 4.8)

## =========================================================
## 6. Family-wise technical screening ------------------------
## =========================================================

run_family_lm <- function(df, tech_var) {
  fit <- lm(as.formula(paste0("cnv_z ~ ", tech_var, " + PC1 + PC2")), data = df)
  td <- broom::tidy(fit)
  gl <- broom::glance(fit)
  term_row <- td %>% filter(term == tech_var)
  
  tibble(
    estimate = term_row$estimate,
    std.error = term_row$std.error,
    statistic = term_row$statistic,
    p.value = term_row$p.value,
    r.squared = gl$r.squared,
    adj.r.squared = gl$adj.r.squared
  )
}

family_results <- map_dfr(tech_selected, function(tv) {
  cnv_long_filt %>%
    group_by(family_id) %>%
    group_modify(~ {
      out <- tryCatch(run_family_lm(.x, tv), error = function(e) NULL)
      if (is.null(out)) {
        return(tibble(
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          r.squared = NA_real_,
          adj.r.squared = NA_real_
        ))
      }
      out
    }) %>%
    ungroup() %>%
    mutate(technical_var = tv)
}) %>%
  group_by(technical_var) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write.csv(family_results,
          file.path(out_dir, "familywise_technical_associations.csv"),
          row.names = FALSE)

family_summary <- family_results %>%
  group_by(technical_var) %>%
  summarise(
    n_tested = sum(!is.na(p.value)),
    n_sig_p05 = sum(p.value < 0.05, na.rm = TRUE),
    n_sig_fdr05 = sum(q.value < 0.05, na.rm = TRUE),
    median_abs_effect = median(abs(estimate), na.rm = TRUE),
    median_adj_r2 = median(adj.r.squared, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_fdr05))

write.csv(family_summary,
          file.path(out_dir, "familywise_technical_summary.csv"),
          row.names = FALSE)

p_bar <- family_summary %>%
  mutate(technical_var = fct_reorder(technical_var, n_sig_fdr05)) %>%
  ggplot(aes(x = technical_var, y = n_sig_fdr05)) +
  geom_col(fill = "#377eb8") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL,
       y = "Number of families (FDR < 0.05)",
       title = "Family-wise associations with Population_name-quality covariates")

ggsave(file.path(out_dir, "Fig_S3_number_of_families_associated_with_technical_covariates.pdf"),
       p_bar, width = 6.5, height = 4.5)

heat_df <- family_results %>%
  filter(q.value < 0.05) %>%
  mutate(
    technical_var = factor(technical_var, levels = tech_selected),
    family_id = fct_reorder(family_id, abs(estimate), .fun = max, .desc = TRUE)
  )

if (nrow(heat_df) > 0) {
  p_heat <- ggplot(heat_df,
                   aes(x = technical_var, y = family_id, fill = estimate)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
    theme_classic(base_size = 11) +
    labs(x = NULL, y = "Gene family",
         fill = "Effect size",
         title = "Effect sizes of significant technical covariates on family CNV")
  
  # ggsave(file.path(out_dir, "Fig_S4_familywise_effect_heatmap_technical_covariates.pdf"),
  #        p_heat, width = 6.8, height = 8.5)
}

flagged_families <- family_results %>%
  group_by(family_id) %>%
  summarise(
    any_sig_fdr05 = any(q.value < 0.05, na.rm = TRUE),
    n_sig_covars = sum(q.value < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(any_sig_fdr05)

write.csv(flagged_families,
          file.path(out_dir, "flagged_families_technical_confounds.csv"),
          row.names = FALSE)

## =========================================================
## 6b. TE content screening: ------------------------
##     Does CNV associate with TE content after controlling
##     for population structure?
## =========================================================

if (!"te_int_frac" %in% names(main_df)) {
  warning("te_int_frac not found in main_df; skipping TE screening")
} else {
  
  ## -------------------------------------------------------
  ## A. Multivariate TE ~ CNV matrix (partial RDA)
  ## -------------------------------------------------------
  meta_te <- main_df %>%
    select(Population_name, PC1, PC2, te_int_frac) %>%
    distinct() %>%
    as.data.frame()
  
  rownames(meta_te) <- meta_te$Population_name
  meta_te$Population_name <- NULL
  
  common_ids_te <- intersect(rownames(cnv_mat), rownames(meta_te))
  cnv_mat_te <- cnv_mat[common_ids_te, , drop = FALSE]
  meta_te <- meta_te[common_ids_te, , drop = FALSE]
  
  # remove rows with missing TE or PCs
  keep_te <- complete.cases(meta_te[, c("PC1", "PC2", "te_int_frac")])
  cnv_mat_te <- cnv_mat_te[keep_te, , drop = FALSE]
  meta_te <- meta_te[keep_te, , drop = FALSE]
  
  if (nrow(meta_te) > 5) {
    rda_te <- vegan::rda(cnv_mat_te ~ te_int_frac + Condition(PC1 + PC2), data = meta_te)
    
    rda_te_overall <- anova.cca(rda_te, permutations = 999)
    rda_te_term <- anova.cca(rda_te, by = "term", permutations = 999)
    rda_te_rsq <- RsquareAdj(rda_te)
    
    write.csv(as.data.frame(rda_te_overall),
              file.path(out_dir, "te_rda_overall.csv"),
              row.names = TRUE)
    
    write.csv(as.data.frame(rda_te_term),
              file.path(out_dir, "te_rda_term_test.csv"),
              row.names = TRUE)
    
    te_term_df <- as.data.frame(rda_te_term)
    te_term_df$term <- rownames(te_term_df)
    rownames(te_term_df) <- NULL
    
    pcol_te <- grep("Pr", names(te_term_df), value = TRUE)[1]
    vcol_te <- if ("Variance" %in% names(te_term_df)) {
      "Variance"
    } else if ("ChiSquare" %in% names(te_term_df)) {
      "ChiSquare"
    } else {
      names(te_term_df)[2]
    }
    
    p_te_rda <- ggplot(te_term_df %>% filter(term != "Residual"),
                       aes(x = reorder(term, .data[[vcol_te]]),
                           y = .data[[vcol_te]],
                           fill = .data[[pcol_te]] < 0.05)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "gray70")) +
      theme_classic(base_size = 12) +
      labs(x = NULL,
           y = "Explained constrained variance",
           fill = "p < 0.05",
           title = "Multivariate association between CNV and TE content",
           subtitle = paste0("Partial RDA controlling for PC1 and PC2; adjusted R² = ",
                             round(rda_te_rsq$adj.r.squared, 3)))
    
    ggsave(file.path(out_dir, "Fig_S4_multivariate_CNV_vs_TE_content.pdf"),
           p_te_rda, width = 6.8, height = 4.5)
  }
  
  ## -------------------------------------------------------
  ## B. Family-wise screening: ------------------------
  ##    cnv_z ~ te_int_frac + PC1 + PC2
  ## -------------------------------------------------------
  run_family_lm_te <- function(df) {
    fit <- lm(cnv_z ~ te_int_frac + PC1 + PC2, data = df)
    td <- broom::tidy(fit)
    gl <- broom::glance(fit)
    term_row <- td %>% filter(term == "te_int_frac")
    
    tibble(
      estimate = term_row$estimate,
      std.error = term_row$std.error,
      statistic = term_row$statistic,
      p.value = term_row$p.value,
      r.squared = gl$r.squared,
      adj.r.squared = gl$adj.r.squared
    )
  }
  
  cnv_long_te <- cnv_long_filt %>%
    left_join(
      main_df %>%
        select(Population_name, te_int_frac) %>%
        distinct(),
      by = "Population_name"
    ) %>%
    filter(!is.na(te_int_frac), !is.na(PC1), !is.na(PC2))
  
  family_results_te <- cnv_long_te %>%
    group_by(family_id) %>%
    group_modify(~ {
      out <- tryCatch(run_family_lm_te(.x), error = function(e) NULL)
      if (is.null(out)) {
        return(tibble(
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          r.squared = NA_real_,
          adj.r.squared = NA_real_
        ))
      }
      out
    }) %>%
    ungroup() %>%
    mutate(q.value = p.adjust(p.value, method = "BH"))
  
  write.csv(family_results_te,
            file.path(out_dir, "familywise_TE_associations.csv"),
            row.names = FALSE)
  
  family_summary_te <- family_results_te %>%
    summarise(
      n_tested = sum(!is.na(p.value)),
      n_sig_p05 = sum(p.value < 0.05, na.rm = TRUE),
      n_sig_fdr05 = sum(q.value < 0.05, na.rm = TRUE),
      median_abs_effect = median(abs(estimate), na.rm = TRUE),
      median_adj_r2 = median(adj.r.squared, na.rm = TRUE)
    )
  
  write.csv(family_summary_te,
            file.path(out_dir, "familywise_TE_summary.csv"),
            row.names = FALSE)
  
  p_te_bar <- family_summary_te %>%
    pivot_longer(cols = c(n_sig_p05, n_sig_fdr05),
                 names_to = "metric",
                 values_to = "n_families") %>%
    ggplot(aes(x = metric, y = n_families)) +
    geom_col(fill = "#984ea3") +
    theme_classic(base_size = 12) +
    labs(x = NULL,
         y = "Number of families",
         title = "Family-wise associations between CNV and TE content") +
    scale_x_discrete(labels = c("n_sig_p05" = "p < 0.05",
                                "n_sig_fdr05" = "FDR < 0.05"))
  
  ggsave(file.path(out_dir, "Fig_S5_number_of_families_associated_with_TE_content.pdf"),
         p_te_bar, width = 4.8, height = 4.2)
  
  heat_df_te <- family_results_te %>%
    filter(q.value < 0.05) %>%
    mutate(family_id = fct_reorder(family_id, abs(estimate), .desc = TRUE))
  
  if (nrow(heat_df_te) > 0) {
    p_heat_te <- ggplot(heat_df_te,
                        aes(x = "te_int_frac", y = family_id, fill = estimate)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
      theme_classic(base_size = 11) +
      labs(x = NULL, y = "Gene family",
           fill = "Effect size",
           title = "Effect sizes of TE content on family CNV")
    
    ggsave(file.path(out_dir, "Fig_S6_familywise_effect_heatmap_TE_content.pdf"),
           p_heat_te, width = 5.5, height = 8)
  }
  
  flagged_families_te <- family_results_te %>%
    filter(q.value < 0.05)
  
  write.csv(flagged_families_te,
            file.path(out_dir, "flagged_families_TE_confounds.csv"),
            row.names = FALSE)
}

## =========================================================
## 7. Filter gene families before CNV-ENV analysis ------------------------
## =========================================================

family_stats <- main_df %>%
  select(all_of(c("Population_name", og_cols))) %>%
  pivot_longer(cols = all_of(og_cols),
               names_to = "family_id",
               values_to = "copy_number") %>%
  group_by(family_id) %>%
  summarise(
    n_individuals = n(),
    mean_copy = mean(copy_number, na.rm = TRUE),
    sd_copy = sd(copy_number, na.rm = TRUE),
    n_states = n_distinct(copy_number),
    n_nonzero = sum(copy_number > 0, na.rm = TRUE),
    zero_fraction = mean(copy_number == 0, na.rm = TRUE),
    min_copy = min(copy_number, na.rm = TRUE),
    max_copy = max(copy_number, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(family_stats,
          file.path(out_dir, "family_filtering_statistics.csv"),
          row.names = FALSE)

min_nonzero <- 5
min_states <- 2
max_zero_fraction <- 0.95

family_stats <- family_stats %>%
  mutate(
    keep_sd = sd_copy > 0,
    keep_nonzero = n_nonzero >= min_nonzero,
    keep_states = n_states >= min_states,
    keep_zero = zero_fraction <= max_zero_fraction,
    keep_final = keep_sd & keep_nonzero & keep_states & keep_zero
  )

families_kept <- family_stats %>% filter(keep_final) %>% pull(family_id)
families_removed <- setdiff(og_cols, families_kept)

write.csv(family_stats,
          file.path(out_dir, "family_filtering_decisions.csv"),
          row.names = FALSE)

write.table(families_kept,
            file.path(out_dir, "families_retained.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(families_removed,
            file.path(out_dir, "families_removed.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

filter_summary <- tibble(
  criterion = c("Invariant (sd = 0)",
                paste0("< ", min_nonzero, " nonzero individuals"),
                paste0("< ", min_states, " copy-number states"),
                paste0("Zero fraction > ", max_zero_fraction),
                "Retained"),
  n_families = c(
    sum(!family_stats$keep_sd),
    sum(!family_stats$keep_nonzero),
    sum(!family_stats$keep_states),
    sum(!family_stats$keep_zero),
    sum(family_stats$keep_final)
  )
)

p_filter <- ggplot(filter_summary,
                   aes(x = fct_reorder(criterion, n_families), y = n_families)) +
  geom_col(fill = "#377eb8") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(x = NULL, y = "Number of gene families",
       title = "Filtering of gene families prior to CNV–environment analysis")

ggsave(file.path(out_dir, "Fig_S5_family_filtering_summary.pdf"),
       p_filter, width = 7, height = 4.5)

## =========================================================
## 8. Environmental / soil / TE variable assessment ------------------------
## =========================================================

env_raw <- main_df %>%
  select(all_of(c("Population_name", env_cols))) %>%
  distinct()

# 8.1 Missingness assessment
env_missing <- env_raw %>%
  summarise(across(all_of(env_cols), ~ mean(is.na(.x)))) %>%
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "missing_fraction") %>%
  arrange(desc(missing_fraction))

write.csv(env_missing,
          file.path(out_dir, "environment_missingness_summary.csv"),
          row.names = FALSE)

p_missing <- env_missing %>%
  mutate(variable = fct_reorder(variable, missing_fraction)) %>%
  ggplot(aes(x = variable, y = missing_fraction)) +
  geom_col(fill = "#8c6bb1") +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Fraction missing",
       title = "Missingness across climate, soil, and TE variables")

ggsave(file.path(out_dir, "Fig_S6_environment_missingness.pdf"),
       p_missing, width = 7, height = 6)

# threshold for environmental missingness
max_env_missing <- 0.25

env_keep_missing <- env_missing %>%
  filter(missing_fraction <= max_env_missing) %>%
  pull(variable)

env_drop_missing <- setdiff(env_cols, env_keep_missing)

# 8.2 Build reduced env table
env_df <- env_raw %>%
  select(all_of(c("Population_name", env_keep_missing)))

# 8.3 Median imputation for remaining NAs
env_imputed <- env_df %>%
  mutate(across(-Population_name, ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)))

# 8.4 Remove near-zero variance variables
nzv_idx <- nearZeroVar(env_imputed %>% select(-Population_name))
env_nzv_removed <- names(env_imputed %>% select(-Population_name))[nzv_idx]

if (length(nzv_idx) > 0) {
  env_df2 <- env_imputed %>% select(-Population_name) %>% .[, -nzv_idx, drop = FALSE]
} else {
  env_df2 <- env_imputed %>% select(-Population_name)
}

# 8.5 Correlation matrix before collinearity filtering
cor_before <- cor(env_df2, use = "pairwise.complete.obs", method = "spearman")
write.csv(cor_before,
          file.path(out_dir, "environment_correlations_before_filtering.csv"),
          row.names = TRUE)

cor_before_long <- as.data.frame(as.table(cor_before))
colnames(cor_before_long) <- c("Var1", "Var2", "rho")

p_cor_before <- ggplot(cor_before_long, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", midpoint = 0) +
  theme_classic(base_size = 10) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Spearman rho",
       title = "Environmental correlation matrix before filtering") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(out_dir, "Fig_S7_environment_correlation_before.pdf"),
       p_cor_before, width = 8, height = 7)

# 8.6 Remove highly correlated variables
cor_cutoff <- 0.8
high_cor_idx <- findCorrelation(cor_before, cutoff = cor_cutoff, names = FALSE)

if (length(high_cor_idx) > 0) {
  cor_removed <- names(env_df2)[high_cor_idx]
  keep_cols <- setdiff(names(env_df2), cor_removed)
  env_df3 <- env_df2[, ..keep_cols]
} else {
  cor_removed <- character(0)
  env_df3 <- env_df2
}

# 8.7 Iterative VIF filtering
vif_cutoff <- 6

calc_vif_manual <- function(df) {
  vars <- names(df)
  sapply(vars, function(v) {
    others <- setdiff(vars, v)
    if (length(others) == 0) return(1)
    f <- as.formula(paste(v, "~", paste(others, collapse = " + ")))
    r2 <- summary(lm(f, data = df))$r.squared
    1 / (1 - r2)
  })
}

env_df4 <- data.frame(env_df3)
vif_removed <- character(0)

repeat {
  vif_vals <- calc_vif_manual(env_df4)
  max_vif <- max(vif_vals, na.rm = TRUE)
  if (max_vif <= vif_cutoff || ncol(env_df4) <= 1) break
  drop_var <- names(which.max(vif_vals))
  vif_removed <- c(vif_removed, drop_var)
  env_df4 <- env_df4[, setdiff(names(env_df4), drop_var), drop = FALSE]
}

final_env_vars <- names(env_df4)

env_filter_summary <- tibble(
  step = c("Initial variables",
           paste0("Removed missingness > ", max_env_missing),
           "Removed near-zero variance",
           paste0("Removed high correlation (|rho| > ", cor_cutoff, ")"),
           paste0("Removed by VIF > ", vif_cutoff),
           "Final retained variables"),
  n_variables = c(length(env_cols),
                  length(env_drop_missing),
                  length(env_nzv_removed),
                  length(cor_removed),
                  length(vif_removed),
                  length(final_env_vars)),
  variables = c(
    paste(env_cols, collapse = ", "),
    paste(env_drop_missing, collapse = ", "),
    paste(env_nzv_removed, collapse = ", "),
    paste(cor_removed, collapse = ", "),
    paste(vif_removed, collapse = ", "),
    paste(final_env_vars, collapse = ", ")
  )
)

write.csv(env_filter_summary,
          file.path(out_dir, "environment_filtering_summary.csv"),
          row.names = FALSE)

write.table(final_env_vars,
            file.path(out_dir, "environment_variables_retained.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

cor_after <- cor(env_df4, use = "pairwise.complete.obs", method = "spearman")
write.csv(cor_after,
          file.path(out_dir, "environment_correlations_after_filtering.csv"),
          row.names = TRUE)

cor_after_long <- as.data.frame(as.table(cor_after))
colnames(cor_after_long) <- c("Var1", "Var2", "rho")

p_cor_after <- ggplot(cor_after_long, aes(Var1, Var2, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", mid = "white", midpoint = 0) +
  theme_classic(base_size = 10) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "Spearman rho",
       title = "Environmental correlation matrix after filtering") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(file.path(out_dir, "Fig_S8_environment_correlation_after.pdf"),
       p_cor_after, width = 7, height = 6)

## =========================================================
## 9. Final datasets for CNV-ENV analysis ------------------------
## =========================================================

# raw table (retains NAs)
main_df_final_raw <- main_df %>%
  select(all_of(c(id_cols, structure_cols, tech_cols2, final_env_vars, families_kept)))

write.csv(main_df_final_raw,
          file.path(out_dir, "main_df_final_CNV_ENV_ready_raw.csv"),
          row.names = FALSE)

# imputed table for multivariate / regression workflows
env_imputed_final <- env_imputed %>%
  select(all_of(c("Population_name", final_env_vars)))

main_df_final_imputed <- main_df %>%
  select(all_of(c(id_cols, structure_cols, tech_cols2, families_kept))) %>%
  left_join(env_imputed_final, by = "Population_name")

write.csv(main_df_final_imputed,
          file.path(out_dir, "main_df_final_CNV_ENV_ready_imputed.csv"),
          row.names = FALSE)

## also save a compact variable-class table
variable_classes <- tibble(
  variable = c(tech_selected, final_env_vars),
  class = c(rep("technical", length(tech_selected)),
            ifelse(final_env_vars %in% bio_cols, "climate",
                   ifelse(final_env_vars %in% soil_cols, "soil", "TE")))
)

write.csv(variable_classes,
          file.path(out_dir, "variable_classes_for_CNV_ENV.csv"),
          row.names = FALSE)




## CNV–ENV PARTIAL RDA ##################

# ====================================================
### 0. Input ------------------------
# ====================================================

in_file <- "outputs/AL/eco/main_df_final_CNV_ENV_ready_imputed.csv"
varclass_file <- "outputs/AL/eco/variable_classes_for_CNV_ENV.csv"
out_dir <- "outputs/AL/eco/cnv_env_rda"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

main_df <- data.frame(fread(in_file))
var_classes <- fread(varclass_file)

#tetraploids
tet<-c("WS1","al1","AL27","PU6")


# ====================================================
### 1. Define variable sets ------------------------
# ====================================================

id_cols <- intersect(c("Population_name","latitude", "longitude"), names(main_df))

structure_cols <- intersect(c("PC1", "PC2"), names(main_df))
if (length(structure_cols) < 2) stop("PC1 and/or PC2 missing")

og_cols <- grep("^OG", names(main_df), value = TRUE)
if (length(og_cols) == 0) stop("No OG columns found")

# remove TE and technical variables from analysis
climate_vars <- var_classes %>%
  filter(class == "climate") %>%
  pull(variable) %>%
  intersect(names(main_df))

soil_vars <- var_classes %>%
  filter(class == "soil") %>%
  pull(variable) %>%
  intersect(names(main_df))

env_vars <- unique(c(climate_vars, soil_vars))

if (length(env_vars) == 0) stop("No climate/soil variables found")
if (length(climate_vars) == 0) warning("No climate variables found")
if (length(soil_vars) == 0) warning("No soil variables found")

message("Climate vars: ", paste(climate_vars, collapse = ", "))
message("Soil vars: ", paste(soil_vars, collapse = ", "))


# ====================================================
### X. Ploidy diagnostic and optional normalization ----
# ====================================================
## add ploidy
main_df$ploidy <- ifelse(main_df$Population_name %in% tet, 4, 2)

rda_ploidy <- rda(main_df[, og_cols] ~ ploidy, data = main_df)
anova(rda_ploidy, permutations = 999)
#====================================================

cnv_wide <- main_df %>%
  select(Population_name, all_of(og_cols)) %>%
  distinct() %>%
  as.data.frame()

rownames(cnv_wide) <- cnv_wide$Population_name
cnv_wide$Population_name <- NULL

# keep only variable families
keep_fams <- apply(cnv_wide, 2, function(x) sd(x, na.rm = TRUE) > 0)
cnv_raw <- cnv_wide[, keep_fams, drop = FALSE]

# z-score within family
cnv_z <- scale(cnv_raw, center = TRUE, scale = TRUE)
cnv_z <- as.data.frame(cnv_z)

# remove any families that became all-NA after scaling
cnv_z <- cnv_z[, colSums(is.na(cnv_z)) == 0, drop = FALSE]

# align raw to z columns
shared_fams <- intersect(colnames(cnv_raw), colnames(cnv_z))
cnv_raw <- cnv_raw[, shared_fams, drop = FALSE]
cnv_z   <- cnv_z[, shared_fams, drop = FALSE]

# ====================================================
### 3. Build predictor matrices ------------------------
# ====================================================

meta_df <- main_df %>%
  select(all_of(c("Population_name", id_cols, structure_cols, "ploidy", climate_vars, soil_vars))) %>%
  distinct() %>%
  as.data.frame()

rownames(meta_df) <- meta_df$Population_name
meta_df$Population_name <- NULL

# scale predictors for RDA
scale_if_present <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (length(vars) > 0) {
    df[, vars] <- lapply(df[, vars, drop = FALSE], function(x) as.numeric(scale(x)))
  }
  df
}

meta_df <- scale_if_present(meta_df, c(structure_cols, climate_vars, soil_vars))
# ploidy intentionally left unscaled

# align rows
common_ids <- Reduce(intersect, list(
  rownames(cnv_raw),
  rownames(cnv_z),
  rownames(meta_df)
))

cnv_raw <- cnv_raw[common_ids, , drop = FALSE]
cnv_z   <- cnv_z[common_ids, , drop = FALSE]
meta_df <- meta_df[common_ids, , drop = FALSE]

# ====================================================
### 4. Helper functions ------------------------
# ====================================================
source("scripts/rda_helpers.R")

# ====================================================
### 5. Partial RDA analyses ------------------------
# ====================================================

# A. combined environment after conditioning on population structure
res_raw_env <- run_partial_rda(
  Y = cnv_raw,
  meta_df = meta_df,
  explanatory = env_vars,
  conditional = structure_cols,
  label = "Raw CNV ~ climate + soil | structure + ploidy"
)

res_z_env <- run_partial_rda(
  Y = cnv_z,
  meta_df = meta_df,
  explanatory = env_vars,
  conditional = c(structure_cols, "ploidy"),
  label = "Z-score CNV ~ climate + soil | structure + ploidy"
)

# B. pure climate after conditioning on soil + structure
res_raw_clim <- if (length(climate_vars) > 0) {
  run_partial_rda(
    Y = cnv_raw,
    meta_df = meta_df,
    explanatory = climate_vars,
    conditional = unique(c(structure_cols, soil_vars)),
    label = "Raw CNV ~ climate | soil + structure + ploidy"
  )
} else NULL

res_z_clim <- if (length(climate_vars) > 0) {
  run_partial_rda(
    Y = cnv_z,
    meta_df = meta_df,
    explanatory = climate_vars,
    conditional = unique(c(structure_cols, soil_vars)),
    label = "Z-score CNV ~ climate | soil + structure + ploidy"
  )
} else NULL

# C. pure soil after conditioning on climate + structure
res_raw_soil <- if (length(soil_vars) > 0) {
  run_partial_rda(
    Y = cnv_raw,
    meta_df = meta_df,
    explanatory = soil_vars,
    conditional = unique(c(structure_cols, climate_vars)),
    label = "Raw CNV ~ soil | climate + structure"
  )
} else NULL

res_z_soil <- if (length(soil_vars) > 0) {
  run_partial_rda(
    Y = cnv_z,
    meta_df = meta_df,
    explanatory = soil_vars,
    conditional = unique(c(structure_cols, climate_vars)),
    label = "Z-score CNV ~ soil | climate + structure"
  )
} else NULL

# save core results
save_rda_outputs <- function(res_obj, prefix) {
  if (is.null(res_obj)) return(NULL)
  write.csv(res_obj$overall, file.path(out_dir, paste0(prefix, "_overall.csv")), row.names = TRUE)
  write.csv(res_obj$terms,   file.path(out_dir, paste0(prefix, "_terms.csv")), row.names = TRUE)
  write.csv(res_obj$r2,      file.path(out_dir, paste0(prefix, "_r2.csv")), row.names = FALSE)
}

save_rda_outputs(res_raw_env,  "rda_raw_env")
save_rda_outputs(res_z_env,    "rda_z_env")
save_rda_outputs(res_raw_clim, "rda_raw_climate")
save_rda_outputs(res_z_clim,   "rda_z_climate")
save_rda_outputs(res_raw_soil, "rda_raw_soil")
save_rda_outputs(res_z_soil,   "rda_z_soil")

# ====================================================
### 6. Variance partitioning ------------------------
# ====================================================

# among climate, soil, structure
vp_raw <- varpart(cnv_raw,
                  meta_df[, climate_vars, drop = FALSE],
                  meta_df[, soil_vars, drop = FALSE],
                  meta_df[, structure_cols, drop = FALSE])

vp_z <- varpart(cnv_z,
                meta_df[, climate_vars, drop = FALSE],
                meta_df[, soil_vars, drop = FALSE],
                meta_df[, structure_cols, drop = FALSE])

vp_raw_df <- extract_varpart_summary(vp_raw, "Raw CNV")
vp_z_df   <- extract_varpart_summary(vp_z, "Z-score CNV")

write.csv(vp_raw_df, file.path(out_dir, "variance_partition_raw.csv"), row.names = FALSE)
write.csv(vp_z_df,   file.path(out_dir, "variance_partition_zscore.csv"), row.names = FALSE)

# explicit model-level summary table
model_summary <- bind_rows(
  res_raw_env$r2 %>% mutate(model = "Raw: climate+soil | structure"),
  res_z_env$r2   %>% mutate(model = "Z-score: climate+soil | structure"),
  if (!is.null(res_raw_clim)) res_raw_clim$r2 %>% mutate(model = "Raw: climate | soil+structure"),
  if (!is.null(res_z_clim))   res_z_clim$r2   %>% mutate(model = "Z-score: climate | soil+structure"),
  if (!is.null(res_raw_soil)) res_raw_soil$r2 %>% mutate(model = "Raw: soil | climate+structure"),
  if (!is.null(res_z_soil))   res_z_soil$r2   %>% mutate(model = "Z-score: soil | climate+structure")
)

write.csv(model_summary, file.path(out_dir, "rda_model_summary.csv"), row.names = FALSE)

# ====================================================
### 7. Figures ------------------------
# ====================================================

# 7.1 compare raw vs z-score for key partial RDA models
p_model_comp <- model_summary %>%
  pivot_longer(c(r2, adj_r2), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = fct_reorder(model, value), y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c("r2" = "#4daf4a", "adj_r2" = "#377eb8"),
                    labels = c("Adjusted R²","R²")) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Explained variance",
       fill = NULL,
       title = "Comparison of partial RDA models")

ggsave(file.path(out_dir, "Fig_1_partial_rda_model_comparison.pdf"),
       p_model_comp, width = 8.2, height = 5.6)

# 7.2 term-wise environmental contributions for combined-environment models
prep_term_plot <- function(tt, label) {
  tt <- as.data.frame(tt)
  tt$term <- rownames(tt)
  rownames(tt) <- NULL
  
  pcol <- grep("Pr", names(tt), value = TRUE)[1]
  vcol <- if ("Variance" %in% names(tt)) {
    "Variance"
  } else if ("ChiSquare" %in% names(tt)) {
    "ChiSquare"
  } else {
    names(tt)[2]
  }
  
  tt %>%
    filter(term != "Residual") %>%
    mutate(
      model = label,
      significant = .data[[pcol]] < 0.05,
      variance = .data[[vcol]]
    ) %>%
    select(term, model, significant, variance)
}

term_df <- bind_rows(
  prep_term_plot(res_raw_env$terms, "Raw CNV"),
  prep_term_plot(res_z_env$terms, "Z-score CNV")
)

write.csv(term_df, file.path(out_dir, "rda_env_term_plot_table.csv"), row.names = FALSE)

all_terms <- unique(term_df$term)
bio_terms <- sort(all_terms[str_detect(all_terms, "^bio")])
other_terms <- sort(all_terms[!str_detect(all_terms, "^bio")])

term_order <- c(bio_terms, other_terms)
term_df <- term_df %>%
  mutate(term = factor(term, levels = term_order))

p_terms <- ggplot(term_df,
                  aes(x = term, y = variance, fill = significant)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~model, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = "#1b9e77", "FALSE" = "gray70")) +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Constrained variance",
       fill = "p < 0.05",
       title = "Environmental terms contributing to CNV structure")

ggsave(file.path(out_dir, "Fig_2_partial_rda_env_term_variance.pdf"),
       p_terms, width = 7, height = 3)


# 7.3 ordinations
al_meta<-read.csv("data/231/AL/genDist/AL_gen_div_covar_svd.csv")
#al_meta$Sample[19]<-"73_3a"

pop_col_use <- if ("pop_name" %in% names(res_raw_env$sites)) {
  "pop_name"
} else if ("Population_name" %in% names(res_raw_env$sites)) {
  "Population_name"
} else NULL

p_ord_raw <- plot_rda_ordination(res_raw_clim, pop_col = pop_col_use, max_arrows = 8,clusters = al_meta)
p_ord_z   <- plot_rda_ordination(res_z_clim,   pop_col = pop_col_use, max_arrows = 8,clusters = al_meta)

ggsave(file.path(out_dir, "Fig_3_partial_rda_ordination_raw.pdf"),
       p_ord_raw, width = 7.5, height = 6)
ggsave(file.path(out_dir, "Fig_4_partial_rda_ordination_zscore.pdf"),
       p_ord_z, width = 7.5, height = 4)
ggsave(file.path(out_dir, "Fig_5_partial_rda_ordination_comparison.pdf"),
       p_ord_raw + p_ord_z, width = 14, height = 6)

# 7.4 variance partitioning diagrams (base plotting)
pdf(file.path(out_dir, "Fig_6_variance_partitioning_raw.pdf"), width = 7, height = 7)
plot(vp_raw, bg = c("#8dd3c7", "#bebada", "#fb8072"),
     Xnames = c("Climate", "Soil", "Structure"),
     id.size = 1)
title("Variance partitioning: Raw CNV")
dev.off()

pdf(file.path(out_dir, "Fig_7_variance_partitioning_zscore.pdf"), width = 5, height = 5)
plot(vp_z, bg = c("#8dd3c7", "#bebada", "#fb8072"),
     Xnames = c("Climate", "Soil", "Structure"),
     id.size = 1)
title("Variance partitioning: Z-score CNV")
dev.off()

# ====================================================
### 8. Concise topline summary ------------------------
# ====================================================

analysis_summary <- tibble(
  analysis = c(
    "Raw: climate+soil | structure",
    "Z-score: climate+soil | structure",
    "Raw: climate | soil+structure",
    "Z-score: climate | soil+structure",
    "Raw: soil | climate+structure",
    "Z-score: soil | climate+structure"
  ),
  adj_r2 = c(
    res_raw_env$r2$adj_r2,
    res_z_env$r2$adj_r2,
    if (!is.null(res_raw_clim)) res_raw_clim$r2$adj_r2 else NA_real_,
    if (!is.null(res_z_clim))   res_z_clim$r2$adj_r2 else NA_real_,
    if (!is.null(res_raw_soil)) res_raw_soil$r2$adj_r2 else NA_real_,
    if (!is.null(res_z_soil))   res_z_soil$r2$adj_r2 else NA_real_
  ),
  r2 = c(
    res_raw_env$r2$r2,
    res_z_env$r2$r2,
    if (!is.null(res_raw_clim)) res_raw_clim$r2$r2 else NA_real_,
    if (!is.null(res_z_clim))   res_z_clim$r2$r2 else NA_real_,
    if (!is.null(res_raw_soil)) res_raw_soil$r2$r2 else NA_real_,
    if (!is.null(res_z_soil))   res_z_soil$r2$r2 else NA_real_
  )
)

write.csv(analysis_summary,
          file.path(out_dir, "analysis_summary_topline.csv"),
          row.names = FALSE)



# ====================================================
## 5. Family-wise association testing ------------------------
#    raw counts only, NB models
#    climate + soil predictors only
#    controlling only for population structure
# ====================================================
################## Z-SCORE CNV–ENV ANALYSIS ##################

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
})

# =========================================================
### 0. INPUT ------------------------
# =========================================================
in_file <- "outputs/AL/eco/main_df_final_CNV_ENV_ready_imputed.csv"
out_dir <- "outputs/AL/eco/cnv_env_univariate"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

main_df <- data.table::fread(in_file)
# tetraploids
tet<-c("WS1","al1","AL27","PU6")
## add ploidy
main_df$ploidy <- ifelse(main_df$Population_name %in% tet, 4, 2)

# =========================================================
### 1. DEFINE COLUMNS ------------------------
# =========================================================
structure_cols <- intersect(c("PC1", "PC2"), names(main_df))
if (length(structure_cols) < 2) stop("PC1 and/or PC2 missing")

og_cols <- grep("^OG", names(main_df), value = TRUE)
if (length(og_cols) == 0) stop("No OG columns found")

bio_cols  <- grep("^bio\\d+$", names(main_df), value = TRUE)
soil_cols <- grep("cm", names(main_df), value = TRUE)

# leave out TE
env_cols <- unique(c(bio_cols, soil_cols))
env_cols <- intersect(env_cols, names(main_df))
if (length(env_cols) == 0) stop("No climate/soil variables found")

# =========================================================
### 2. BUILD LONG TABLE + FAMILY-WISE Z-SCORE ------------------------
# =========================================================
cnv_long_z <- main_df %>%
  select(all_of(c("Population_name", "ploidy", structure_cols, env_cols, og_cols))) %>%
  pivot_longer(
    cols = all_of(og_cols),
    names_to = "family_id",
    values_to = "copy_number"
  ) %>%
  mutate(
    copy_number = as.numeric(copy_number),
    ploidy = factor(ploidy)
  )

# keep only variable families
family_stats <- cnv_long_z %>%
  group_by(family_id) %>%
  summarise(
    mean_copy = mean(copy_number, na.rm = TRUE),
    sd_copy = sd(copy_number, na.rm = TRUE),
    n_states = n_distinct(copy_number),
    .groups = "drop"
  )

families_keep <- family_stats %>%
  filter(sd_copy > 0, n_states > 1) %>%
  pull(family_id)

cnv_long_z <- cnv_long_z %>%
  filter(family_id %in% families_keep) %>%
  left_join(family_stats, by = "family_id") %>%
  mutate(
    cnv_z = ifelse(sd_copy > 0, (copy_number - mean_copy) / sd_copy, 0)
  ) %>%
  select(-mean_copy, -sd_copy, -n_states)

# standardize predictors for comparable coefficients
cnv_long_z <- cnv_long_z %>%
  mutate(across(all_of(c(structure_cols, env_cols)), ~ as.numeric(scale(.x))))

write.csv(
  family_stats,
  file.path(out_dir, "family_zscore_baseline_statistics.csv"),
  row.names = FALSE
)

# =========================================================
### 3. HELPER FUNCTIONS ------------------------
# =========================================================
run_family_lm_z <- function(df, env_var, covars = c("PC1", "PC2", "ploidy")) {
  rhs <- paste(c(env_var, covars), collapse = " + ")
  form <- as.formula(paste0("cnv_z ~ ", rhs))
  
  fit <- lm(form, data = df)
  td <- broom::tidy(fit)
  gl <- broom::glance(fit)
  
  tr <- td %>% filter(term == env_var)
  if (nrow(tr) == 0) return(tibble())
  
  tibble(
    estimate = tr$estimate,
    std.error = tr$std.error,
    statistic = tr$statistic,
    p.value = tr$p.value,
    r.squared = gl$r.squared,
    adj.r.squared = gl$adj.r.squared
  )
}
# =========================================================
### 4. FAMILY-WISE SCREENING WITH Z-SCORED CNV ------------------------
# =========================================================
family_env_z_results <- purrr::map_dfr(env_cols, function(ev) {
  cnv_long_z %>%
    filter(!is.na(.data[[ev]]), !is.na(cnv_z), !is.na(PC1), !is.na(PC2), !is.na(ploidy)) %>%
    group_by(family_id) %>%
    group_modify(~ {
      out <- tryCatch(
        run_family_lm_z(.x, env_var = ev, covars = c(structure_cols, "ploidy")),
        error = function(e) tibble()
      )
      out
    }) %>%
    ungroup() %>%
    mutate(env_var = ev)
}) %>%
  group_by(env_var) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup()

## HERE ------------------
ogList<-data.frame(fread("data/ogList_with_functions_v2.csv"))
function_short<-data.frame(fread("data/og_fam_function_short.tsv",h=T))
function_short<-function_short[match(ogList$Family.ID,function_short$Family.ID),]
ogList$func_short<-function_short$short_title

og_matched<-ogList[match(family_env_z_results$family_id,ogList$Family.ID),c("func_short","FunctionalCategory")]
family_env_z_results<-data.frame(short_func=og_matched[,1],func_cat=og_matched[,2],family_env_z_results)

write.csv(
  family_env_z_results,
  file.path(out_dir, "familywise_zscore_environment_associations.csv"),
  row.names = FALSE
)

family_env_z_summary <- family_env_z_results %>%
  group_by(env_var) %>%
  summarise(
    n_tested = sum(!is.na(p.value)),
    n_sig_p05 = sum(p.value < 0.05, na.rm = TRUE),
    n_sig_fdr05 = sum(q.value < 0.05, na.rm = TRUE),
    median_abs_effect = median(abs(estimate), na.rm = TRUE),
    median_adj_r2 = median(adj.r.squared, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_fdr05), desc(n_sig_p05))

write.csv(
  family_env_z_summary,
  file.path(out_dir, "familywise_zscore_environment_summary.csv"),
  row.names = FALSE
)

# =========================================================
### 5. FIGURE 1: SUMMARY BARPLOT ------------------------
# =========================================================

all_env_vars <- unique(family_env_z_summary$env_var)

bio_vars <- sort(all_env_vars[str_detect(all_env_vars, "^bio")])
other_vars <- sort(all_env_vars[!str_detect(all_env_vars, "^bio")])
var_order <- c(bio_vars, other_vars)
family_env_z_summary_ordered <- family_env_z_summary %>%
  mutate(env_var = factor(env_var, levels = var_order))

p_env_bar_z <- family_env_z_summary_ordered %>%
  pivot_longer(
    cols = c(n_sig_p05, n_sig_fdr05),
    names_to = "metric",
    values_to = "n_families"
  ) %>%
  ggplot(aes(x = env_var, y = n_families, fill = metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  coord_flip() +
  scale_fill_manual(
    values = c("n_sig_p05" = "#984ea3", "n_sig_fdr05" = "#4daf4a"),
    labels = c("n_sig_p05" = "p < 0.05", "n_sig_fdr05" = "FDR < 0.05")
  ) +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "Number of families",
    fill = NULL,
    title = "Family-wise CNV–environment associations using z-scored CNV"
  )

p_env_bar_z

ggsave(
  file.path(out_dir, "Fig_1_zscore_familywise_hits.pdf"),
  p_env_bar_z, width = 4, height = 6
)

# =========================================================
### 6. FIGURE 2: HEATMAP OF SIGNIFICANT / STRONG FAMILIES ------------------------
# =========================================================
# use FDR hits if present, otherwise fall back to nominal hits

library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)

heat_df_z <- family_env_z_results %>%
  filter(!is.na(p.value)) %>%
  mutate(
    sig_class = case_when(
      q.value < 0.05 ~ "FDR",
      p.value < 0.05 ~ "p05",
      TRUE ~ "NS"
    )
  )

fam_rank <- heat_df_z %>%
  filter(p.value < 0.05) %>%
  group_by(family_id) %>%
  summarise(
    best_p = min(p.value, na.rm = TRUE),
    best_q = min(q.value, na.rm = TRUE),
    max_abs_effect = max(abs(estimate), na.rm = TRUE),
    .groups = "drop"
  )

heat_subtitle <- "Stars indicate FDR < 0.05; families shown have nominal p < 0.05"

top_fams_heat <- fam_rank %>%
  arrange(best_p, best_q, desc(max_abs_effect)) %>%
  slice(1:min(50, n())) %>%
  pull(family_id)

func_levels <- sort(unique(heat_df_z$func_cat))
func_cols <- setNames(hue_pal(l = 65, c = 100)(length(func_levels)), func_levels)

fam_meta <- heat_df_z %>%
  distinct(family_id, short_func, func_cat) %>%
  filter(family_id %in% top_fams_heat) %>%
  left_join(fam_rank, by = "family_id") %>%
  arrange(func_cat, best_p, best_q, desc(max_abs_effect), short_func, family_id) %>%
  mutate(
    family_id = factor(family_id, levels = rev(family_id)),
    y_lab = paste0("<span style='color:", func_cols[func_cat], ";'>", short_func, "</span>")
  )

heat_df_plot <- heat_df_z %>%
  filter(family_id %in% top_fams_heat) %>%
  left_join(fam_meta %>% select(family_id, short_func, func_cat, y_lab), by = "family_id") %>%
  mutate(
    env_var = factor(env_var, levels = unique(heat_df_z$env_var)),
    family_id = factor(family_id, levels = levels(fam_meta$family_id))
  )

n_env <- nlevels(heat_df_plot$env_var)

right_lab_df <- fam_meta %>%
  mutate(xpos = n_env + 0.35)

legend_df <- data.frame(
  func_cat = factor(func_levels, levels = func_levels),
  x = Inf,
  y = Inf
)

p_heat_z <- ggplot(heat_df_plot, aes(x = as.numeric(env_var), y = family_id, fill = estimate)) +
  geom_tile(color = "white") +
  geom_point(
    data = subset(heat_df_plot, sig_class == "FDR"),
    shape = 8, size = 1.7, color = "black"
  ) +
  geom_text(
    data = right_lab_df,
    aes(x = xpos, y = family_id, label = as.character(family_id)),
    inherit.aes = FALSE, hjust = 0, size = 2.8
  ) +
  geom_point(
    data = legend_df,
    aes(x = x, y = y, color = func_cat),
    inherit.aes = FALSE, alpha = 0, size = 3, show.legend = TRUE
  ) +
  scale_x_continuous(
    breaks = seq_len(n_env),
    labels = levels(heat_df_plot$env_var),
    expand = expansion(mult = c(0, 0), add = c(0, 1.8))
  ) +
  scale_y_discrete(
    labels = setNames(as.character(fam_meta$y_lab), as.character(fam_meta$family_id))
  ) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    name = "Effect (Z-score)"
  ) +
  scale_color_manual(values = func_cols, name = "Functional category") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = ggtext::element_markdown(),
    plot.margin = margin(5.5, 140, 5.5, 5.5),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    color = guide_legend(order = 2, override.aes = list(alpha = 1, size = 4))
  ) +
  labs(
    x = NULL,
    y = "Gene family",
    title = "Family × environment associations using z-scored CNV",
    subtitle = heat_subtitle
  )

ggsave(
  file.path(out_dir, "Fig_2_zscore_family_environment_heatmap_v2.1.pdf"),
  p_heat_z, width = 8, height = 10
)

topfams<-heat_df_plot[heat_df_plot$sig_class!="NS",]
topfams<-topfams[!duplicated(topfams$family_id),]
topfams<-topfams[order(topfams$func_cat.y),]

write.csv(
  topfams,
  file.path(out_dir,"AL_top_families_summary.csv"),
  row.names = F
)


# =========================================================
### 7. FIGURE 3: REACTION-NORM STYLE PLOTS ------------------------
# =========================================================
# choose strongest env variables
top_envs <- family_env_z_summary %>%
  slice(1:min(3, n())) %>%
  pull(env_var)

# choose strongest family × env combinations
top_hits_z <- family_env_z_results %>%
  filter(env_var %in% top_envs) %>%
  arrange(q.value, p.value) %>%
  slice(1:min(12, n())) %>%
  select(family_id, env_var) %>%
  distinct()


fam_lookup <- family_env_z_results %>%
  distinct(family_id, short_func, func_cat)

reaction_plot_list <- purrr::pmap(
  list(top_hits_z$family_id, top_hits_z$env_var),
  function(fid, ev) {
    
    meta <- fam_lookup %>% filter(family_id == fid) %>% slice(1)
    
    dat <- cnv_long_z %>%
      filter(family_id == fid) %>%
      filter(!is.na(.data[[ev]]), !is.na(cnv_z), !is.na(PC1), !is.na(PC2))
    
    ggplot(dat, aes(x = .data[[ev]], y = cnv_z)) +
      geom_point(size = 2, alpha = 0.75) +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      theme_classic(base_size = 10) +
      labs(
        title = paste0(meta$short_func, " (", fid, ")"),
        subtitle = meta$func_cat,
        x = ev,
        y = "CNV (z-score)"
      )
  }
)

if (length(reaction_plot_list) > 0) {
  pdf(file.path(out_dir, "Fig_3_zscore_reaction_norms_top_hits.pdf"), width = 10, height = 12)
  print(patchwork::wrap_plots(reaction_plot_list, ncol = 2))
  dev.off()
}

# =========================================================
### TOPLINE COMPARISON TABLE ------------------------
# =========================================================
analysis_summary_z <- family_env_z_summary %>%
  transmute(
    env_var,
    n_tested,
    n_sig_p05,
    n_sig_fdr05,
    median_abs_effect,
    median_adj_r2
  )

write.csv(
  analysis_summary_z,
  file.path(out_dir, "analysis_summary_zscore_familywise.csv"),
  row.names = FALSE
)


