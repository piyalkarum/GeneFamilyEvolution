###################### TE ANALYSIS ######################
#-------------------------------------------

# THALIANA ---------------------------------------

library(dplyr)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(stringr)
library(purrr)
library(forcats)

# =========================
# SETTINGS
# =========================
out_dir <- "outputs/AT/TE_assessment/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

te_classes <- c(
  "DNA/Helitron","DNA/MULE-MuDR","LTR/Copia", "LTR/Gypsy",
  "SINE/tRNA", "SINE/unknown","LINE/L1", "LINE/unknown",
  "LTR/Ty3", "LTR/unknown"
)

windows <- c(0, 1000, 10000, 25000, 50000)

fam_file <- "data/231/AT/cnv_within/all_231_blast_duprm_v1_joined.tsv"
nonfam_file <- "data/AT_non_family_genes_blast_groups.csv"
rm_dir <- "data/TE_profiling/RepeatMasker_results/AT/v2"

# ====================================
# helper functions are in TE_helper.R
# ====================================
source("scripts/TE_helpers.R")

# =========================
# INPUT GENE TABLES
# =========================

fam_genes <- fread(fam_file) %>%
  as.data.frame() %>%
  transmute(
    assembly,
    gene = gene,
    scaffold,
    start = as.integer(start),
    end = as.integer(end),
    group = "Family"
  )

nonfam_genes <- read.csv(nonfam_file) %>%
  filter(copy < 2) %>%
  transmute(
    assembly,
    gene,
    scaffold,
    start = as.integer(start),
    end = as.integer(end),
    group = "Single-copy"
  )

all_genes <- bind_rows(fam_genes, nonfam_genes)

# =========================
# RUN PER ASSEMBLY
# =========================

fls <- list.files(rm_dir, full.names = TRUE, pattern = "\\.fna\\.out$|\\.out$")
all_metrics <- list()
all_classes <- list()

pb <- txtProgressBar(max = length(fls), style = 3, width = 50)

for (i in seq_along(fls)) {
  setTxtProgressBar(pb, i)
  
  te_df <- read_repeatmasker(fls[i])
  assnam <- gsub("-", "_", gsub("\\.fna\\.out$|\\.out$", "", basename(fls[i])))
  
  genes_ass <- all_genes %>% filter(assembly == assnam)
  if (nrow(genes_ass) == 0) next
  
  for (w in windows) {
    met <- summarize_te_metrics(genes_ass, te_df, buffer = w)
    cls <- summarize_te_classes(genes_ass, te_df, buffer = w)
    
    all_metrics[[paste0(assnam, "_", w)]] <- met
    if (!is.null(cls)) all_classes[[paste0(assnam, "_", w)]] <- cls
  }
}
close(pb)

metrics_df <- bind_rows(all_metrics)
class_df   <- bind_rows(all_classes)

write.csv(metrics_df, file.path(out_dir, "TE_metrics_nested_windows.csv"), row.names = FALSE)
write.csv(class_df,   file.path(out_dir, "TE_class_counts_nested_windows.csv"), row.names = FALSE)

# =========================
# SUMMARY TABLES
# =========================
summary_df <- metrics_df %>%
  group_by(group, window) %>%
  summarise(
    mean_te_count = mean(n_te, na.rm = TRUE),
    mean_te_bp = mean(te_bp, na.rm = TRUE),
    mean_te_fraction = mean(te_fraction, na.rm = TRUE),
    mean_nearest_te = mean(nearest_te_dist, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, file.path(out_dir, "TE_metrics_summary_by_group_window.csv"), row.names = FALSE)

# =========================
# STATISTICS
# =========================

stats_df <- metrics_df %>%
  group_by(window) %>%
  summarise(
    p_te_count = wilcox.test(n_te ~ group)$p.value,
    p_te_fraction = wilcox.test(te_fraction ~ group)$p.value,
    p_nearest_te = wilcox.test(nearest_te_dist ~ group)$p.value,
    .groups = "drop"
  )

write.csv(stats_df, file.path(out_dir, "TE_metrics_stats_by_window.csv"), row.names = FALSE)

# =========================
# VISUALIZATIONS
# =========================
metrics_df<-read.csv(file.path(out_dir, "TE_metrics_nested_windows.csv"))
class_df <-read.csv(file.path(out_dir, "TE_class_counts_nested_windows.csv"))

# 1. TE landscape profile across windows
p_profile <- summary_df %>%
  mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb"))) %>%
  ggplot(aes(x = window, y = mean_te_fraction, group = group, color = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Mean TE fraction",
    title = "TE landscape around genes across nested windows",
    subtitle = "Comparison of duplicated-family genes and single-copy genes"
  )

ggsave(file.path(out_dir, "TE_profile_fraction_by_window.pdf"), p_profile, width = 5, height = 3)

# 2. Distribution of TE fraction per locus
library(dplyr)
library(rstatix)

# Calculate significance for each window
signif_data <- metrics_df %>%
  filter(group %in% c("Family", "Single-copy")) %>%
  group_by(window) %>%
  summarise(
    max_y = max(te_fraction, na.rm = TRUE),
    y_pos = log1p(max_y * 1.1)  # 10% above max value, then transformed
  ) %>%
  left_join(
    metrics_df %>%
      filter(group %in% c("Family", "Single-copy")) %>%
      group_by(window) %>%
      wilcox_test(te_fraction ~ group) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance() %>%
      mutate(asterisk = case_when(
        p.adj.signif == "****" ~ "****",
        p.adj.signif == "***" ~ "***",
        p.adj.signif == "**" ~ "**",
        p.adj.signif == "*" ~ "*",
        TRUE ~ ""
      )),
    by = "window"
  ) %>% mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb")))

p_violin <- metrics_df %>%
  mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb"))) %>%
  ggplot(aes(x = group, y = te_fraction, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  facet_wrap(~window, scales = "free_y", nrow = 1) +
  scale_y_continuous(
    trans = "log1p",
    limits = c(0, 3)
  ) +
  # Add significance asterisks
  geom_text(
    data = signif_data,
    aes(x = 1.5, y = y_pos, label = asterisk),  
    size = 6,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "TE fraction per locus (log1p transformed)",
    title = "Distribution of TE occupancy around loci"
  ) + 
  theme(legend.position = "none")

ggsave(file.path(out_dir, "TE_fraction_violin_nested_windows_log.pdf"), p_violin, width = 10, height = 3)

##################### add a truncated y-axis #####################
library(dplyr)
library(rstatix)
library(ggplot2)

window_levels <- c("gene_body", "pm1kb", "pm10kb", "pm25kb", "pm50kb")

metrics_plot <- metrics_df %>%
  filter(group %in% c("Family", "Single-copy")) %>%
  mutate(
    window = factor(window, levels = window_levels),
    te_fraction_log = log1p(te_fraction)
  )

# split/compress y-axis from Q3 to 2
q3_y <- quantile(metrics_plot$te_fraction_log, 0.95, na.rm = TRUE)

break_low  <- q3_y
break_high <- 2.8
gap_size   <- 0.10

compress_y <- function(y) {
  ifelse(
    y <= break_low,
    y,
    ifelse(
      y < break_high,
      break_low + ((y - break_low) / (break_high - break_low)) * gap_size,
      y - (break_high - break_low) + gap_size
    )
  )
}

metrics_plot <- metrics_plot %>%
  mutate(te_fraction_plot = compress_y(te_fraction_log))

# test significance per window
signif_data <- metrics_plot %>%
  group_by(window) %>%
  wilcox_test(te_fraction ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(
    asterisk = case_when(
      p.adj.signif == "****" ~ "****",
      p.adj.signif == "***"  ~ "***",
      p.adj.signif == "**"   ~ "**",
      p.adj.signif == "*"    ~ "*",
      TRUE ~ ""
    ),
    y_pos = compress_y(2.85)
  )

p_violin <- metrics_plot %>%
  ggplot(aes(x = group, y = te_fraction_plot, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  facet_wrap(~window, scales = "fixed", nrow = 1) +
  scale_y_continuous(
    limits = c(0, compress_y(3)),
    breaks = c(0, break_low, compress_y(2), compress_y(3)),
    labels = c("0", paste0("Q3=", round(break_low, 2)), "2", "3")
  ) +
  geom_hline(yintercept = break_low, linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = break_low + gap_size, linetype = "dashed", linewidth = 0.25) +
  geom_text(
    data = signif_data,
    aes(x = 1.5, y = y_pos, label = asterisk),
    size = 6,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "TE fraction per locus (log1p transformed; compressed above Q3)",
    title = "Distribution of TE occupancy around loci"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "TE_fraction_violin_nested_windows_log_q3_compressed.pdf"),
  p_violin,
  width = 10,
  height = 3
)

p_violin
##################################################################

# 3. TE class enrichment heatmap
class_summary <- class_df %>%
  group_by(group, window, repeat_class) %>%
  summarise(mean_count = mean(n_te_class, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_count, values_fill = 0) %>%
  dplyr::rename(Single_copy = `Single-copy`) %>%
  mutate(effect = log2((Family + 0.01) / (Single_copy + 0.01)))

p_heat <- class_summary %>%
  mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb")),
         repeat_class = fct_reorder(repeat_class, effect, .fun = max)) %>%
  ggplot(aes(x = window, y = repeat_class, fill = effect)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "TE class",
    fill = "log2(Fam./Sing)",
    title = "TE class enrichment around duplicated-family genes"
  )  + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "TE_class_enrichment_heatmap.pdf"), p_heat, width = 5, height = 3)

# 4. Nearest TE distance
p_dist <- metrics_df %>%
  mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb"))) %>%
  ggplot(aes(x = window, y = nearest_te_dist, color = group, group = group)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.2, position = position_dodge(width = 0.1)) +
  stat_summary(fun = median, geom = "point", size = 3, position = position_dodge(width = 0.1)) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Median distance to nearest TE (bp)",
    title = "Distance to nearest TE across nested windows"
  ) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(file.path(out_dir, "TE_nearest_distance_profile.pdf"), p_dist, width = 5, height = 3)



# LYRATA -----------------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(stringr)
library(purrr)
library(forcats)

# =========================
# SETTINGS
# =========================
out_dir <- "outputs/AL/TE_assessment/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

te_classes <- c(
  "DNA/Helitron","DNA/MULE-MuDR","LTR/Copia", "LTR/Gypsy",
  "SINE/tRNA", "SINE/unknown","LINE/L1", "LINE/unknown",
  "LTR/Ty3", "LTR/unknown"
)

windows <- c(0, 1000, 10000, 25000, 50000)

fam_file <- "data/231/all_231_blast_duprm_v1_joined.tsv"
nonfam_file <- "data/AL_non_family_genes_blast_groups.csv"
rm_dir <- "data/TE_profiling/RepeatMasker_results/AL"

# ====================================
# helper functions are in TE_helper.R
# ====================================
source("scripts/TE_helpers.R")

# =========================
# INPUT GENE TABLES
# =========================

fam_genes <- fread(fam_file) %>%
  as.data.frame() %>%
  transmute(
    assembly,
    gene = gene,
    scaffold,
    start = as.integer(start),
    end = as.integer(end),
    group = "Family"
  )

nonfam_genes <- read.csv(nonfam_file) %>%
  filter(copy < 2) %>%
  transmute(
    assembly,
    gene,
    scaffold,
    start = as.integer(start),
    end = as.integer(end),
    group = "Single-copy"
  )

all_genes <- bind_rows(fam_genes, nonfam_genes)

# =========================
# RUN PER ASSEMBLY
# =========================

fls <- list.files(rm_dir, full.names = TRUE, pattern = "\\.fa\\.out$|\\.out$")
all_metrics <- list()
all_classes <- list()

pb <- txtProgressBar(max = length(fls), style = 3, width = 50)

for (i in seq_along(fls)) {
  setTxtProgressBar(pb, i)
  
  te_df <- read_repeatmasker(fls[i])
  assnam <- gsub("-", "_", gsub("\\.fa\\.out$|\\.out$", "", basename(fls[i])))
  
  genes_ass <- all_genes %>% filter(assembly == assnam)
  if (nrow(genes_ass) == 0) next
  
  for (w in windows) {
    met <- summarize_te_metrics(genes_ass, te_df, buffer = w)
    cls <- summarize_te_classes(genes_ass, te_df, buffer = w)
    
    all_metrics[[paste0(assnam, "_", w)]] <- met
    if (!is.null(cls)) all_classes[[paste0(assnam, "_", w)]] <- cls
  }
}
close(pb)

metrics_df <- bind_rows(all_metrics)
class_df   <- bind_rows(all_classes)

write.csv(metrics_df, file.path(out_dir, "TE_metrics_nested_windows.csv"), row.names = FALSE)
write.csv(class_df,   file.path(out_dir, "TE_class_counts_nested_windows.csv"), row.names = FALSE)

# =========================
# SUMMARY TABLES
# =========================

summary_df <- metrics_df %>%
  group_by(group, window) %>%
  summarise(
    mean_te_count = mean(n_te, na.rm = TRUE),
    mean_te_bp = mean(te_bp, na.rm = TRUE),
    mean_te_fraction = mean(te_fraction, na.rm = TRUE),
    mean_nearest_te = mean(nearest_te_dist, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, file.path(out_dir, "TE_metrics_summary_by_group_window.csv"), row.names = FALSE)

# =========================
# STATISTICS
# =========================

stats_df <- metrics_df %>%
  group_by(window) %>%
  summarise(
    p_te_count = wilcox.test(n_te ~ group)$p.value,
    p_te_fraction = wilcox.test(te_fraction ~ group)$p.value,
    p_nearest_te = wilcox.test(nearest_te_dist ~ group)$p.value,
    .groups = "drop"
  )

write.csv(stats_df, file.path(out_dir, "TE_metrics_stats_by_window.csv"), row.names = FALSE)

# =========================
# VISUALIZATIONS
# =========================
metrics_df<-read.csv(file.path(out_dir, "TE_metrics_nested_windows.csv"))
class_df <-read.csv(file.path(out_dir, "TE_class_counts_nested_windows.csv"))

# 1. TE landscape profile across windows
p_profile <- summary_df %>%
  mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb"))) %>%
  ggplot(aes(x = window, y = mean_te_fraction, group = group, color = group)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Mean TE fraction",
    title = "TE landscape around genes across nested windows",
    subtitle = "Comparison of duplicated-family genes and single-copy genes"
  ) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(file.path(out_dir, "TE_profile_fraction_by_window.pdf"), p_profile, width = 5, height = 3)

# 2. Distribution of TE fraction per locus
library(dplyr)
library(rstatix)

# Calculate significance for each window
signif_data <- metrics_df %>%
  filter(group %in% c("Family", "Single-copy")) %>%
  group_by(window) %>%
  summarise(
    max_y = max(te_fraction, na.rm = TRUE),
    y_pos = log1p(max_y * 1.1)  # 10% above max value, then transformed
  ) %>%
  left_join(
    metrics_df %>%
      filter(group %in% c("Family", "Single-copy")) %>%
      group_by(window) %>%
      wilcox_test(te_fraction ~ group) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance() %>%
      mutate(asterisk = case_when(
        p.adj.signif == "****" ~ "****",
        p.adj.signif == "***" ~ "***",
        p.adj.signif == "**" ~ "**",
        p.adj.signif == "*" ~ "*",
        TRUE ~ ""
      )),
    by = "window"
  ) %>% mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb")))

p_violin <- metrics_df %>%
  mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb"))) %>%
  ggplot(aes(x = group, y = te_fraction, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  facet_wrap(~window, scales = "free_y", nrow = 1) +
  scale_y_continuous(
    trans = "log1p",
    limits = c(0, 3)
  ) +
  # Add significance asterisks
  geom_text(
    data = signif_data,
    aes(x = 1.5, y = y_pos, label = asterisk),  
    size = 6,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "TE fraction per locus (log1p transformed)",
    title = "Distribution of TE occupancy around loci"
  ) + 
  theme(legend.position = "none")

ggsave(file.path(out_dir, "TE_fraction_violin_nested_windows_log.pdf"), p_violin, width = 10, height = 3)

##################### add a truncated y-axis #####################
library(dplyr)
library(rstatix)
library(ggplot2)

window_levels <- c("gene_body", "pm1kb", "pm10kb", "pm25kb", "pm50kb")

metrics_plot <- metrics_df %>%
  filter(group %in% c("Family", "Single-copy")) %>%
  mutate(
    window = factor(window, levels = window_levels),
    te_fraction_log = log1p(te_fraction)
  )

# split/compress y-axis from Q3 to 2
q3_y <- quantile(metrics_plot$te_fraction_log, 0.95, na.rm = TRUE)

break_low  <- q3_y
break_high <- 2.8
gap_size   <- 0.10

compress_y <- function(y) {
  ifelse(
    y <= break_low,
    y,
    ifelse(
      y < break_high,
      break_low + ((y - break_low) / (break_high - break_low)) * gap_size,
      y - (break_high - break_low) + gap_size
    )
  )
}

metrics_plot <- metrics_plot %>%
  mutate(te_fraction_plot = compress_y(te_fraction_log))

# test significance per window
signif_data <- metrics_plot %>%
  group_by(window) %>%
  wilcox_test(te_fraction ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(
    asterisk = case_when(
      p.adj.signif == "****" ~ "****",
      p.adj.signif == "***"  ~ "***",
      p.adj.signif == "**"   ~ "**",
      p.adj.signif == "*"    ~ "*",
      TRUE ~ ""
    ),
    y_pos = compress_y(2.85)
  )

p_violin <- metrics_plot %>%
  ggplot(aes(x = group, y = te_fraction_plot, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
  facet_wrap(~window, scales = "fixed", nrow = 1) +
  scale_y_continuous(
    limits = c(0, compress_y(3)),
    breaks = c(0, break_low, compress_y(2), compress_y(3)),
    labels = c("0", paste0("Q3=", round(break_low, 2)), "2", "3")
  ) +
  geom_hline(yintercept = break_low, linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = break_low + gap_size, linetype = "dashed", linewidth = 0.25) +
  geom_text(
    data = signif_data,
    aes(x = 1.5, y = y_pos, label = asterisk),
    size = 6,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "TE fraction per locus (log1p transformed; compressed above Q3)",
    title = "Distribution of TE occupancy around loci"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "TE_fraction_violin_nested_windows_log_q3_compressed.pdf"),
  p_violin,
  width = 10,
  height = 3
)

p_violin
##################################################################

# 3. TE class enrichment heatmap
class_summary <- class_df %>%
  group_by(group, window, repeat_class) %>%
  summarise(mean_count = mean(n_te_class, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_count, values_fill = 0) %>%
  rename(Single_copy = `Single-copy`) %>%
  mutate(effect = log2((Family + 0.01) / (Single_copy + 0.01)))

p_heat <- class_summary %>%
  mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb")),
         repeat_class = fct_reorder(repeat_class, effect, .fun = max)) %>%
  ggplot(aes(x = window, y = repeat_class, fill = effect)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = "TE class",
    fill = "log2(Fam./Single)",
    title = "TE class enrichment around duplicated-family genes"
  ) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(file.path(out_dir, "TE_class_enrichment_heatmap.pdf"), p_heat, width = 5, height = 3)

# 4. Nearest TE distance
p_dist <- metrics_df %>%
  mutate(window = factor(window, levels = c("gene_body","pm1kb","pm10kb","pm25kb","pm50kb"))) %>%
  ggplot(aes(x = window, y = nearest_te_dist, color = group, group = group)) +
  stat_summary(fun = median, geom = "line", linewidth = 1.2, position = position_dodge(width = 0.1)) +
  stat_summary(fun = median, geom = "point", size = 3, position = position_dodge(width = 0.1)) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0, 9000)) + 
  labs(
    x = NULL,
    y = "Median distance to nearest TE (bp)",
    title = "Distance to nearest TE across nested windows"
  ) + theme(axis.text.x = element_text(angle = 45,hjust = 1))

ggsave(file.path(out_dir, "TE_nearest_distance_profile.pdf"), p_dist, width = 5, height = 3)

