library(tidyverse)
library(patchwork)

#--------------------------------------------------
# Input
#--------------------------------------------------
at_file <- "data/231/AT/cnv_within/Athal231_sig_fam_gain_loss_v1.csv"
al_file <- "data/231/AL/cnv_within/Alyrata231_sig_fam_gain_loss_v1.csv"
out_pdf <- "plots/AT_AL_CNV_variable_families_multipage_vertical_mirror.pdf"
out_pdf2 <- "plots/AT_AL_top30_dispersion_families_multipage_vertical_mirror.pdf"

AT_COL <- "#355C7D"
AL_COL <- "#C06C84"

#--------------------------------------------------
# Read data
#--------------------------------------------------
at <- read.csv(at_file, stringsAsFactors = FALSE)
al <- read.csv(al_file, stringsAsFactors = FALSE)

## add function short description
function_short<-data.frame(fread("data/og_fam_function_short.tsv",h=T))
at$func_short<-function_short[match(at$family_id,function_short$Family.ID),"short_title"]
al$func_short<-function_short[match(al$family_id,function_short$Family.ID),"short_title"]

at$species <- "AT"
al$species <- "AL"

keep_cols <- c("family_id", "assembly", "family_size", "gene_count",
               "gain", "loss", "gain_loss_count", "total_gene_count", "species","func_short")

at <- at[, intersect(keep_cols, names(at))]
al <- al[, intersect(keep_cols, names(al))]

df <- bind_rows(at, al)

num_cols <- c("family_size", "gene_count", "gain", "loss", "gain_loss_count", "total_gene_count")
for (cc in intersect(num_cols, names(df))) {
  df[[cc]] <- as.numeric(df[[cc]])
}

#--------------------------------------------------
# Family-level summary stats
#--------------------------------------------------
fam_stats <- df %>%
  group_by(family_id, species,func_short) %>%
  summarise(
    n = n(),
    mean_count = mean(total_gene_count, na.rm = TRUE),
    sd_count = sd(total_gene_count, na.rm = TRUE),
    var_count = var(total_gene_count, na.rm = TRUE),
    cv = ifelse(mean_count > 0, sd_count / mean_count, NA_real_),
    dispersion = ifelse(mean_count > 0, var_count / mean_count, NA_real_),
    min_count = min(total_gene_count, na.rm = TRUE),
    max_count = max(total_gene_count, na.rm = TRUE),
    .groups = "drop"
  )


# keep families variable in at least one species
variable_fams <- fam_stats %>%
  group_by(family_id) %>%
  summarise(
    variable_any = any(var_count > 0, na.rm = TRUE),
    max_cv = max(cv, na.rm = TRUE),
    max_disp = max(dispersion, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(variable_any) %>%
  arrange(desc(max_cv), desc(max_disp)) %>%
  pull(family_id)

#--------------------------------------------------
# Select top 30 families by dispersion across species
#--------------------------------------------------
### rerun from this part to plot only top 30 fams by dispersion 
top_n_fams <- 40

variable_fams <- fam_stats %>%
  group_by(family_id) %>%
  summarise(
    variable_any = any(var_count > 0, na.rm = TRUE),
    max_disp = max(dispersion, na.rm = TRUE),
    max_cv = max(cv, na.rm = TRUE),
    mean_disp = mean(dispersion, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(variable_any, is.finite(max_disp)) %>%
  arrange(desc(max_disp), desc(max_cv), desc(mean_disp)) %>%
  slice_head(n = top_n_fams) %>%
  pull(family_id)


#--------------------------------------------------
# Distribution table
#--------------------------------------------------
plot_df <- df %>%
  filter(family_id %in% variable_fams) %>%
  count(family_id, species,total_gene_count, name = "freq")

all_counts <- plot_df %>%
  group_by(family_id) %>%
  summarise(
    min_copy = min(total_gene_count, na.rm = TRUE),
    max_copy = max(total_gene_count, na.rm = TRUE),
    .groups = "drop"
  )

fam_annot <- df %>%
  distinct(family_id, func_short)

plot_df_full <- all_counts %>%
  rowwise() %>%
  do({
    tibble(
      family_id = .$family_id,
      total_gene_count = seq(.$min_copy, .$max_copy)
    )
  }) %>%
  ungroup() %>%
  crossing(species = c("AT", "AL")) %>%
  left_join(plot_df, by = c("family_id", "species", "total_gene_count")) %>%
  left_join(fam_annot, by = "family_id") %>%
  mutate(
    freq = replace_na(freq, 0),
    freq_plot = ifelse(species == "AL", -freq, freq)
  )

fam_stats_wide <- fam_stats %>%
  filter(family_id %in% variable_fams) %>%
  select(family_id, species, cv, dispersion) %>%
  pivot_wider(
    names_from = species,
    values_from = c(cv, dispersion),
    names_sep = "_"
  )

#--------------------------------------------------
# One mirrored plot per family
# AT upward, AL downward
# Each species scaled independently within family
#--------------------------------------------------
plot_one_family_vertical_mirror <- function(fid, plot_df_full, fam_stats_wide,
                                            at_col = AT_COL, al_col = AL_COL) {
  
  tm <- plot_df_full %>%
    filter(family_id == fid)
  
  func_lab <- tm %>%
    filter(!is.na(func_short)) %>%
    distinct(func_short) %>%
    pull(func_short)
  
  func_lab <- ifelse(length(func_lab) == 0, fid, func_lab[1])
  
  st <- fam_stats_wide %>% filter(family_id == fid)
  
  # Get total number of assemblies for each species
  n_at <- n_distinct(df$assembly[df$species == "AT"])
  n_al <- n_distinct(df$assembly[df$species == "AL"])
  
  # Convert frequencies to percentages
  tm <- tm %>%
    group_by(species) %>%
    mutate(
      freq_pct = freq / ifelse(species == "AT", n_at, n_al) * 100,
      freq_plot = ifelse(species == "AL", -freq_pct, freq_pct)
    ) %>%
    ungroup()
  
  max_at <- max(tm$freq_pct[tm$species == "AT"], na.rm = TRUE)
  max_al <- max(tm$freq_pct[tm$species == "AL"], na.rm = TRUE)
  
  # avoid zero-range problems
  if (!is.finite(max_at) || max_at == 0) max_at <- 1
  if (!is.finite(max_al) || max_al == 0) max_al <- 1
  
  subtitle_txt <- paste0(
    "AT: CV=", sprintf("%.2f", st$cv_AT),
    ", Disp=", sprintf("%.2f", st$dispersion_AT),
    "    |    AL: CV=", sprintf("%.2f", st$cv_AL),
    ", Disp=", sprintf("%.2f", st$dispersion_AL)
  )
  
  ggplot(tm, aes(x = factor(total_gene_count), y = freq_plot, fill = species)) +
    geom_col(width = 0.8, color = NA) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
    scale_fill_manual(values = c("AT" = at_col, "AL" = al_col)) +
    scale_y_continuous(
      labels = function(x) abs(x),
      limits = c(-max(max_al * 1.18, max_at * 0.2), 
                 max(max_at * 1.18, max_al * 0.2))
    ) +
    labs(
      title = paste0(func_lab, " (", fid, ")"), #fid,
      subtitle = subtitle_txt,
      x = "Copy number",
      y = "Assemblies (%)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, hjust = 0.5),
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

#--------------------------------------------------
# Build family plots
#--------------------------------------------------
family_plots <- lapply(variable_fams, function(fid) {
  plot_one_family_vertical_mirror(fid, plot_df_full, fam_stats_wide)
})
names(family_plots) <- variable_fams

#--------------------------------------------------
# Split into pages: 15 families per page
# layout = 5 rows x 3 cols
#--------------------------------------------------
split_into_chunks <- function(x, chunk_size = 15) {
  split(x, ceiling(seq_along(x) / chunk_size))
}

plot_chunks <- split_into_chunks(family_plots, chunk_size = 20)

#--------------------------------------------------
# Output multipage PDF
#--------------------------------------------------
pdf(out_pdf, width = 12, height = 16, onefile = TRUE)

for (i in seq_along(plot_chunks)) {
  chunk <- plot_chunks[[i]]
  
  if (length(chunk) < 10) {
    n_missing <- 10 - length(chunk)
    blank_plot <- ggplot() + theme_void()
    chunk <- c(chunk, replicate(n_missing, blank_plot, simplify = FALSE))
  }
  
  page_plot <- wrap_plots(chunk, ncol = 3, nrow = 5)
  print(page_plot)
}

dev.off()

## only for top 30 families 
pdf(out_pdf2, width = 12, height = 16, onefile = TRUE)

for (i in length(plot_chunks):1) {
  chunk <- plot_chunks[[i]]
  
  if (length(chunk) < 20) {
    n_missing <- 20 - length(chunk)
    blank_plot <- ggplot() + theme_void()
    chunk <- c(chunk, replicate(n_missing, blank_plot, simplify = FALSE))
  }
  
  page_plot <- wrap_plots(chunk, ncol = 4, nrow = 5)
  print(page_plot)
}

dev.off()






#### RDA plots ########

## helpers ##
match_sample_names <- function(site_names, cluster_names) {
  # Create a lookup table with cleaned names
  lookup <- data.frame(
    original = cluster_names,
    cleaned = clean_sample_name(cluster_names),
    stringsAsFactors = FALSE
  )
  
  # Match each site name to cluster name
  matches <- sapply(site_names, function(site_name) {
    site_clean <- clean_sample_name(site_name)
    
    # Try exact match on cleaned names
    match_idx <- which(lookup$cleaned == site_clean)
    
    if (length(match_idx) == 1) {
      return(lookup$original[match_idx])
    }
    
    # Try partial matching (for cases like "Cvi-0" vs "Cvi-0.6911")
    match_idx <- grep(paste0("^", site_clean), lookup$cleaned)
    if (length(match_idx) == 1) {
      return(lookup$original[match_idx])
    }
    
    # Try matching without special characters
    site_simple <- gsub("[[:punct:]]", "", site_clean)
    match_idx <- which(gsub("[[:punct:]]", "", lookup$cleaned) == site_simple)
    if (length(match_idx) == 1) {
      return(lookup$original[match_idx])
    }
    
    return(NA)
  })
  
  return(matches)
}

clean_sample_name <- function(name) {
  # Remove parentheses and their contents
  name <- gsub("\\([^)]*\\)", "", name)
  # Remove spaces before numbers/suffixes
  name <- gsub("\\s+(\\d)", "\\1", name)
  # Standardize dots and spaces
  name <- gsub("[[:space:]]+", " ", name)
  name <- trimws(name)
  # Handle specific cases
  name <- gsub("^Bön\\s*1$", "Bön", name)
  name <- gsub("^Cvi-0\\s*", "Cvi-0", name)
  name <- gsub("\\.\\d+$", "", name)  # Remove .6911 suffixes
  return(name)
}

plot_rda_ordination <- function(res_obj, pop_col = NULL, clusters=NULL, max_arrows = 8) {
  site_df <- res_obj$sites
  bp_df <- res_obj$biplot
  
  if (!is.null(bp_df) && nrow(bp_df) > 0) {
    bp_df <- bp_df %>%
      mutate(arrow_len = sqrt(RDA1^2 + RDA2^2)) %>%
      arrange(desc(arrow_len)) %>%
      slice(1:min(max_arrows, n()))
  }
  
  site_df$clean_name <- clean_sample_name(unlist(site_df[pop_col]))
  
  # exact match between assembly and cluster sample name
  cl_sub <- clusters[match(site_df$clean_name, clusters$Sample), , drop = FALSE]
  
  site_df$geo_group <- cl_sub$Group
  site_df$gen_cluster <- factor(cl_sub$LR_cluster)
  
  # fallback labels for unmatched rows
  site_df$geo_group[is.na(site_df$geo_group)] <- "Unknown"
  site_df$gen_cluster[is.na(site_df$gen_cluster)] <- "Unknown"
  site_df$gen_cluster <- factor(site_df$gen_cluster)
  
  # variance explained
  constrained_eig <- res_obj$model$CCA$eig
  unconstrained_eig <- res_obj$model$CA$eig
  total_variance <- sum(constrained_eig, unconstrained_eig)
  
  # Percent variance for first two RDA axes
  var_exp_rda1 <- round(constrained_eig[1] / total_variance * 100, 1)
  var_exp_rda2 <- round(constrained_eig[2] / total_variance * 100, 1)
  
  p <- ggplot(site_df, aes(RDA1, RDA2)) +
    geom_point(
      aes(
        color = geo_group,
        shape = gen_cluster
      ),
      size = 1.5,
      stroke = 1
    ) + 
    theme_classic(base_size = 11) +
    labs(
      title = res_obj$label,
      x = paste0("RDA1 (", var_exp_rda1, "%)"),
      y = paste0("RDA2 (", var_exp_rda2, "%)")
    ) +
    scale_shape_manual(
      values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
      name = "Genetic Cluster"
    ) +
    scale_color_brewer(
      palette = "Set1",
      name = "Geographic Group"
    ) +
    ggrepel::geom_text_repel(
      aes(label = clean_name),
      size = 3,
      color = "black",
      max.overlaps = 10,
      box.padding = 0.35,
      point.padding = 0.2,
      segment.color = "grey60",
      segment.size = 0.3
    )
  
  if (!is.null(bp_df) && nrow(bp_df) > 0) {
    mult <- 1.15
    p <- p +
      geom_segment(
        data = bp_df,
        aes(x = 0, y = 0, xend = RDA1 * mult, yend = RDA2 * mult),
        inherit.aes = FALSE,
        arrow = arrow(length = unit(0.18, "cm")),
        color = "black"
      ) +
      geom_text(
        data = bp_df,
        aes(x = RDA1 * mult, y = RDA2 * mult, label = variable),
        inherit.aes = FALSE,
        size = 3,
        color=2,
        hjust = 0.5, vjust = -0.3
      )
  }
  
  p
}


## gene family functional category plot ----------------------
################################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)}
)


# -----------------------------
# read data
# -----------------------------
fun_added <- data.frame(fread("data/ogList_with_functions_v2.csv"))

species_cols <- c("Aalpina", "Ahalleri", "Alyrata", "Athaliana")

# -----------------------------
# 1. family-category percentages
# -----------------------------
fam_summary <- fun_added %>%
  count(FunctionalCategory, name = "FamilyCount") %>%
  mutate(
    FamilyPercent = 100 * FamilyCount / sum(FamilyCount)
  )

# -----------------------------
# 2. gene counts per species per category
# -----------------------------
bar_df <- fun_added[, c(species_cols, "FunctionalCategory")]

bar_long <- bar_df %>%
  pivot_longer(
    cols = all_of(species_cols),
    names_to = "Species",
    values_to = "GeneCount"
  )

gene_summary <- bar_long %>%
  group_by(FunctionalCategory, Species) %>%
  summarise(TotalGenes = sum(GeneCount, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 3. secondary category summary
#    include NA in denominator,
#    but exclude NA from plotting only
# -----------------------------
sec_all <- fun_added %>%
  count(FunctionalCategory, FunctionalCategory1, name = "n_secondary")

# denominator includes ALL (including NA)
sec_denoms <- sec_all %>%
  group_by(FunctionalCategory) %>%
  summarise(total_secondary = sum(n_secondary), .groups = "drop")

sec_summary <- sec_all %>%
  left_join(sec_denoms, by = "FunctionalCategory") %>%
  mutate(
    prop_secondary = n_secondary / total_secondary
  ) %>%
  # remove only NA from plotting
  filter(!is.na(FunctionalCategory1)) %>%
  group_by(FunctionalCategory) %>%
  arrange(desc(n_secondary), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  ungroup()


# -----------------------------
# 4. order categories by family percentage
# -----------------------------
cat_levels <- fam_summary %>%
  arrange(FamilyPercent) %>%
  pull(FunctionalCategory)

fam_summary$FunctionalCategory <- factor(fam_summary$FunctionalCategory, levels = cat_levels)
gene_summary$FunctionalCategory <- factor(gene_summary$FunctionalCategory, levels = cat_levels)
sec_summary$FunctionalCategory <- factor(sec_summary$FunctionalCategory, levels = cat_levels)

fam_summary <- fam_summary %>%
  mutate(y = as.numeric(FunctionalCategory))

# -----------------------------
# 5. identify top contributing species per category
# -----------------------------
top_species <- gene_summary %>%
  group_by(FunctionalCategory) %>%
  slice_max(order_by = TotalGenes, n = 1, with_ties = FALSE) %>%
  ungroup()

# -----------------------------
# 6. rescale gene counts to same x-axis
# -----------------------------
max_percent <- max(fam_summary$FamilyPercent, na.rm = TRUE)
max_genes   <- max(gene_summary$TotalGenes, na.rm = TRUE)

gene_summary <- gene_summary %>%
  mutate(
    GeneScaled = TotalGenes / max_genes * max_percent
  )

top_species <- top_species %>%
  mutate(
    GeneScaled = TotalGenes / max_genes * max_percent
  )

# -----------------------------
# 7. manual y-offsets for lollipops
# -----------------------------
species_levels <- c("Aalpina", "Ahalleri", "Alyrata", "Athaliana")
species_offsets <- c(-0.24, -0.08, 0.08, 0.24)

gene_summary <- gene_summary %>%
  mutate(
    Species = factor(Species, levels = species_levels),
    y_base = as.numeric(FunctionalCategory),
    y = y_base + species_offsets[as.numeric(Species)]
  )

top_species <- top_species %>%
  mutate(
    Species = factor(Species, levels = species_levels),
    y_base = as.numeric(FunctionalCategory),
    y = y_base + species_offsets[as.numeric(Species)]
  )

# -----------------------------
# 8. build right-side secondary-category tile strip
# -----------------------------
# place tiles just to the right of the % labels
tile_x_start <- max_percent + 6
tile_width <- 1.2

sec_summary <- sec_summary %>%
  group_by(FunctionalCategory) %>%
  arrange(desc(n_secondary), .by_group = TRUE) %>%
  mutate(tile_id = row_number()) %>%
  ungroup() %>%
  mutate(
    y = as.numeric(FunctionalCategory),
    xmin = tile_x_start + (tile_id - 1) * tile_width,
    xmax = tile_x_start + tile_id * tile_width - 0.1,
    ymin = y - 0.28,
    ymax = y + 0.28
  )

# -----------------------------
# 9. plot
# -----------------------------
cat_lebels<-gsub("and","&",cat_levels)
cat_lebels<-gsub("& ","&\n",cat_lebels)

p <- ggplot() +
  # background horizontal bars = percentage of gene families
  geom_rect(
    data = fam_summary,
    aes(
      xmin = 0,
      xmax = FamilyPercent,
      ymin = y - 0.36,
      ymax = y + 0.36
    ),
    fill = "grey85",
    color = NA
  ) +
  # family percentage labels
  geom_text(
    data = fam_summary,
    aes(
      x = FamilyPercent + 0.8,
      y = y,
      label = paste0(round(FamilyPercent, 1), "%")
    ),
    size = 3.2,
    hjust = 0
  ) +
  # lollipop stems
  geom_segment(
    data = gene_summary,
    aes(x = 0, xend = GeneScaled, y = y, yend = y, color = Species),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  # lollipop heads
  geom_point(
    data = gene_summary,
    aes(x = GeneScaled, y = y, color = Species),
    size = 2.8
  ) +
  # label top species contribution per category
  geom_text(
    data = top_species,
    aes(x = GeneScaled, y = y, label = TotalGenes, color = Species),
    size = 3,
    fontface = "bold",
    nudge_x = 0.8,
    show.legend = FALSE
  ) +
  # secondary-category tile strip
  geom_rect(
    data = sec_summary,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = FunctionalCategory1
    ),
    color = "white",
    linewidth = 0.2
  ) +
  # optional percentages inside tiles
  geom_text(
    data = sec_summary,
    aes(
      x = (xmin + xmax) / 2,
      y = y,
      label = paste0(round(prop_secondary * 100), "%")
    ),
    size = 3,
    color = "black",
    angle = 90
  ) +
  scale_y_continuous(
    breaks = seq_along(cat_levels),
    labels = cat_lebels,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_x_continuous(
    name = "Gene families (%)",
    expand = expansion(mult = c(0, 0.02)),
    limits = c(0, tile_x_start + 3.8),
    sec.axis = sec_axis(
      trans = ~ . / max_percent * max_genes,
      name = "Number of genes"
    )
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Set3", na.translate = FALSE) +
  labs(
    y = NULL,
    color = "Species",
    fill = "Secondary function",
    title = "Functional composition of gene families and species-specific gene contribution",
    subtitle = "Grey bars: primary family category frequency\nLollipops: gene counts per species\nRight tiles: dominant secondary categories"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 40, 10, 10)
  )

print(p)

ggsave(
  "plots/functional_categories_combined_plot_with_secondary_v1.1.pdf",
  p,
  width = 9,
  height = 9
)



### GO enrichment assessment v2 ----------------------
library(clusterProfiler)
library(org.At.tair.db)

fun_added <- data.frame(fread("data/ogList_with_functions_v2.csv"))
my_genes <- stringr::str_trim(stringr::str_split_fixed(unlist(strsplit(fun_added$at_genes,split=",")),pattern = "\\.",n=2)[,1],"both")

# Run GO enrichment
ego <- enrichGO(gene          = my_genes,
                universe      = keys(org.At.tair.db),
                OrgDb         = org.At.tair.db,
                keyType       = "TAIR",
                ont           = "ALL",         # All Processes (Biological, cellular, Molecular)
                pAdjustMethod = "BH",         # Benjamini-Hochberg correction
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.10)

# Extract results and filter for defense/stress
ego_df <- as.data.frame(ego)

# View results
dotplot(ego, showCategory = 50)

# search and add functional categories to GO table
go_ids<-stringr::str_split(ego_df$geneID,pattern = "\\/")

fun_cats<-lapply(go_ids,function(x){
  tm<-NULL
  for(i in seq_along(x)){
    tm<-c(tm,fun_added[grep(x[i],fun_added$at_genes),"FunctionalCategory"])
    vv<-table(tm)
    tm<-names(which.max(vv))
  }
  return(tm)
})


library(dplyr)
library(ggplot2)
library(tidytext)

ego_df <- as.data.frame(ego)
ego_df$Functional_category<-unlist(fun_cats)


tv<-data.frame(stringr::str_split_fixed(ego_df$GeneRatio,"\\/",n=2))
tv$X1<-as.numeric(tv$X1);tv$X2<-as.numeric(tv$X2)
ego_df$GeneRatio<-tv$X1/tv$X2
write.csv(ego_df,"data/GO_raw_output.csv",row.names=F)

# import GO short form added table
ego_df0<-read.csv("data/GO_raw_output_with_short_descriptions.csv")

ego_df$Short_description<-ego_df0$Short_description

# Create a numbered ID for each GO term (e.g., "GO_1", "GO_2")
ego_df <- ego_df %>%
  group_by(Functional_category) %>%
  arrange(p.adjust) %>%  # Sort by significance within each category
  mutate(
    GO_short = Short_description,  # Short label
    Description_full = Description            # Store full description
  ) %>%
  ungroup()

# Maintain order for plotting
ego_df$GO_short <- factor(ego_df$GO_short, levels = unique(ego_df$GO_short))

ego_df2 <- ego_df %>%
  mutate(GO_lab = as.character(GO_short))


facet_labels <- c(
  "Defense & Stress response" = "DS",
  "Genetic Information & Regulation" = "GI",
  "Biosynthesis & Metabolism" = "BM",
  "Growth & Reproduction" = "GR",
  "Protein Processing and Degradation" = "PP",
  "Multifunction" = "MF"
)

pp<-ggplot(
  ego_df2,
  aes(
    x = GeneRatio,
    y = reorder_within(GO_lab, -p.adjust, Functional_category),
    color = p.adjust,
    size = Count
  )
) +
  geom_point() +
  facet_grid(
    Functional_category ~ .,
    scales = "free_y",
    space = "free_y",
    labeller = labeller(Functional_category = facet_labels)
  ) +
  # facet_grid(Functional_category ~ ., scales = "free_y", space = "free_y") +
  scale_y_reordered() +
  scale_color_gradient(
    low = "red", high = "blue",
    name = "Adjusted\n p-value",
    trans = "log10"
  ) +
  scale_size_continuous(name = "Gene count") +
  labs(x = "Gene Ratio", y = "GO Term") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(angle = 90)
  ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "plots/GO_enrichment_with_short_descriptions_v2.pdf",
  pp,
  width = 4.75,
  height = 8
)

#==================================================================

## Pseudogenization/functional changes plot -----------------------
### THALIANA --------------

suppressPackageStartupMessages({library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(tidyr)})

# -----------------------------
# input
# -----------------------------
ps_dir  <- "outputs/AT/function_fate/step5"
infile  <- file.path(ps_dir, "gene_status_proportions.tsv")
out_pdf <- file.path(ps_dir, "Fig_familywise_boxplots_PFU_corrected.pdf")
out_tsv <- file.path(ps_dir, "familywise_boxplot_medians_corrected.tsv")

df <- read.delim(infile, check.names = FALSE)
ogList<-data.frame(fread("data/ogList_with_functions_v2.csv"))

df<- df%>% left_join(ogList %>% select(Family.ID, FunctionalCategory),
                     by=c("family_id"="Family.ID"))

# First, aggregate to get proportions per gene family
family_summary <- df %>%
  group_by(family_id, final_status) %>%
  summarise(total_n = sum(n), .groups = 'drop') %>%
  group_by(family_id) %>%
  mutate(proportion = total_n / sum(total_n))

# Color palette for statuses
status_colors <- c("P" = "#E41A1C", "F" = "#377EB8", "U" = "#4DAF4A", "A" = "grey70")

# Calculate proportions per gene per status
gene_level_data <- df %>%
  group_by(family_id, gene, final_status) %>%
  summarise(prop = unique(prop), .groups = 'drop')

# Or aggregate to family-level proportions
family_props <- df %>%
  group_by(family_id, final_status) %>%
  summarise(total_genes = sum(n), .groups = 'drop') %>%
  group_by(family_id) %>%
  mutate(proportion = total_genes / sum(total_genes))

p_distribution <- ggplot(family_props, aes(x = final_status, y = proportion, fill = final_status)) +
  geom_violin(trim = FALSE, alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.8, color = "black") +
  scale_fill_manual(values = status_colors,
                    labels = c("P" = "Pseudogenization (P)",
                               "F" = "Functional changes (F)",
                               "U" = "No functional change (U)",
                               "A" = "Other (A)")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Functional status", 
       y = "Proportion per gene family",
       title = "Distribution of functional status across all gene families",
       subtitle = "Each point represents one gene family") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.major.y = element_line(color = "gray90", size = 0.3))

# Add statistics
library(ggpubr)
p_distribution <- p_distribution + 
  stat_compare_means(comparisons = list(c("P", "F"), c("P", "U"), c("F", "U")),
                     method = "wilcox.test",
                     label = "p.signif",
                     tip_length = 0.01)

ggsave("plots/AT_functional_status_distribution.pdf", p_distribution, width = 6, height = 4)


#### stacked bars with size comparison ###########

# Calculate correct family size (unique genes per family)
family_sizes <- df %>%
  group_by(family_id) %>%
  summarise(family_size = n_distinct(gene), .groups = 'drop')

# Aggregate to get proportions per gene family
family_summary <- df %>%
  group_by(family_id, final_status) %>%
  summarise(total_n = sum(n), .groups = 'drop') %>%
  group_by(family_id) %>%
  mutate(proportion = total_n / sum(total_n)) %>%
  left_join(family_sizes, by = "family_id")


write.table(family_summary,file.path(out_tsv),quote=F,sep="\t",row.names = F)

# Color palette for statuses
status_colors <- c("P" = "#E41A1C", "F" = "#377EB8", "U" = "#4DAF4A", "A" = "grey70")

# Plot 1: Family sizes (ascending order)
p_size <- ggplot(family_sizes, aes(x = reorder(family_id, family_size), y = family_size)) +
  geom_bar(stat = "identity", fill = "gray50", color = "black", width = 0.7, size = 0.2) +
  labs(x = NULL, 
       y = "Number of unique genes",
       title = "Gene family sizes") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray90", size = 0.3))

# Plot 2: Stacked proportions (same order as top plot)
p_stacked <- ggplot(family_summary, aes(x = reorder(family_id, family_size), 
                                        y = proportion, 
                                        fill = final_status)) +
  geom_bar(stat = "identity", width = 0.7, color = NA, size = 0.2) +
  scale_fill_manual(values = status_colors,
                    labels = c("P" = "Pseudogenization",
                               "F" = "Functional changes",
                               "U" = "No functional change",
                               "A" = "Ambiguous")) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Gene families (ordered by size, smallest to largest)", 
       y = "Proportion of genes",
       title = "Distribution of functional status",
       fill = "Status") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "gray90", size = 0.3)) +
  guides(fill = guide_legend(nrow = 2))


# plot 3: sorted by functional category then by pseudogenization
family_pseudo_for_sorting <- df %>%
  group_by(family_id) %>%
  summarise(
    pseudo_prop = sum(n[final_status == "P"]) / sum(n),
    .groups = 'drop'
  )
family_order <- df %>%
  distinct(family_id, FunctionalCategory) %>%
  left_join(family_pseudo_for_sorting, by = "family_id") %>%
  arrange(FunctionalCategory, desc(pseudo_prop)) %>%
  pull(family_id)

p_stacked_by_category <- family_summary %>%
  left_join(df %>% distinct(family_id, FunctionalCategory), by = "family_id") %>%
  mutate(
    family_id = factor(family_id, levels = family_order),
    # Create combined label for x-axis
    x_label = paste0(FunctionalCategory, "\n(", family_id, ")")
  ) %>%
  ggplot(aes(x = family_id, y = proportion, fill = final_status)) +
  geom_bar(stat = "identity", width = 0.7, color = NA, size = 0.2) +
  # Add subtle vertical lines between categories
  geom_vline(
    data = . %>% 
      group_by(FunctionalCategory) %>% 
      summarise(x = max(as.numeric(family_id)) + 0.5, .groups = 'drop') %>%
      filter(!is.na(x)),
    aes(xintercept = x),
    linetype = "dotted", color = "gray20", alpha = 1
  ) +
  scale_fill_manual(values = status_colors,
                    labels = c("P" = "Pseudogenization",
                               "F" = "Functional changes",
                               "U" = "No functional change",
                               "A" = "Ambiguous")) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Gene families (sorted by functional category, then by pseudogenization rate)", 
       y = "Proportion of genes",
       title = "Distribution of functional status by category",
       fill = "Status") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "gray90", size = 0.3)) +
  guides(fill = guide_legend(nrow = 2))


# Combine both plots
combined_plot <-  p_stacked / p_stacked_by_category+ #p_size /
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

# Save
ggsave("plots/AT_functional_status_with_sizes_v2.pdf", combined_plot, width = 8, height = 4)


### LYRATA --------------

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(tidyr)

# -----------------------------
# input
# -----------------------------
ps_dir  <- "outputs/AL/function_fate/step5"
infile  <- file.path(ps_dir, "gene_status_proportions.tsv")
out_pdf <- file.path(ps_dir, "AL_Fig_familywise_boxplots_PFU_corrected.pdf")
out_tsv <- file.path(ps_dir, "AL_familywise_boxplot_medians_corrected.tsv")

df <- read.delim(infile, check.names = FALSE)
ogList<-data.frame(fread("data/ogList_with_functions_v2.csv"))

df<- df%>% left_join(ogList %>% select(Family.ID, FunctionalCategory),
                     by=c("family_id"="Family.ID"))


family_summary <- df %>%
  group_by(family_id, final_status) %>%
  summarise(total_n = sum(n), .groups = 'drop') %>%
  group_by(family_id) %>%
  mutate(proportion = total_n / sum(total_n))

# Color palette for statuses
status_colors <- c("P" = "#E41A1C", "F" = "#377EB8", "U" = "#4DAF4A", "A" = "grey70")

# Calculate proportions per gene per status
gene_level_data <- df %>%
  group_by(family_id, gene, final_status) %>%
  summarise(prop = unique(prop), .groups = 'drop')

# Or aggregate to family-level proportions
family_props <- df %>%
  group_by(family_id, final_status) %>%
  summarise(total_genes = sum(n), .groups = 'drop') %>%
  group_by(family_id) %>%
  mutate(proportion = total_genes / sum(total_genes))

p_distribution <- ggplot(family_props, aes(x = final_status, y = proportion, fill = final_status)) +
  geom_violin(trim = FALSE, alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.8, color = "black") +
  scale_fill_manual(values = status_colors,
                    labels = c("P" = "Pseudogenization (P)",
                               "F" = "Functional changes (F)",
                               "U" = "No functional change (U)",
                               "A" = "Other (A)")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Functional status", 
       y = "Proportion per gene family",
       title = "Distribution of functional status across all gene families",
       subtitle = "Each point represents one gene family") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.major.y = element_line(color = "gray90", size = 0.3))

# Add statistics
library(ggpubr)
p_distribution <- p_distribution + 
  stat_compare_means(comparisons = list(c("P", "F"), c("P", "U"), c("F", "U")),
                     method = "wilcox.test",
                     label = "p.signif",
                     tip_length = 0.01)

ggsave("plots/AL_functional_status_distribution.pdf", p_distribution, width = 6, height = 4)


#### stacked bars with size comparison ###########

# Calculate correct family size (unique genes per family)
family_sizes <- df %>%
  group_by(family_id) %>%
  summarise(family_size = n_distinct(gene), .groups = 'drop')

# Aggregate to get proportions per gene family
family_summary <- df %>%
  group_by(family_id, final_status) %>%
  summarise(total_n = sum(n), .groups = 'drop') %>%
  group_by(family_id) %>%
  mutate(proportion = total_n / sum(total_n)) %>%
  left_join(family_sizes, by = "family_id")

write.table(family_summary,file.path(out_tsv),quote=F,sep="\t",row.names = F)

# Color palette for statuses
status_colors <- c("P" = "#E41A1C", "F" = "#377EB8", "U" = "#4DAF4A", "A" = "grey70")

# Plot 1: Family sizes (ascending order)
p_size <- ggplot(family_sizes, aes(x = reorder(family_id, family_size), y = family_size)) +
  geom_bar(stat = "identity", fill = "gray50", color = "black", width = 0.7, size = 0.2) +
  labs(x = NULL, 
       y = "Number of unique genes",
       title = "Gene family sizes") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray90", size = 0.3))

# Plot 2: Stacked proportions (same order as top plot)
p_stacked <- ggplot(family_summary, aes(x = reorder(family_id, family_size), 
                                        y = proportion, 
                                        fill = final_status)) +
  geom_bar(stat = "identity", width = 0.7, color = NA, size = 0.2) +
  scale_fill_manual(values = status_colors,
                    labels = c("P" = "Pseudogenization",
                               "F" = "Functional changes",
                               "U" = "No functional change",
                               "A" = "Ambiguous")) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Gene families (ordered by size, smallest to largest)", 
       y = "Proportion of genes",
       title = "Distribution of functional status",
       fill = "Status") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "gray90", size = 0.3)) +
  guides(fill = guide_legend(nrow = 2))

# plot 3: sorted by functional category then by pseudogenization
family_pseudo_for_sorting <- df %>%
  group_by(family_id) %>%
  summarise(
    pseudo_prop = sum(n[final_status == "P"]) / sum(n),
    .groups = 'drop'
  )
family_order <- df %>%
  distinct(family_id, FunctionalCategory) %>%
  left_join(family_pseudo_for_sorting, by = "family_id") %>%
  arrange(FunctionalCategory, desc(pseudo_prop)) %>%
  pull(family_id)

p_stacked_by_category <- family_summary %>%
  left_join(df %>% distinct(family_id, FunctionalCategory), by = "family_id") %>%
  mutate(
    family_id = factor(family_id, levels = family_order),
    # Create combined label for x-axis
    x_label = paste0(FunctionalCategory, "\n(", family_id, ")")
  ) %>%
  ggplot(aes(x = family_id, y = proportion, fill = final_status)) +
  geom_bar(stat = "identity", width = 0.7, color = NA, size = 0.2) +
  # Add subtle vertical lines between categories
  geom_vline(
    data = . %>% 
      group_by(FunctionalCategory) %>% 
      summarise(x = max(as.numeric(family_id)) + 0.5, .groups = 'drop') %>%
      filter(!is.na(x)),
    aes(xintercept = x),
    linetype = "dotted", color = "gray20", alpha = 1
  ) +
  scale_fill_manual(values = status_colors,
                    labels = c("P" = "Pseudogenization",
                               "F" = "Functional changes",
                               "U" = "No functional change",
                               "A" = "Ambiguous")) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Gene families (sorted by functional category, then by pseudogenization rate)", 
       y = "Proportion of genes",
       title = "Distribution of functional status by category",
       fill = "Status") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "gray90", size = 0.3)) +
  guides(fill = guide_legend(nrow = 2))


# Combine both plots
combined_plot <- p_stacked / p_stacked_by_category + 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

# Save
ggsave("plots/AL_functional_status_with_sizes_v2.pdf", combined_plot, width = 8, height = 4)

#========================================================================================



