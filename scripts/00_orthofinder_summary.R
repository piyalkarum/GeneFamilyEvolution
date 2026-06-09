############################################################
## OrthoFinder comparative genomics summary figure
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ape)
  library(ggtree)
  library(patchwork)
  library(scales)
})

out_dir <- "outputs/orthofinder_comparative_figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Input files
# -----------------------------
in_dir<-"data/orthofinder_outputs/"
tree_file <- "data/orthofinder_outputs/SpeciesTree_rooted_node_labels.txt"

dup_node_file <- "data/orthofinder_outputs/Duplications_per_Species_Tree_Node.tsv"
dup_og_file   <- "data/orthofinder_outputs/Duplications_per_Orthogroup.tsv"
dup_events_file <- "data/orthofinder_outputs/Duplications.tsv"

stats_species_file <- "data/orthofinder_outputs/Statistics_PerSpecies.tsv"
orth_total_file <- "data/orthofinder_outputs/OrthologuesStats_Totals.tsv"
orth_1to1_file  <- "data/orthofinder_outputs/OrthologuesStats_one-to-one.tsv"
orth_m2m_file   <- "data/orthofinder_outputs/OrthologuesStats_many-to-many.tsv"
orth_1tom_file  <- "data/orthofinder_outputs/OrthologuesStats_one-to-many.tsv"
orth_mto1_file  <- "data/orthofinder_outputs/OrthologuesStats_many-to-one.tsv"
overlap_file    <- "data/orthofinder_outputs/Orthogroups_SpeciesOverlaps.tsv"

# -----------------------------
# Species name map
# -----------------------------
species_map <- c(
  "Aalpina_V4_priTranscripts" = "A. alpina",
  "Ahalleri_v2.1.0.priTranscripts" = "A. halleri",
  "Alyrata_v2.1.priTranscripts" = "A. lyrata",
  "Athaliana_Araport11.priTranscripts" = "A. thaliana",
  "Cpapayav_0.4_priTranscriptOnly" = "C. papaya"
)

species_order <- c("A. alpina", "A. halleri", "A. lyrata", "A. thaliana", "C. papaya")

clean_species <- function(x) {
  out <- species_map[x]
  ifelse(is.na(out), x, out)
}

# ============================================================
# 1. Read species-level statistics
# ============================================================

read_stats_first_block <- function(file) {
  lines <- readLines(file)
  blank <- which(lines == "")[1]
  if (is.na(blank)) blank <- length(lines) + 1
  txt <- paste(lines[1:(blank - 1)], collapse = "\n")
  fread(text = txt, header = TRUE)
}

stats_sp_raw <- read_stats_first_block(stats_species_file)

stats_sp <- stats_sp_raw %>%
  as.data.frame() %>%
  rename(metric = 1) %>%
  pivot_longer(-metric, names_to = "species_raw", values_to = "value") %>%
  mutate(
    species = clean_species(species_raw),
    value = as.numeric(value)
  ) %>%
  select(species_raw, species, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

dup_node <- fread(dup_node_file) %>%
  rename(
    node = `Species Tree Node`,
    dup_all = `Duplications (all)`,
    dup_50 = `Duplications (50% support)`
  )

tip_dup <- dup_node %>%
  filter(node %in% names(species_map)) %>%
  mutate(species = clean_species(node))

stats_sp <- stats_sp %>%
  left_join(tip_dup %>% select(species, terminal_duplications = dup_50), by = "species") %>%
  mutate(
    species = factor(species, levels = species_order),
    pct_genes_in_ogs = `Percentage of genes in orthogroups`,
    pct_unassigned = `Percentage of unassigned genes`,
    pct_species_specific = `Percentage of genes in species-specific orthogroups`,
    terminal_dup_per_1000_genes = terminal_duplications / `Number of genes` * 1000
  )

# ============================================================
# 2. Tree with duplication burden and species metrics
# ============================================================

tr <- read.tree(tree_file)

p_tree_base <- ggtree(tr, size = 0.7) +
  geom_tiplab(aes(label = clean_species(label)), size = 3.5, offset = 0.01) +
  theme_tree2()

tree_dat <- p_tree_base$data

# Add duplication counts for tips and internal nodes
node_dup <- dup_node %>%
  mutate(label = node)

tree_ann <- tree_dat %>%
  left_join(dup_node, by = c("label" = "node")) %>%
  mutate(
    dup_50 = ifelse(is.na(dup_50), 0, dup_50),
    is_dup_node = dup_50 > 0
  )

p_tree <- p_tree_base %<+% tree_ann +
  geom_point(
    data = subset(tree_ann, is_dup_node),
    aes(x = x, y = y, size = dup_50, fill = dup_50),
    shape = 21, color = "black", alpha = 0.85
  ) +
  geom_text(
    data = subset(tree_ann, is_dup_node & !isTip),
    aes(x = x, y = y, label = dup_50),
    nudge_y = 0.22, size = 3
  ) +
  scale_size_continuous(name = "Duplications\n(50% support)", range = c(3, 10)) +
  scale_fill_viridis_c(name = "Duplications\n(50% support)", option = "C") +
  labs(title = "A. Species tree with inferred gene duplication burden") +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

# ============================================================
# 3. Species-level orthogroup composition
# ============================================================

species_comp <- stats_sp %>%
  select(
    species,
    genes = `Number of genes`,
    genes_in_ogs = `Number of genes in orthogroups`,
    unassigned = `Number of unassigned genes`,
    species_specific_genes = `Number of genes in species-specific orthogroups`
  ) %>%
  mutate(
    conserved_or_shared_og_genes = genes_in_ogs - species_specific_genes
  ) %>%
  select(species, conserved_or_shared_og_genes, species_specific_genes, unassigned) %>%
  pivot_longer(-species, names_to = "class", values_to = "n_genes") %>%
  mutate(
    class = recode(
      class,
      conserved_or_shared_og_genes = "Shared/other OG genes",
      species_specific_genes = "Species-specific OG genes",
      unassigned = "Unassigned genes"
    )
  )

p_species_comp <- ggplot(species_comp, aes(x = species, y = n_genes, fill = class)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  theme_classic(base_size = 10) +
  labs(
    title = "B. Gene assignment and species-specific content",
    x = NULL,
    y = "Number of genes",
    fill = NULL
  )

# ============================================================
# 4. Terminal duplication burden per species
# ============================================================

p_dup_bar <- stats_sp %>%
  ggplot(aes(x = species, y = terminal_dup_per_1000_genes)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_classic(base_size = 10) +
  labs(
    title = "C. Terminal duplications normalized by gene number",
    x = NULL,
    y = "Duplications per 1,000 genes"
  )

# ============================================================
# 5. Pairwise orthology heatmaps
# ============================================================

read_pairwise <- function(file, value_name) {
  fread(file) %>%
    as.data.frame() %>%
    rename(sp1_raw = 1) %>%
    pivot_longer(-sp1_raw, names_to = "sp2_raw", values_to = value_name) %>%
    mutate(
      sp1 = clean_species(sp1_raw),
      sp2 = clean_species(sp2_raw),
      !!value_name := as.numeric(.data[[value_name]])
    ) %>%
    filter(sp1 != sp2)
}

orth_total <- read_pairwise(orth_total_file, "total_orthologues")
orth_1to1  <- read_pairwise(orth_1to1_file, "one_to_one")

orth_heat <- orth_total %>%
  left_join(orth_1to1 %>% select(sp1, sp2, one_to_one), by = c("sp1", "sp2")) %>%
  mutate(
    prop_one_to_one = one_to_one / total_orthologues,
    sp1 = factor(sp1, levels = species_order),
    sp2 = factor(sp2, levels = species_order)
  )

p_orth_heat <- ggplot(orth_heat, aes(x = sp1, y = sp2, fill = one_to_one)) +
  geom_tile(color = "white") +
  geom_text(aes(label = comma(round(one_to_one))), size = 2.7) +
  scale_fill_viridis_c(option = "H", labels = comma) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "D. Pairwise one-to-one orthologs",
    x = NULL,
    y = NULL,
    fill = "One-to-one\northologs"
  )

## only the upper part
orth_heat_upper <- orth_heat %>%
  filter(as.integer(sp1) <= as.integer(sp2))

p_orth_heat <- ggplot(orth_heat_upper, aes(x = sp1, y = sp2, fill = one_to_one)) +
  geom_tile(color = "white") +
  geom_text(aes(label = comma(round(one_to_one))), size = 2.7) +
  scale_fill_viridis_c(option = "H", labels = comma) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "D. Pairwise one-to-one orthologs",
    x = NULL,
    y = NULL,
    fill = "One-to-one\northologs"
  )


# ============================================================
# 6. Orthologue class composition
# ============================================================

orth_m2m  <- read_pairwise(orth_m2m_file, "many_to_many")
orth_1tom <- read_pairwise(orth_1tom_file, "one_to_many")
orth_mto1 <- read_pairwise(orth_mto1_file, "many_to_one")

orth_classes <- orth_1to1 %>%
  left_join(orth_1tom, by = c("sp1", "sp2")) %>%
  left_join(orth_mto1, by = c("sp1", "sp2")) %>%
  left_join(orth_m2m, by = c("sp1", "sp2")) %>%
  pivot_longer(
    cols = c(one_to_one, one_to_many, many_to_one, many_to_many),
    names_to = "orthology_class",
    values_to = "count"
  ) %>%
  mutate(
    comparison = paste(sp1, sp2, sep = " vs "),
    orthology_class = recode(
      orthology_class,
      one_to_one = "1:1",
      one_to_many = "1:many",
      many_to_one = "many:1",
      many_to_many = "many:many"
    )
  ) %>%
  group_by(comparison) %>%
  mutate(prop = count / sum(count, na.rm = TRUE)) %>%
  ungroup()

# keep unique unordered pairs for cleaner plot
orth_classes_unique <- orth_classes %>%
  mutate(pair_id = map2_chr(as.character(sp1), as.character(sp2), ~ paste(sort(c(.x, .y)), collapse = " vs "))) %>%
  group_by(pair_id, orthology_class) %>%
  summarise(count = mean(count, na.rm = TRUE), .groups = "drop") %>%
  group_by(pair_id) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

p_orth_comp <- ggplot(orth_classes_unique, aes(x = pair_id, y = prop, fill = orthology_class)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = percent_format()) +
  theme_classic(base_size = 9) +
  labs(
    title = "E. Orthologue class composition",
    x = NULL,
    y = "Proportion of pairwise orthology",
    fill = "Class"
  )

# ============================================================
# 7. Duplication distribution across orthogroups
# ============================================================

dup_og <- fread(dup_og_file) %>%
  rename(
    orthogroup = Orthogroup,
    dup_all = `Duplications (all)`,
    dup_50 = `Duplications (50% support)`
  )

p_dup_dist <- ggplot(dup_og, aes(x = dup_50)) +
  geom_histogram(bins = 60) +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans= "log1p") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "F. Duplication events are concentrated in few orthogroups",
    x = "Duplications per orthogroup (log1p scale)",
    y = "Number of orthogroups (log1p scale)"
  )
## better with log scaled both axes 
# pre-compute bin midpoints and counts to filter zeros
breaks <- seq(min(log1p(dup_og$dup_50)), max(log1p(dup_og$dup_50)), length.out = 61)
h <- hist(log1p(dup_og$dup_50), breaks = breaks, plot = FALSE)

dup_binned <- data.frame(
  xmid = h$mids,
  count = h$counts
) %>%
  filter(count > 0)

p_dup_dist <- ggplot(dup_og, aes(x = log1p(dup_50))) +
  geom_histogram(bins = 60, fill = "grey70", color = "white") +
  geom_smooth(data = dup_binned, aes(x = xmid, y = count), color = 2, linewidth = 0.6) +
  geom_point(data = dup_binned, aes(x = xmid, y = count), color = "black", size = 1.5) +
  scale_x_continuous(
    breaks = log1p(c(0, 1, 2, 5, 10, 20, 50, 100,300)),
    labels = c(0, 1, 2, 5, 10, 20, 50, 100,300)
  ) +
  scale_y_continuous(
    trans = "log1p",
    limits = c(0,9000),
    n.breaks = 50
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "F. Duplication events are concentrated in few orthogroups",
    x = "Duplications per orthogroup (log1p scale)",
    y = "Number of orthogroups\n (log1p scale)"
  )



# ============================================================
# 8. Duplication event support and node type from Duplications.tsv
# ============================================================

dup_events <- fread(dup_events_file, select = c("Orthogroup", "Species Tree Node", "Support", "Type")) %>%
  rename(
    orthogroup = Orthogroup,
    node = `Species Tree Node`,
    support = Support,
    type = Type
  ) %>%
  mutate(
    node_clean = clean_species(node),
    node_clean = ifelse(node_clean %in% names(species_map), clean_species(node), node),
    high_support = support >= 0.5
  )

dup_events_summary <- dup_events %>%
  group_by(node, type) %>%
  summarise(
    n_duplications = n(),
    n_high_support = sum(high_support, na.rm = TRUE),
    median_support = median(support, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  dup_events_summary,
  file.path(out_dir, "duplication_event_summary_by_node.csv"),
  row.names = FALSE
)

p_dup_type <- dup_events_summary %>%
  mutate(
    node_label = clean_species(node),
    node_label = factor(node_label, levels = c(species_order, "N0", "N1", "N2", "N3"))
  ) %>%
  ggplot(aes(x = node_label, y = n_high_support, fill = type)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  theme_classic(base_size = 9) +
  labs(
    title = "G. Duplication events by species-tree node",
    x = NULL,
    y = "Duplications with support ≥ 0.5",
    fill = "Event type"
  )

# ============================================================
# 9. Orthogroup overlap heatmap
# ============================================================

overlap <- fread(overlap_file) %>%
  as.data.frame() %>%
  rename(sp1_raw = 1) %>%
  pivot_longer(-sp1_raw, names_to = "sp2_raw", values_to = "shared_ogs") %>%
  mutate(
    sp1 = factor(clean_species(sp1_raw), levels = species_order),
    sp2 = factor(clean_species(sp2_raw), levels = species_order),
    shared_ogs = as.numeric(shared_ogs)
  )

p_overlap <- ggplot(overlap, aes(x = sp1, y = sp2, fill = shared_ogs)) +
  geom_tile(color = "white") +
  geom_text(aes(label = comma(shared_ogs)), size = 2.6) +
  scale_fill_viridis_c(option = "C", labels = comma) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "H. Shared orthogroups between species",
    x = NULL,
    y = NULL,
    fill = "Shared OGs"
  )

# ============================================================
# 10. Save individual plots
# ============================================================

ggsave(file.path(out_dir, "A_tree_duplication_burden.pdf"), p_tree, width = 6, height = 4)
ggsave(file.path(out_dir, "B_species_gene_assignment.pdf"), p_species_comp, width = 7, height = 4)
ggsave(file.path(out_dir, "C_terminal_duplications_per_gene.pdf"), p_dup_bar, width = 4, height = 2.5)
ggsave(file.path(out_dir, "D_pairwise_one_to_one_heatmap_v1.1.pdf"), p_orth_heat, width = 4, height = 2.5)
ggsave(file.path(out_dir, "E_orthologue_class_composition.pdf"), p_orth_comp, width = 7, height = 5)
ggsave(file.path(out_dir, "F_duplications_per_orthogroup_distribution_v1.1.pdf"), p_dup_dist, width = 4, height = 2.5)
ggsave(file.path(out_dir, "G_duplication_events_by_node.pdf"), p_dup_type, width = 6, height = 4)
ggsave(file.path(out_dir, "H_species_orthogroup_overlap.pdf"), p_overlap, width = 6, height = 5)

# ============================================================
# 11. Combined manuscript figure
# ============================================================

combined_fig <- 
  (p_tree | p_species_comp) /
  (p_dup_bar | p_orth_heat) /
  (p_orth_comp | p_dup_dist) +
  plot_annotation(
    title = "Comparative orthogroup structure and gene duplication dynamics",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(
  file.path(out_dir, "Fig_orthofinder_comparative_summary.pdf"),
  combined_fig,
  width = 13,
  height = 14
)

combined_fig

# ============================================================
# 12. Compact summary table for manuscript
# ============================================================

manuscript_table <- stats_sp %>%
  transmute(
    Species = as.character(species),
    Genes = round(`Number of genes`),
    `Genes in OGs (%)` = round(pct_genes_in_ogs, 1),
    `Species-specific genes (%)` = round(pct_species_specific, 1),
    `OGs containing species` = round(`Number of orthogroups containing species`),
    `Terminal duplications` = round(terminal_duplications),
    `Duplications / 1000 genes` = round(terminal_dup_per_1000_genes, 1)
  )

write.csv(
  manuscript_table,
  file.path(out_dir, "orthofinder_species_summary_for_manuscript.csv"),
  row.names = FALSE
)

print(manuscript_table)


### to add the time calibarted tree 
# ============================================================
# 2. Tree with duplication burden and species metrics
# ============================================================

# ============================================================
# Time-calibrate species tree using A. thaliana split, N2 = 5 Ma
# ============================================================

tr <- read.tree(tree_file)

# find node number for N2
n2_node <- which(tr$node.label == "N2") + length(tr$tip.label)

if (length(n2_node) != 1) {
  stop("Could not uniquely identify node N2 in the tree.")
}

# calibration table: N2 fixed at 5 million years
calib <- data.frame(
  node = n2_node,
  age.min = 5,
  age.max = 5
)

# time-calibrate tree
tr_time <- chronos(
  tr,
  lambda = 1,
  calibration = calib,
  quiet = TRUE
)

# check node ages
branching.times(tr_time)

# ============================================================
# Species tree with time scale and duplication burden
# ============================================================

p_tree_base <- ggtree(tr_time, size = 0.7) +
  geom_tiplab(
    aes(label = clean_species(label)),
    size = 3.5,
    offset = 0.3
  ) +
  theme_tree2() +
  labs(
    x = "Time before present (Ma)"
  )

# reverse time axis so tips are at 0 Ma
p_tree_base <- revts(p_tree_base)

tree_dat <- p_tree_base$data

# clean duplication node names if needed
dup_node <- fread(dup_node_file) %>%
  rename(
    node = `Species Tree Node`,
    dup_all = `Duplications (all)`,
    dup_50 = `Duplications (50% support)`
  ) %>%
  mutate(node = sub("_.*", "", node))

tree_ann <- tree_dat %>%
  left_join(dup_node, by = c("label" = "node")) %>%
  mutate(
    dup_50 = ifelse(is.na(dup_50), 0, dup_50),
    is_dup_node = dup_50 > 0
  )

p_tree <- p_tree_base +
  geom_point(
    data = subset(tree_ann, is_dup_node),
    aes(x = x, y = y, size = dup_50, fill = dup_50),
    shape = 21,
    color = "black",
    alpha = 0.85
  ) +
  geom_text(
    data = subset(tree_ann, is_dup_node & !isTip),
    aes(x = x, y = y, label = dup_50),
    nudge_y = 0.22,
    size = 3
  ) +
  scale_size_continuous(
    name = "Duplications\n(50% support)",
    range = c(3, 10)
  ) +
  scale_fill_viridis_c(
    name = "Duplications\n(50% support)",
    option = "C"
  ) +
  labs(
    title = "A. Time-calibrated species tree with inferred gene duplication burden",
    x = "Time before present (Ma)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

p_tree