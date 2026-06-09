
convert_long <- function(df, species_name) {
  df %>%
    pivot_longer(-family_id, names_to = "assembly", values_to = "copy") %>%
    mutate(species = species_name)
}

make_full_long <- function(sp1, sp2, species1_name, species2_name) {
  df1 <- convert_long(sp1, species1_name)
  df2 <- convert_long(sp2, species2_name)
  
  all_families <- union(sp1$family_id, sp2$family_id)
  all_assemblies_sp1 <- setdiff(colnames(sp1), "family_id")
  all_assemblies_sp2 <- setdiff(colnames(sp2), "family_id")
  
  template1 <- expand_grid(
    family_id = all_families,
    assembly = all_assemblies_sp1
  ) %>%
    mutate(species = species1_name)
  
  template2 <- expand_grid(
    family_id = all_families,
    assembly = all_assemblies_sp2
  ) %>%
    mutate(species = species2_name)
  
  df_full <- bind_rows(template1, template2) %>%
    left_join(bind_rows(df1, df2), by = c("family_id", "assembly", "species")) %>%
    mutate(copy = replace_na(copy, 0))
  
  list(
    df = bind_rows(df1, df2),
    df_full = df_full
  )
}

make_wide_matrix <- function(df_long) {
  wide <- df_long %>%
    pivot_wider(names_from = family_id, values_from = copy, values_fill = 0)
  
  mat <- wide %>%
    select(-assembly, -species) %>%
    as.matrix()
  
  keep_cols <- apply(mat, 2, var, na.rm = TRUE) > 0
  mat <- mat[, keep_cols, drop = FALSE]
  
  list(wide = wide, mat = mat)
}

make_pca_plot0 <- function(pca_obj, meta_df, clusters = NULL, max_overlap=20,title_text) {
  var_exp <- round(100 * summary(pca_obj)$importance[2, 1:2], 1)
  
  pca_df <- data.frame(pca_obj$x[, 1:2]) %>%
    mutate(
      assembly = meta_df$assembly,
      species = meta_df$species
    )
  
  # fallback if no cluster metadata is supplied
  if (is.null(clusters)) {
    return(
      ggplot(pca_df, aes(PC1, PC2, color = species)) +
        geom_point(size = 3) +
        geom_text_repel(
          aes(label = assembly),
          size = 3,
          max.overlaps = 20,
          box.padding = 0.5,
          point.padding = 0.3
        ) +
        theme_minimal() +
        labs(
          title = title_text,
          x = paste0("PC1 (", var_exp[1], "%)"),
          y = paste0("PC2 (", var_exp[2], "%)")
        )
    )
  }
  
  # exact match between assembly and cluster sample name
  cl_sub <- clusters[match(pca_df$assembly, clusters$Sample), , drop = FALSE]
  
  pca_df$geo_group <- cl_sub$Group
  pca_df$gen_cluster <- factor(cl_sub$LR_cluster)
  
  # fallback labels for unmatched rows
  pca_df$geo_group[is.na(pca_df$geo_group)] <- "Unknown"
  pca_df$gen_cluster[is.na(pca_df$gen_cluster)] <- "Unknown"
  pca_df$gen_cluster <- factor(pca_df$gen_cluster)
  
  ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(
      aes(
        color = geo_group,
        shape = gen_cluster
      ),
      size = 3,
      stroke = 1
    ) +
    scale_shape_manual(
      values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
      name = "Genetic Cluster"
    ) +
    scale_color_brewer(
      palette = "Set1",
      name = "Geographic Group"
    ) +
    geom_text_repel(
      aes(label = assembly),
      size = 3,
      color = "black",
      max.overlaps = max_overlap,
      box.padding = 0.5,
      point.padding = 0.3
    ) +
    theme_minimal() +
    labs(
      title = title_text,
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    guides(
      color = guide_legend(override.aes = list(shape = 15)),
      shape = guide_legend(override.aes = list(fill = "white"))
    )
}


make_pca_plot <- function(pca_obj, meta_df, clusters = NULL, max_overlap=20, title_text) {
  var_exp <- round(100 * summary(pca_obj)$importance[2, 1:2], 1)
  
  pca_df <- data.frame(pca_obj$x[, 1:2]) %>%
    mutate(
      assembly = meta_df$assembly,
      species = meta_df$species
    )
  
  # fallback if no cluster metadata is supplied
  if (is.null(clusters)) {
    return(
      ggplot(pca_df, aes(PC1, PC2, color = species)) +
        geom_point(size = 3) +
        geom_text_repel(
          aes(x = PC1, y = PC2, label = assembly),  # Explicit x and y
          size = 3,
          color = "black",
          max.overlaps = 20,
          box.padding = 0.5,
          point.padding = 0.3,
          inherit.aes = FALSE
        ) +
        theme_minimal() +
        labs(
          title = title_text,
          x = paste0("PC1 (", var_exp[1], "%)"),
          y = paste0("PC2 (", var_exp[2], "%)")
        )
    )
  }
  
  # exact match between assembly and cluster sample name
  cl_sub <- clusters[match(pca_df$assembly, clusters$Sample), , drop = FALSE]
  
  pca_df$geo_group <- cl_sub$Group
  pca_df$gen_cluster <- factor(cl_sub$LR_cluster)
  
  # fallback labels for unmatched rows
  pca_df$geo_group[is.na(pca_df$geo_group)] <- "Unknown"
  pca_df$gen_cluster[is.na(pca_df$gen_cluster)] <- "Unknown"
  pca_df$gen_cluster <- factor(pca_df$gen_cluster)
  
  ggplot(pca_df, aes(PC1, PC2)) +
    geom_point(
      aes(
        color = geo_group,
        shape = gen_cluster
      ),
      size = 3,
      stroke = 1
    ) +
    scale_shape_manual(
      values = c("1" = 21, "2" = 22, "3" = 23, "4" = 24),
      name = "Genetic Cluster"
    ) +
    scale_color_brewer(
      palette = "Set1",
      name = "Geographic Group"
    ) +
    geom_text_repel(
      aes(x = PC1, y = PC2, label = assembly),  
      size = 3,
      color = "black",
      max.overlaps = max_overlap,
      box.padding = 0.5,
      point.padding = 0.3,
      inherit.aes = FALSE
    ) +
    theme_minimal() +
    labs(
      title = title_text,
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    guides(
      color = guide_legend(override.aes = list(shape = 15)),
      shape = guide_legend(override.aes = list(fill = "white"))
    )
}

run_species_pca <- function(sp_df, species_name, clusters = NULL,max_overlap=10) {
  
  df_sp <- convert_long(sp_df, species_name)
  
  wide_sp <- df_sp %>%
    pivot_wider(names_from = family_id, values_from = copy, values_fill = 0)
  
  mat_sp <- wide_sp %>%
    select(-assembly, -species) %>%
    as.matrix()
  
  keep_cols <- apply(mat_sp, 2, var, na.rm = TRUE) > 0
  mat_sp <- mat_sp[, keep_cols, drop = FALSE]
  
  # Raw PCA
  pca_raw <- prcomp(log1p(mat_sp), scale. = TRUE)
  
  # Z-score PCA
  mat_z <- t(scale(t(log1p(mat_sp))))
  pca_z <- prcomp(mat_z, scale. = FALSE)
  
  # Family-centered PCA
  mat_fc <- scale(log1p(mat_sp), center = TRUE, scale = FALSE)
  pca_fc <- prcomp(mat_fc, scale. = FALSE)
  
  list(
    df = df_sp,
    wide = wide_sp,
    mat = mat_sp,
    mat_z = mat_z,
    mat_fc = mat_fc,
    pca_raw = pca_raw,
    pca_z = pca_z,
    pca_fc = pca_fc,
    p_raw = make_pca_plot(
      pca_obj = pca_raw,
      meta_df = wide_sp,
      clusters = clusters,
      max_overlap = max_overlap,
      title_text = paste0(species_name, " PCA - raw")
    ),
    p_z = make_pca_plot(
      pca_obj = pca_z,
      meta_df = wide_sp,
      clusters = clusters,
      max_overlap = max_overlap,
      title_text = paste0(species_name, " PCA - z-score normalized")
    ),
    p_fc = make_pca_plot(
      pca_obj = pca_fc,
      meta_df = wide_sp,
      clusters = clusters,
      max_overlap = max_overlap,
      title_text = paste0(species_name, " PCA - family-centered")
    )
  )
}

make_mirror_plot <- function(df, col_sp1, col_sp2, xlab, title_text,
                             species1_label, species2_label) {
  plot_df <- df %>%
    transmute(
      family_id = family_id,
      left_value = -1 * .data[[col_sp1]],
      right_value = .data[[col_sp2]]
    ) %>%
    pivot_longer(
      cols = c(left_value, right_value),
      names_to = "side",
      values_to = "value"
    ) %>%
    mutate(
      species = ifelse(side == "left_value", species1_label, species2_label)
    )
  
  max_abs <- max(abs(plot_df$value), na.rm = TRUE)
  
  ggplot(plot_df, aes(x = value, y = family_id, fill = species)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    scale_x_continuous(
      limits = c(-max_abs * 1.05, max_abs * 1.05),
      labels = function(x) abs(x)
    ) +
    labs(x = xlab, y = "Gene family", title = title_text) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 7)
    )
}
