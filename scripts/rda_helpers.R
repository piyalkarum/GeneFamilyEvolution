####### RDA helper functions ##########


make_formula <- function(response_name, explanatory, conditional = NULL) {
  rhs1 <- if (length(explanatory) > 0) paste(explanatory, collapse = " + ") else "1"
  rhs2 <- if (!is.null(conditional) && length(conditional) > 0) {
    paste0(" + Condition(", paste(conditional, collapse = " + "), ")")
  } else {
    ""
  }
  as.formula(paste0(response_name, " ~ ", rhs1, rhs2), env = parent.frame())
}

extract_rda_results <- function(mod, permutations = 999) {
  overall <- as.data.frame(anova.cca(mod, permutations = permutations))
  terms   <- as.data.frame(anova.cca(mod, by = "term", permutations = permutations))
  
  rr <- RsquareAdj(mod)
  
  if (is.list(rr)) {
    r2_val <- rr$r.squared
    adj_r2_val <- rr$adj.r.squared
  } else {
    rr_num <- as.numeric(rr)
    r2_val <- rr_num[1]
    adj_r2_val <- rr_num[2]
  }
  
  list(
    overall = overall,
    terms = terms,
    r2 = tibble(
      r2 = r2_val,
      adj_r2 = adj_r2_val
    )
  )
}
run_partial_rda <- function(Y, meta_df, explanatory, conditional, label) {
  Y <- as.matrix(Y)
  form <- make_formula("Y", explanatory, conditional)
  mod <- rda(form, data = meta_df)
  
  rr <- extract_rda_results(mod)
  
  site_scores <- scores(mod, display = "sites", choices = 1:2)
  site_df <- as.data.frame(site_scores)
  site_df$Population_name <- rownames(site_df)
  site_df <- left_join(
    site_df,
    meta_df %>% tibble::rownames_to_column("Population_name"),
    by = "Population_name"
  )
  
  bp <- tryCatch(scores(mod, display = "bp", choices = 1:2), error = function(e) NULL)
  bp_df <- NULL
  if (!is.null(bp)) {
    bp_df <- as.data.frame(bp)
    bp_df$variable <- rownames(bp_df)
  }
  
  list(
    label = label,
    model = mod,
    overall = rr$overall,
    terms = rr$terms,
    r2 = rr$r2,
    sites = site_df,
    biplot = bp_df
  )
}

plot_rda_ordination0 <- function(res_obj, pop_col = NULL, max_arrows = 8) {
  site_df <- res_obj$sites
  bp_df <- res_obj$biplot
  
  if (!is.null(bp_df) && nrow(bp_df) > 0) {
    bp_df <- bp_df %>%
      mutate(arrow_len = sqrt(RDA1^2 + RDA2^2)) %>%
      arrange(desc(arrow_len)) %>%
      slice(1:min(max_arrows, n()))
  }
  
  p <- ggplot(site_df, aes(RDA1, RDA2)) +
    geom_point(size = 1.7, color = "black", alpha = 0.65) +
    theme_classic(base_size = 11) +
    labs(
      title = res_obj$label,
      x = "RDA1",
      y = "RDA2"
    )
  
  if (!is.null(pop_col) && pop_col %in% names(site_df)) {
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = .data[[pop_col]]),
        size = 3,
        color = "black",
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.2,
        segment.color = "grey60",
        segment.size = 0.3
      )
  }
  
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
        hjust = 0.5, vjust = -0.3
      )
  }
  
  p
}


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



plot_rda_ordination1 <- function(res_obj, pop_col = NULL, clusters=NULL, max_arrows = 8) {
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


plot_rda_ordination <- function(res_obj, pop_col = NULL, clusters=NULL, 
                                max_arrows = 8, arrow_scale = 1, 
                                auto_scale = TRUE) {
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
    # Calculate scaling factor
    mult <- arrow_scale
    
    if (auto_scale) {
      # Scale arrows to be comparable to site scores range
      site_range <- max(abs(c(site_df$RDA1, site_df$RDA2)))
      arrow_max <- max(bp_df$arrow_len)
      if (arrow_max > 0) {
        mult <- site_range / arrow_max * 0.8  # 0.8 leaves some margin
      }
    }
    
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
        color = 2,
        hjust = 0.5, vjust = -0.3
      )
  }
  
  p
}

extract_varpart_summary <- function(vp_obj, label) {
  frac <- as.data.frame(vp_obj$part$fract)
  frac$fraction <- rownames(frac)
  rownames(frac) <- NULL
  frac$model <- label
  frac
}
