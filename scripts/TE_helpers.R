
# =========================
# TE ANNALYSIS HELPERS
# =========================

read_repeatmasker <- function(f) {
  tm <- fread(f, skip = 3, fill = Inf)
  tm <- as.data.frame(tm)
  
  # keep only first 15 standard RepeatMasker columns
  tm <- tm[, 1:15]
  
  colnames(tm) <- c(
    "SW_score", "perc_div", "perc_del", "perc_ins",
    "query_sequence", "query_start", "query_end", "query_left",
    "strand", "repeat_name", "repeat_class",
    "repeat_start", "repeat_end", "repeat_left", "ID"
  )
  
  tm %>%
    filter(repeat_class %in% te_classes) %>%
    transmute(
      scaffold = query_sequence,
      te_start = as.integer(query_start),
      te_end = as.integer(query_end),
      repeat_name,
      repeat_class
    )
}

make_gene_windows <- function(df, buffer = 0) {
  df %>%
    mutate(
      win_start = pmax(start - buffer, 1),
      win_end   = end + buffer,
      window_bp = win_end - win_start + 1,
      window = ifelse(buffer == 0, "gene_body", paste0("pm", buffer/1000, "kb"))
    )
}

to_granges_genes <- function(df) {
  GRanges(
    seqnames = df$scaffold,
    ranges = IRanges(df$win_start, df$win_end),
    gene = df$gene,
    assembly = df$assembly,
    group = df$group,
    window = df$window,
    window_bp = df$window_bp
  )
}

to_granges_tes <- function(df) {
  GRanges(
    seqnames = df$scaffold,
    ranges = IRanges(df$te_start, df$te_end),
    repeat_class = df$repeat_class,
    repeat_name = df$repeat_name
  )
}

# calculate per-gene TE metrics
summarize_te_metrics <- function(gene_df, te_df, buffer) {
  
  gene_win <- make_gene_windows(gene_df, buffer)
  gr_gene <- to_granges_genes(gene_win)
  gr_te   <- to_granges_tes(te_df)
  
  # nearest TE distance
  nearest_hits <- distanceToNearest(gr_gene, gr_te, ignore.strand = TRUE)
  nearest_df <- data.frame(
    queryHits = queryHits(nearest_hits),
    nearest_te_dist = mcols(nearest_hits)$distance
  )
  
  # overlaps
  ov <- findOverlaps(gr_gene, gr_te, ignore.strand = TRUE)
  ov_df <- data.frame(
    gene_idx = queryHits(ov),
    te_idx = subjectHits(ov)
  )
  
  # if no overlaps
  if (nrow(ov_df) == 0) {
    out <- gene_win %>%
      mutate(
        n_te = 0,
        te_bp = 0,
        te_fraction = 0,
        te_classes = NA_character_,
        nearest_te_dist = NA_real_
      )
    
    if (nrow(nearest_df) > 0) {
      out$nearest_te_dist[nearest_df$queryHits] <- nearest_df$nearest_te_dist
    }
    
    return(out)
  }
  
  ov_join <- ov_df %>%
    mutate(
      gene = mcols(gr_gene)$gene[gene_idx],
      assembly = mcols(gr_gene)$assembly[gene_idx],
      group = mcols(gr_gene)$group[gene_idx],
      window = mcols(gr_gene)$window[gene_idx],
      window_bp = mcols(gr_gene)$window_bp[gene_idx],
      g_start = start(gr_gene)[gene_idx],
      g_end   = end(gr_gene)[gene_idx],
      te_start = start(gr_te)[te_idx],
      te_end   = end(gr_te)[te_idx],
      repeat_class = mcols(gr_te)$repeat_class[te_idx],
      repeat_name  = mcols(gr_te)$repeat_name[te_idx]
    ) %>%
    mutate(
      ov_start = pmax(g_start, te_start),
      ov_end   = pmin(g_end, te_end),
      ov_bp    = pmax(0, ov_end - ov_start + 1)
    )
  
  te_bp_df <- ov_join %>%
    group_by(gene, assembly, group, window, window_bp) %>%
    summarise(
      n_te = n(),
      te_classes = paste(sort(unique(repeat_class)), collapse = ", "),
      .groups = "drop"
    )
  
  te_bp_reduce <- ov_join %>%
    group_by(gene, assembly, group, window, window_bp) %>%
    group_split() %>%
    purrr::map_dfr(function(x) {
      ir <- IRanges::IRanges(start = x$ov_start, end = x$ov_end)
      tibble(
        gene = x$gene[1],
        assembly = x$assembly[1],
        group = x$group[1],
        window = x$window[1],
        window_bp = x$window_bp[1],
        te_bp = sum(IRanges::width(IRanges::reduce(ir)))
      )
    })
  
  out <- gene_win %>%
    left_join(te_bp_df, by = c("gene", "assembly", "group", "window", "window_bp")) %>%
    left_join(te_bp_reduce, by = c("gene", "assembly", "group", "window", "window_bp")) %>%
    mutate(
      n_te = dplyr::coalesce(n_te, 0L),
      te_bp = dplyr::coalesce(te_bp, 0),
      te_fraction = te_bp / window_bp
    )
  
  out$nearest_te_dist <- NA_real_
  if (nrow(nearest_df) > 0) {
    out$nearest_te_dist[nearest_df$queryHits] <- nearest_df$nearest_te_dist
  }
  
  out
}
# class-specific counts per gene/window
summarize_te_classes <- function(gene_df, te_df, buffer) {
  gene_win <- make_gene_windows(gene_df, buffer)
  gr_gene <- to_granges_genes(gene_win)
  gr_te   <- to_granges_tes(te_df)
  
  ov <- findOverlaps(gr_gene, gr_te, ignore.strand = TRUE)
  if (length(ov) == 0) return(NULL)
  
  tibble(
    gene = mcols(gr_gene)$gene[queryHits(ov)],
    assembly = mcols(gr_gene)$assembly[queryHits(ov)],
    group = mcols(gr_gene)$group[queryHits(ov)],
    window = mcols(gr_gene)$window[queryHits(ov)],
    repeat_class = mcols(gr_te)$repeat_class[subjectHits(ov)]
  ) %>%
    count(gene, assembly, group, window, repeat_class, name = "n_te_class")
}
