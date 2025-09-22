# Quick test for DAVID clustering using packaged example data
# Run from anywhere after installing/loading richCluster

suppressPackageStartupMessages(library(richCluster))

# ---- 1) Load example data from the package ----
path <- system.file("extdata", "cluster_result.rds", package = "richCluster")
if (identical(path, "")) {
  stop("Example data not found. Is 'richCluster' installed with extdata?")
}
cluster_result <- readRDS(path)

# Expecting list elements similar to what tests use
# cluster_result$df_list         # list of enrichment data.frames
# cluster_result$df_names        # optional names for df_list
# cluster_result$merged_df       # combined tidy frame (optional for this script)

# ---- 2) Build enrichment_results in the format david_* expect ----
# We will try to pick standard column names if present.
pick_cols <- function(df) {
  # Flexible selectors for common column names
  pick <- function(cands) {
    hit <- intersect(cands, names(df))
    if (length(hit)) hit[1] else NA_character_
  }
  term_col <- pick(c("Term","term","ID","id","Name","name"))
  gene_col <- pick(c("GeneID","gene_id","Genes","genes","GENE","Gene"))
  p_col    <- pick(c("Pvalue","pvalue","pval","P.Value","p"))
  padj_col <- pick(c("Padj","padj","FDR","qvalue","q"))
  
  # Minimal requirement: Term and GeneID
  if (is.na(term_col) || is.na(gene_col)) {
    stop("Could not find Term and GeneID-like columns in one of df_list entries.")
  }
  
  out <- data.frame(
    Term  = df[[term_col]],
    GeneID = df[[gene_col]],
    stringsAsFactors = FALSE
  )
  # Add Pvalue / Padj if available to keep downstream happy
  if (!is.na(p_col))    out$Pvalue <- df[[p_col]]
  if (!is.na(padj_col)) out$Padj   <- df[[padj_col]]
  out
}

stopifnot(is.list(cluster_result$df_list), length(cluster_result$df_list) >= 1)
enrichment_results <- lapply(cluster_result$df_list, pick_cols)
names(enrichment_results) <- cluster_result$df_names %||% paste0("df", seq_along(enrichment_results))

cat("Loaded example enrichment results from:", path, "\n")
cat("Datasets found:", paste(names(enrichment_results), collapse = ", "), "\n\n")

# ---- 3) Parameters used in your original snippet ----
similarity_threshold      <- 0.3
initial_group_membership  <- 2
final_group_membership    <- 2
multiple_linkage_threshold<- 0.5

# ---- 4) Helper to print cluster summaries ----
print_clusters <- function(res, label) {
  if (!is.list(res) || is.null(res$clusters)) {
    cat("No 'clusters' element returned by", label, "\n")
    return(invisible(NULL))
  }
  cl <- res$clusters
  n  <- tryCatch(nrow(cl), error = function(e) NA_integer_)
  cat("Number of clusters found by", label, ":", n, "\n")
  if (is.data.frame(cl) && n > 0) {
    print(utils::head(cl, 10))
  }
  invisible(NULL)
}

# ---- 5) Run DAVID implementations on the example data ----
cat("Testing DAVID clustering implementations on packaged example data...\n\n")

## C++ implementation
cat("=== C++ implementation (david_cluster) ===\n")
res_cpp <- tryCatch({
  david_cluster(
    enrichment_results = enrichment_results,
    similarity_threshold = similarity_threshold,
    initial_group_membership = initial_group_membership,
    final_group_membership = final_group_membership,
    multiple_linkage_threshold = multiple_linkage_threshold
  )
}, error = function(e) {
  cat("C++ implementation failed:", conditionMessage(e), "\n")
  NULL
})
if (!is.null(res_cpp)) {
  cat("✓ C++ implementation succeeded\n")
  print_clusters(res_cpp, "C++")
}
cat("\n")
