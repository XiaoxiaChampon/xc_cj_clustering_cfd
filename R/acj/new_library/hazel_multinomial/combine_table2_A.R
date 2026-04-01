# ============================================================
# combine_table2_A.R
# Prepare Table 2 inputs for Scenario A | n=100 and n=500, t=2000
#
# Strategy (in priority order):
#   1. If a complete 100-replica file already exists in
#      outputs/PaperTables/Table2Files/A_100_multinomial_hazel_table2/
#      -> copy it to the standard output folder and use it directly.
#   2. If batch outputs exist in
#      outputs/clustersims/A_10_multinomial_hazel_table2/batch_1..10/
#      -> pool 10 batches x 10 replicas = 100 total (using per-replica vectors).
#   3. Otherwise -> stop with an informative error.
#
# Note: n=1000, t=2000 for Scenario A is handled by combine_table1_A.R.
#
# Output folder (both strategies write here):
#   outputs/clustersims/A_100_multinomial_hazel_table2/
#     ClusterSim_{n}_2000_A_100_multinomial_FALSE_neworder.RData
# ============================================================

num_batches            <- 10
num_replicas_per_batch <- 10
scenario               <- "A"
tt                     <- 2000
est_choice             <- "multinomial"
n_sizes                <- c(100, 500)   # n=1000 handled by combine_table1_A.R

# Paths
papertables_src <- file.path(
  "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd",
  "outputs", "PaperTables", "Table2Files",
  "A_100_multinomial_hazel_table2"
)
batch_root <- file.path("outputs", "clustersims",
                        "A_10_multinomial_hazel_table2")
out_folder <- file.path("outputs", "clustersims",
                        "A_100_multinomial_hazel_table2")
if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

suffix_existing <- "FALSE_neworder.RData"
N_total         <- num_batches * num_replicas_per_batch   # 100

for (nn in n_sizes) {
  cat("\n=== Scenario A | n=", nn, "t=", tt, "===\n")

  out_fname <- paste("ClusterSim", nn, tt, scenario,
                     N_total, est_choice, suffix_existing, sep = "_")
  out_path  <- file.path(out_folder, out_fname)

  # ---------- Strategy 1: existing complete file ----------
  existing_path <- file.path(papertables_src, out_fname)

  if (file.exists(existing_path)) {
    cat("  Found existing 100-replica file. Copying to output folder.\n")
    cat("  Source:", existing_path, "\n")
    file.copy(existing_path, out_path, overwrite = TRUE)
    e <- new.env(); load(out_path, envir = e); ev <- e$est_values
    cat("  dbscan ARI :", round(ev$cluster_table_est["dbscan ARI"], 4),
        "  SE:", round(ev$cluster_table_est_se["dbscan ARI"], 4), "\n")
    cat("  dbscan RI  :", round(ev$cluster_table_est["dbscan RI"],  4),
        "  SE:", round(ev$cluster_table_est_se["dbscan RI"],  4), "\n")
    flush(stdout())
    next
  }

  # ---------- Strategy 2: pool batches ----------
  batch_fname      <- paste("ClusterSim", nn, tt, scenario,
                             num_replicas_per_batch, est_choice, suffix_existing, sep = "_")
  first_batch_path <- file.path(batch_root, "batch_1", batch_fname)
  if (!file.exists(first_batch_path)) {
    stop("Neither complete file nor batch outputs found for n=", nn,
         ".\n  Expected complete: ", existing_path,
         "\n  Expected batch:    ", first_batch_path)
  }

  cat("  Pooling", num_batches, "batches x", num_replicas_per_batch, "replicas\n")
  ari_all <- numeric(0)
  ri_all  <- numeric(0)
  last_ev <- NULL

  for (b in seq_len(num_batches)) {
    fpath <- file.path(batch_root, paste0("batch_", b), batch_fname)
    cat("  batch", b, ":", batch_fname, "\n")
    e <- new.env(); load(fpath, envir = e); ev <- e$est_values

    ari_all <- c(ari_all, ev$est_dbscan_ari_reps)
    ri_all  <- c(ri_all,  ev$est_dbscan_ri_reps)
    last_ev <- ev
  }

  stopifnot(length(ari_all) == N_total, length(ri_all) == N_total)

  combined_est    <- last_ev$cluster_table_est
  combined_est_se <- last_ev$cluster_table_est_se
  combined_est["dbscan ARI"]    <- mean(ari_all)
  combined_est["dbscan RI"]     <- mean(ri_all)
  combined_est_se["dbscan ARI"] <- sd(ari_all) / sqrt(N_total)
  combined_est_se["dbscan RI"]  <- sd(ri_all)  / sqrt(N_total)

  est_values <- list(
    cluster_table_est    = combined_est,
    cluster_table_est_se = combined_est_se
  )
  save(est_values, file = out_path)
  cat("  Saved:", out_fname, "\n")
  cat("  dbscan ARI :", round(combined_est["dbscan ARI"], 4),
      "  SE:", round(combined_est_se["dbscan ARI"], 4), "\n")
  cat("  dbscan RI  :", round(combined_est["dbscan RI"],  4),
      "  SE:", round(combined_est_se["dbscan RI"],  4), "\n")
  flush(stdout())
}

cat("\nDone. Output folder:", out_folder, "\n")
cat("n=1000 t=2000 for Scenario A -> run combine_table1_A.R instead.\n")