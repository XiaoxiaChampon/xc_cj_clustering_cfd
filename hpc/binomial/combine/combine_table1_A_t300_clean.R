# ============================================================
# combine_table1_A_t300_clean.R
#
# Re-combine t=300 batches for Table 1 — Scenario A | n=1000
# with outlier-replica exclusion.
#
# Strategy:
#   Load batch_0 (original 100 reps) plus any new reruns
#   (batch_1..batch_5).  Exclude any replica where the
#   per-replica RMSE for Z2 in any cluster exceeds OUTLIER_THRESH.
#   Take the first TARGET_N clean reps and recompute the combined
#   ClusterSim RData file.
#
# Outlier evidence: batch_0 reps 42 & 70 have z2/cluster-3 RMSE ~251.
# OUTLIER_THRESH = 20 is a conservative cut (next highest is ~11).
#
# Output overwrites:
#   outputs/clustersims/A_100_binomial_hazel_table1/
#     ClusterSim_1000_300_A_100_binomial_TRUE_neworder.RData
# ============================================================

OUTLIER_THRESH <- 20   # RMSE threshold for z2 outlier detection
TARGET_N       <- 100  # desired clean replicas in final combined file

scenario      <- "A"
num_indvs     <- 1000
est_choice    <- "binomial"
suffix        <- "TRUE_neworder.RData"

batch_root  <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/RAW_HAZEL_BINOMIAL/A_100_binomial_hazel_table1"
out_folder  <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/A_100_binomial_hazel_table1"

fname <- paste("ClusterSim", num_indvs, 300, scenario, 100, est_choice, suffix, sep = "_")

# Batch IDs to try: 0 (original) then 1-5 (reruns)
batch_ids <- c(0, 1:5)

# ---- Accumulate clean replicas ----------------------------------------
ari_clean       <- numeric(0)
ri_clean        <- numeric(0)
rmse_clean      <- matrix(0, nrow = 0, ncol = 6)   # will rbind 2x3 rows
hellinger_clean <- matrix(0, nrow = 0, ncol = 9)   # will rbind 3x3 rows
last_ev         <- NULL
total_loaded    <- 0
total_excluded  <- 0

for (b in batch_ids) {
  fpath <- file.path(batch_root, paste0("batch_", b), fname)
  if (!file.exists(fpath)) {
    cat("  batch", b, ": not found, skipping\n")
    next
  }
  e <- new.env()
  load(fpath, envir = e)
  ev <- e[["est_values"]]

  rr  <- ev[["rmse_reps"]]       # [reps x 2 x 3]
  hr  <- ev[["hellinger_reps"]]  # [reps x 3 x 3]
  ari <- ev[["est_dbscan_ari_reps"]]
  ri  <- ev[["est_dbscan_ri_reps"]]
  n_reps <- dim(rr)[1]

  # Outlier flag: any z2 (row 2) RMSE > threshold
  z2_max <- apply(rr[, 2, , drop = FALSE], 1, max)  # max z2 across 3 clusters
  keep   <- z2_max <= OUTLIER_THRESH
  n_keep <- sum(keep)
  n_drop <- sum(!keep)

  cat(sprintf("  batch %d : %d reps loaded, %d outliers excluded (z2 RMSE > %.0f), %d kept\n",
              b, n_reps, n_drop, OUTLIER_THRESH, n_keep))
  if (n_drop > 0) {
    cat("    excluded replica indices:", which(!keep), "\n")
    cat("    their z2 max values:", round(z2_max[!keep], 2), "\n")
  }

  # Append clean reps
  ari_clean       <- c(ari_clean, ari[keep])
  ri_clean        <- c(ri_clean,  ri[keep])
  # Flatten each rep's 2x3 / 3x3 into a row for storage
  for (i in which(keep)) {
    rmse_clean      <- rbind(rmse_clean,      as.vector(rr[i, , ]))
    hellinger_clean <- rbind(hellinger_clean, as.vector(hr[i, , ]))
  }

  total_loaded   <- total_loaded + n_reps
  total_excluded <- total_excluded + n_drop
  last_ev <- ev

  if (length(ari_clean) >= TARGET_N) {
    cat(sprintf("  Reached %d clean reps — stopping batch loading.\n", TARGET_N))
    break
  }
}

N_clean <- length(ari_clean)
cat(sprintf("\nTotal loaded: %d reps across batches, excluded: %d, clean available: %d\n",
            total_loaded, total_excluded, N_clean))

if (N_clean < TARGET_N) {
  warning(sprintf("Only %d clean reps available (need %d). Run more array batches.", N_clean, TARGET_N))
}

# Trim to TARGET_N if we have more
use_n <- min(N_clean, TARGET_N)
ari_use       <- ari_clean[1:use_n]
ri_use        <- ri_clean[1:use_n]
rmse_use      <- rmse_clean[1:use_n, , drop = FALSE]
hellinger_use <- hellinger_clean[1:use_n, , drop = FALSE]

# Reconstruct [use_n x 2 x 3] and [use_n x 3 x 3] arrays for mean computation
rmse_arr      <- array(rmse_use,      dim = c(use_n, 2, 3))
hellinger_arr <- array(hellinger_use, dim = c(use_n, 3, 3))

# ---- Build combined est_values ----------------------------------------
combined_est    <- last_ev[["cluster_table_est"]]
combined_est_se <- last_ev[["cluster_table_est_se"]]

combined_est["dbscan ARI"]    <- mean(ari_use)
combined_est["dbscan RI"]     <- mean(ri_use)
combined_est_se["dbscan ARI"] <- sd(ari_use) / sqrt(use_n)
combined_est_se["dbscan RI"]  <- sd(ri_use)  / sqrt(use_n)

est_values <- list(
  cluster_table_est    = combined_est,
  cluster_table_est_se = combined_est_se,
  mse                  = apply(rmse_arr,      c(2, 3), mean),
  hellinger            = apply(hellinger_arr, c(2, 3), mean)
)

cat(sprintf("\nFinal: %d clean reps used\n", use_n))
cat("MSE matrix:\n");       print(round(est_values[["mse"]],       4))
cat("Hellinger matrix:\n"); print(round(est_values[["hellinger"]], 4))
cat(sprintf("dbscan ARI: %.4f (SE: %.4f)\n",
            est_values[["cluster_table_est"]]["dbscan ARI"],
            est_values[["cluster_table_est_se"]]["dbscan ARI"]))

# ---- Save ---------------------------------------------------------------
out_fname <- paste("ClusterSim", num_indvs, 300, scenario,
                   use_n, est_choice, suffix, sep = "_")
save(est_values, file = file.path(out_folder, out_fname))
cat("\nSaved:", file.path(out_folder, out_fname), "\n")
