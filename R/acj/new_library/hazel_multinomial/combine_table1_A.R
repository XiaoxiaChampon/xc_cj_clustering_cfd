# ============================================================
# combine_table1_A.R
# Combine batch outputs for Table 1 — Scenario A | n=1000, t=300/750/2000
#
# Batch configs (reps_per_batch x num_batches = 100 reps total):
#   t=300 : 10 reps/batch x 10 batches -> A_10_multinomial_hazel_table1/
#   t=750 : 5 reps/batch  x 20 batches -> A_5_multinomial_hazel_table1/
#   t=2000: 5 reps/batch  x 20 batches -> A_5_multinomial_hazel_table1/
#
# This script pools all 100 replicas per t-length and computes:
#   - cluster_table_est / _se  : mean and SE of dbscan ARI and RI over 100 reps
#   - mse, hellinger           : grand-mean 2x3 and 3x3 matrices
#
# Output folder:
#   outputs/clustersims/A_100_multinomial_hazel_table1/
#     ClusterSim_1000_{t}_A_100_multinomial_TRUE_neworder.RData
# ============================================================

scenario      <- "A"
num_indvs     <- 1000
est_choice    <- "multinomial"
t_lengths     <- c(300, 750, 2000)
run_hellinger <- TRUE

# Per-t batch configuration: reps_per_batch x num_batches = 100 total
t_configs <- list(
  "300"  = list(num_batches = 10, reps_per_batch = 10, folder = "A_10_multinomial_hazel_table1"),
  "750"  = list(num_batches = 20, reps_per_batch =  5, folder = "A_5_multinomial_hazel_table1"),
  "2000" = list(num_batches = 20, reps_per_batch =  5, folder = "A_5_multinomial_hazel_table1")
)

batch_root_base <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/Batched_Hazel_Table1"

out_folder <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/A_100_multinomial_hazel_table1"
if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

suffix  <- "TRUE_neworder.RData"
N_total <- 100   # always 100 regardless of batch config

for (tt in t_lengths) {
  cfg                    <- t_configs[[as.character(tt)]]
  num_batches            <- cfg$num_batches
  num_replicas_per_batch <- cfg$reps_per_batch
  batch_root             <- file.path(batch_root_base, cfg$folder)

  cat("\n=== Combining batches: Scenario A | n=", num_indvs, "t=", tt,
      "(", num_batches, "batches x", num_replicas_per_batch, "reps) ===\n")

  ari_all <- numeric(0)
  ri_all  <- numeric(0)
  mse_sum        <- matrix(0, nrow = 2, ncol = 3)
  hellinger_sum  <- matrix(0, nrow = 3, ncol = 3)
  batches_used   <- 0
  last_ev        <- NULL

  for (b in seq_len(num_batches)) {
    fname <- paste("ClusterSim", num_indvs, tt, scenario,
                   num_replicas_per_batch, est_choice, suffix, sep = "_")
    fpath <- file.path(batch_root, paste0("batch_", b), fname)
    if (!file.exists(fpath)) {
      warning("  MISSING batch ", b, " for t=", tt, " -- skipping: ", fpath)
      next
    }
    cat("  batch", b, ":", fname, "\n")
    e <- new.env()
    load(fpath, envir = e)
    ev <- e$est_values

    ari_all <- c(ari_all, ev$est_dbscan_ari_reps)
    ri_all  <- c(ri_all,  ev$est_dbscan_ri_reps)

    # Average rmse/hellinger per batch (each is [10,2,3] or [10,3,3])
    mse_sum       <- mse_sum       + apply(ev$rmse_reps,      c(2, 3), mean)
    hellinger_sum <- hellinger_sum + apply(ev$hellinger_reps, c(2, 3), mean)

    batches_used <- batches_used + 1
    last_ev <- ev
  }

  stopifnot(length(ari_all) == length(ri_all))
  N_actual <- length(ari_all)
  cat("  Batches used:", batches_used, "/ ", num_batches,
      "-- Total replicas:", N_actual, "\n")
  if (N_actual < num_batches * num_replicas_per_batch)
    warning("  Only ", N_actual, " of ", num_batches * num_replicas_per_batch,
            " expected replicas for t=", tt)

  # Build combined est_values
  combined_est    <- last_ev$cluster_table_est
  combined_est_se <- last_ev$cluster_table_est_se

  combined_est["dbscan ARI"]    <- mean(ari_all)
  combined_est["dbscan RI"]     <- mean(ri_all)
  combined_est_se["dbscan ARI"] <- sd(ari_all) / sqrt(N_actual)
  combined_est_se["dbscan RI"]  <- sd(ri_all)  / sqrt(N_actual)

  est_values <- list(
    cluster_table_est    = combined_est,
    cluster_table_est_se = combined_est_se,
    mse                  = mse_sum       / batches_used,
    hellinger            = hellinger_sum / batches_used
  )

  out_fname <- paste("ClusterSim", num_indvs, tt, scenario,
                     N_actual, est_choice, suffix, sep = "_")
  save(est_values, file = file.path(out_folder, out_fname))
  cat("  Saved:", out_fname, "\n")

  cat("  dbscan ARI :", round(est_values$cluster_table_est["dbscan ARI"], 4),
      "  SE:", round(est_values$cluster_table_est_se["dbscan ARI"], 4), "\n")
  cat("  MSE:\n");       print(round(est_values$mse,       4))
  cat("  Hellinger:\n"); print(round(est_values$hellinger, 4))
  flush(stdout())
}

cat("\nDone. Combined files in:", out_folder, "\n")
