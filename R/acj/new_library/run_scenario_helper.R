# ============================================================
# Shared helper for scenario runner scripts.
# Requires ClusterSimulation() to already be loaded.
#
# Usage:
#   results <- RunScenarioSims(scenario, n_values, t_values,
#                              num_replicas, est_choice, temp_folder)
#   results$true_table   — RI/ARI on true cluster labels
#   results$est_table    — RI/ARI on estimated cluster labels
#   results$est_table_se — standard errors of above
# ============================================================

RunScenarioSims <- function(scenario, n_values, t_values,
                             num_replicas, est_choice, temp_folder,
                             run_univfpca = TRUE, run_kmeans = TRUE,
                             run_fadp = TRUE, run_dbscan = TRUE,
                             run_cfda = TRUE) {

  results  <- list()

  for (n in n_values) {
    for (t in t_values) {
      key <- paste0("n", n, "t", t)
      cat("\n--- Running:", key, "---\n")
      results[[key]] <- ClusterSimulation(
        n, t, scenario, num_replicas, est_choice, TRUE, temp_folder,
        run_univfpca, run_kmeans, run_fadp, run_dbscan, run_cfda
      )
    }
  }

  # Compile tables from the ordered result list
  true_table   <- do.call(rbind, lapply(results, `[[`, "cluster_table_true"))
  est_table    <- do.call(rbind, lapply(results, `[[`, "cluster_table_est"))
  est_table_se <- do.call(rbind, lapply(results, `[[`, "cluster_table_est_se"))
  rownames(true_table) <- rownames(est_table) <- rownames(est_table_se) <- names(results)

  # Print RI / ARI columns
  ri_cols <- grep("RI|ARI", colnames(true_table))
  cat("\n===== Scenario", scenario, ": TRUE labels — RI / ARI =====\n")
  print(true_table[, ri_cols])
  cat("\n===== Scenario", scenario, ": ESTIMATED labels — RI / ARI =====\n")
  print(est_table[, ri_cols])
  cat("\n===== Scenario", scenario, ": ESTIMATED labels — RI / ARI (SE) =====\n")
  print(est_table_se[, ri_cols])

  # Save
  out_file <- file.path(temp_folder, paste0("scenario", scenario, "_results.RData"))
  save(true_table, est_table, est_table_se, file = out_file)
  cat("\nResults saved to", out_file, "\n")

  invisible(list(
    true_table   = true_table,
    est_table    = est_table,
    est_table_se = est_table_se
  ))
}
