# ============================================================
# Compare RI / ARI clustering results: old code vs new code
# Loads existing .RData outputs and prints tables side by side.
# ============================================================
# setwd("D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd")

scenario     <- "A"
num_replicas <- 2
est_choice   <- "multinomial"
t_lengths    <- c(300, 750, 2000)
n_values     <- c(100)           # adjust if more n were run

old_folder <- file.path("outputs", "clustersims",
                        paste(scenario, num_replicas, est_choice,
                              "to_compare_old_code", sep = "_"))
new_folder <- file.path("outputs", "clustersims",
                        paste(scenario, num_replicas, est_choice,
                              "to_compare_new_code", sep = "_"))

# Helper: load est_values for a given n and t
load_est <- function(folder, n, t, suffix) {
  fname <- paste0("ClusterSim_", n, "_", t, "_", scenario, "_",
                  num_replicas, "_", est_choice, "_TRUE", suffix, ".RData")
  path  <- file.path(folder, fname)
  if (!file.exists(path)) stop("File not found: ", path)
  e <- new.env()
  load(path, envir = e)
  e$est_values
}

ri_ari_cols <- c("cfda-db RI", "cfda-db ARI",
                 "fadp RI",    "fadp ARI",
                 "kmeans RI",  "kmeans ARI",
                 "dbscan RI",  "dbscan ARI")

for (n in n_values) {
  cat("\n", strrep("=", 70), "\n")
  cat(" n =", n, "— Estimated cluster labels: RI and ARI\n")
  cat(strrep("=", 70), "\n")

  old_rows <- lapply(t_lengths, function(t) {
    ev <- load_est(old_folder, n, t, "")
    ev$cluster_table_est[ri_ari_cols]
  })
  new_rows <- lapply(t_lengths, function(t) {
    ev <- load_est(new_folder, n, t, "_neworder")
    ev$cluster_table_est[ri_ari_cols]
  })

  old_mat <- do.call(rbind, old_rows)
  new_mat <- do.call(rbind, new_rows)
  rownames(old_mat) <- rownames(new_mat) <- paste0("n", n, "t", t_lengths)

  cat("\n--- OLD code ---\n");  print(round(old_mat, 4))
  cat("\n--- NEW code ---\n");  print(round(new_mat, 4))

  diff_mat <- new_mat - old_mat
  cat("\n--- Difference (new - old) ---\n");  print(round(diff_mat, 4))
}
