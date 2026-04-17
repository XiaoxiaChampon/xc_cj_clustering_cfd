# ============================================================
# Hazel HPC -- Scenario A | n=1000, t=2000
# Cluster fractions:  75% / 22% / 3%
# Method:             catFDA-dbscan only
# Replicas:           100
# ============================================================

# setwd("/path/to/submission_revision")
if (!exists("ClusterSimulation")) source("../catfda_cluster_lib.R")

batch_id_env <- Sys.getenv("LSB_JOBINDEX", unset = "")
if (nchar(batch_id_env) > 0) {
  batch_id <- as.integer(batch_id_env)
  set.seed(123 + batch_id * 1000)
} else {
  set.seed(123)
}

scenario     <- "A"
num_replicas <- 50
est_choice   <- "binomial"

run_univfpca <- TRUE
run_kmeans   <- FALSE
run_fadp     <- FALSE
run_dbscan   <- TRUE
run_cfda     <- FALSE

temp_folder <- if (nchar(batch_id_env) > 0) {
  file.path("outputs", "clustersims",
            paste0("A_", num_replicas, "_", est_choice, "_hazel_table1"),
            paste0("batch_", batch_id))
} else {
  file.path("outputs", "clustersims",
            paste0("A_", num_replicas, "_", est_choice, "_hazel_table1"))
}
if (!dir.exists(temp_folder)) dir.create(temp_folder, recursive = TRUE)
cat("Output folder:", temp_folder, "\n")

start_time <- Sys.time()

n1000t2000 <- ClusterSimulation(
  1000, 2000, scenario, num_replicas, est_choice,
  run_hellinger = TRUE, temp_folder,
  run_univfpca, run_kmeans, run_fadp, run_dbscan, run_cfda,
  save_curves = TRUE
)

end_time <- Sys.time()
print("run_hazel_table1_A_n1000t2000.R completed")
cat("\nTotal time:\n")
print(end_time - start_time)
flush(stdout())

cat("\n=== catFDA-dbscan ARI -- Scenario A | n=1000, t=2000 ===\n")
cat(sprintf("  Mean ARI : %.4f\n", n1000t2000[[est_choice]]$cluster_table_est["dbscan ARI"]))
cat(sprintf("  SE (ARI) : %.4f\n", n1000t2000[[est_choice]]$cluster_table_est_se["dbscan ARI"]))
cat(sprintf("  Replicas : %d\n", num_replicas))
flush(stdout())

if (run_parallel) {
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
