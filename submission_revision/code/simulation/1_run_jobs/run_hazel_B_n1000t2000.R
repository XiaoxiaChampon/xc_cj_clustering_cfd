# ============================================================
# Hazel HPC — Scenario B | n=1000, t=2000
# Cluster fractions:  50% / 30% / 20%
# mu_2 (setting 1):   -0.5 + exp(2t)  [closer means, harder separation]
# Method:             catFDA-dbscan only
# Replicas:            10
# ============================================================

# setwd("/path/to/submission_revision")
if (!exists("ClusterSimulation")) source("../catfda_cluster_lib.R")

set.seed(123)

scenario     <- "B"
num_replicas <- 100
est_choice   <- "binomial"

run_univfpca <- TRUE
run_kmeans   <- FALSE
run_fadp     <- FALSE
run_dbscan   <- TRUE
run_cfda     <- FALSE

temp_folder <- file.path("outputs", "clustersims",
                         paste0("B_", num_replicas, "_", est_choice, "_hazel_table2"))
if (!dir.exists(temp_folder)) dir.create(temp_folder, recursive = TRUE)
cat("Output folder:", temp_folder, "\n")

start_time <- Sys.time()

n1000t2000 <- ClusterSimulation(
  1000, 2000, scenario, num_replicas, est_choice,
  run_hellinger = FALSE, temp_folder,
  run_univfpca, run_kmeans, run_fadp, run_dbscan, run_cfda,
  save_curves = FALSE
)

end_time <- Sys.time()
print("run_hazel_table2_B_n1000t2000.R completed")
cat("\nTotal time:\n")
print(end_time - start_time)
flush(stdout())

cat("\n=== catFDA-dbscan ARI — Scenario B | n=1000, t=2000 ===\n")
cat(sprintf("  Mean ARI : %.4f\n", n1000t2000[[est_choice]]$cluster_table_est["dbscan ARI"]))
cat(sprintf("  SE (ARI) : %.4f\n", n1000t2000[[est_choice]]$cluster_table_est_se["dbscan ARI"]))
cat(sprintf("  Replicas : %d\n", num_replicas))
flush(stdout())

if (run_parallel) {
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
