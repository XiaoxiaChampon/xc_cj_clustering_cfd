# ============================================================
# Hazel HPC — Scenario A | n=100, t=2000
# Cluster fractions:  75% / 22% / 3%
# Method:             catFDA-dbscan only
# Replicas:           100
# ============================================================

# setwd("D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd")
source("catfda_cluster_lib_hazel.R")

set.seed(123)

scenario     <- "A"
num_replicas <- 100
est_choice   <- "multinomial"

run_univfpca <- TRUE
run_kmeans   <- FALSE
run_fadp     <- FALSE
run_dbscan   <- TRUE
run_cfda     <- FALSE

temp_folder <- file.path("outputs", "clustersims",
                         paste0("A_", num_replicas, "_", est_choice, "_hazel_table2"))
if (!dir.exists(temp_folder)) dir.create(temp_folder, recursive = TRUE)
cat("Output folder:", temp_folder, "\n")

start_time <- Sys.time()

n100t2000 <- ClusterSimulation(
  100, 2000, scenario, num_replicas, est_choice,
  run_hellinger = FALSE, temp_folder,
  run_univfpca, run_kmeans, run_fadp, run_dbscan, run_cfda
)

end_time <- Sys.time()
print("run_hazel_A_n100t2000.R completed")
cat("\nTotal time:\n")
print(end_time - start_time)
flush(stdout())

cat("\n=== catFDA-dbscan ARI — Scenario A | n=100, t=2000 ===\n")
cat(sprintf("  Mean ARI : %.4f\n", n100t2000[[est_choice]]$cluster_table_est["dbscan ARI"]))
cat(sprintf("  SE (ARI) : %.4f\n", n100t2000[[est_choice]]$cluster_table_est_se["dbscan ARI"]))
cat(sprintf("  Replicas : %d\n", num_replicas))
flush(stdout())

if (run_parallel) {
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
