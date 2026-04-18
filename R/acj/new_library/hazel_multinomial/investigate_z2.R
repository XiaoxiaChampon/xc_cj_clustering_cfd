# investigate_z2_s3_t300.R
# Investigates the anomalous z2 MSE for cluster 3 (s3), t=300

data_folder <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/A_100_multinomial_hazel_table1"
raw_folder  <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/RAW_HAZEL/A_10_multinomial_hazel_table1"

e <- new.env()
load(file.path(data_folder, "ClusterSim_1000_300_A_100_multinomial_TRUE_neworder.RData"), envir = e)
obj <- e$est_values

cat("=== est_values names ===\n")
print(names(obj))
cat("\n=== Full MSE matrix ===\n")
print(obj$mse)

# Per-replica MSE if available
cat("\n=== mse_reps structure ===\n")
if (!is.null(obj$mse_reps)) {
  cat("dim:", dim(obj$mse_reps), "\n")  # expect 100 x 2 x 3 or similar
  print(str(obj$mse_reps))
}
if (!is.null(obj$rmse_reps)) {
  cat("rmse_reps dim:", paste(dim(obj$rmse_reps), collapse="x"), "\n")
  # z2 is row 2, cluster 3 is col 3
  z2_s3 <- obj$rmse_reps[, 2, 3]^2   # square RMSE to get MSE
  cat("z2 s3 per-replica MSE (first 20):\n")
  print(round(head(z2_s3, 20), 3))
  cat("Summary:\n")
  print(summary(z2_s3))
  
  # Check which batch each replica came from
  outlier_reps <- which(z2_s3 > 10)
  cat("\nOutlier replicas (z2_s3 MSE > 10):", outlier_reps, "\n")
  cat("Their MSE values:", round(z2_s3[outlier_reps], 2), "\n")
}

# Also check hellinger_reps
cat("\n=== hellinger_reps ===\n")
if (!is.null(obj$hellinger_reps)) {
  cat("dim:", paste(dim(obj$hellinger_reps), collapse="x"), "\n")
}
