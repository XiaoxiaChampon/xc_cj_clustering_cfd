data_folder <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/A_100_binomial_hazel_table1"
e <- new.env()
load(file.path(data_folder, "ClusterSim_1000_300_A_100_binomial_TRUE_neworder.RData"), envir = e)
obj <- e[["est_values"]]

cat("mse matrix:\n")
print(obj[["mse"]])
cat("\ndim rmse_reps:", dim(obj[["rmse_reps"]]), "\n")

z2_s3 <- obj[["rmse_reps"]][, 2, 3]^2
cat("\nz2 s3 - mean:", mean(z2_s3), "\n")
cat("z2 s3 - summary:\n")
print(summary(z2_s3))
cat("\noutlier replicas (>5):\n")
print(which(z2_s3 > 5))
cat("\ntop 10 values:\n")
print(sort(z2_s3, decreasing = TRUE)[1:10])
