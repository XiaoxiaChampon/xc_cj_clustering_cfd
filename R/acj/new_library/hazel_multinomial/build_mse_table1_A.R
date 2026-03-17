# Build MSE table for Table 1 results from Hazel HPC.
#
# Reads three ClusterSimulation outputs from:
#   outputs/clustersims/A_100_multinomial_hazel_table1/
# (produced by combine_table1_A.R after all 10 array-job batches complete)
#
# Files expected:
#   ClusterSim_1000_300_A_100_multinomial_TRUE_neworder.RData
#   ClusterSim_1000_750_A_100_multinomial_TRUE_neworder.RData
#   ClusterSim_1000_2000_A_100_multinomial_TRUE_neworder.RData
#
# Each .RData file contains a list named `est_values` with:
#   $mse       : 2×3 matrix  – rows = Z1/Z2 component, cols = true cluster 1/2/3
#   $hellinger : 3×3 matrix  – rows = p1/p2/p3 category, cols = true cluster 1/2/3
#
# The resulting table is a 9×5 matrix:
#   rows  : s1/s2/s3 (true cluster) × t300/t750/t2000 (time series length)
#   cols  : z1, z2, p1, p2, p3

num_indvs    <- 1000
scenario     <- "A"
num_replicas <- 100
est_choice   <- "multinomial"
suffix       <- "TRUE_neworder.RData"

data_folder <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/A_100_multinomial_hazel_table1"

t_lengths <- c(300, 750, 2000)

print(paste("Building MSE table for num_indvs =", num_indvs))
print(paste("Data folder:", data_folder))
print(paste("Time series lengths:", paste(t_lengths, collapse = ", ")))

# Load each file into its own environment to avoid name collisions.
load_est_values <- function(n, t) {
  fname <- paste("ClusterSim", n, t, scenario, num_replicas, est_choice, suffix, sep = "_")
  print(paste("Loading:", fname))
  e <- new.env()
  load(file.path(data_folder, fname), envir = e)
  stopifnot(!is.null(e$est_values$mse), !is.null(e$est_values$hellinger))
  e$est_values
}

results <- lapply(t_lengths, function(t) load_est_values(num_indvs, t))
names(results) <- paste0("t", t_lengths)

mse_row <- function(r, cluster_col) {
  c(r$mse[1, cluster_col],        # z1 MSE for that cluster
    r$mse[2, cluster_col],        # z2 MSE for that cluster
    r$hellinger[1, cluster_col],  # p1 Hellinger for that cluster
    r$hellinger[2, cluster_col],  # p2 Hellinger for that cluster
    r$hellinger[3, cluster_col])  # p3 Hellinger for that cluster
}

mse_table <- do.call(rbind, lapply(1:3, function(cluster_col) {
  do.call(rbind, lapply(results, function(r) mse_row(r, cluster_col)))
}))

rownames(mse_table) <- as.vector(outer(
  paste0("s", 1:3, "n", num_indvs),
  paste0("t", t_lengths),
  paste0
))
colnames(mse_table) <- c("z1", "z2", "p1", "p2", "p3")

ta <- t(mse_table)
ta <- ta[c("p1", "p2", "p3", "z1", "z2"), ]
print(round(ta, 2))
