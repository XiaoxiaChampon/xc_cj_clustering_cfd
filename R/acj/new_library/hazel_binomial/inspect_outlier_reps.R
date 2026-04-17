e <- new.env()
load("D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/RAW_HAZEL_BINOMIAL/A_100_binomial_hazel_table1/batch_0/ClusterSim_1000_300_A_100_binomial_TRUE_neworder.RData", envir = e)
obj <- e[["est_values"]]
reps <- obj[["rmse_reps"]]

cat("Replica 42 - rmse_reps [z1/z2 x cluster]:\n")
print(reps[42, , ])
cat("\nReplica 70 - rmse_reps [z1/z2 x cluster]:\n")
print(reps[70, , ])
cat("\nHellinger for rep 42:\n")
print(obj[["hellinger_reps"]][42, , ])
cat("\nHellinger for rep 70:\n")
print(obj[["hellinger_reps"]][70, , ])
