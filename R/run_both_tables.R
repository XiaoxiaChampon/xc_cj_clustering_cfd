# setwd("D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd")

cat("\n===== OLD (inline code) =====\n")
identifier_str <- "to_compare_old_code"
suffix <- "TRUE.RData"
source("R/build_mse_tableC.R")

cat("\n===== NEW (catfda library) =====\n")
rm(list = c("identifier_str", "suffix", "sub_folder", "data_folder"))
identifier_str <- "to_compare_new_code"
suffix <- "TRUE_neworder.RData"
source("R/build_mse_tableC.R")
