# ============================================================
# Scenario A — Clustering Simulation
# Cluster fractions:  75% / 22% / 3%
# mu_2 (setting 1):   -2.5 + exp(2t)
# Score variances:    1, 1/2, 1/4
# Users (n):          100, 500, 1000
# Time series (t):    300, 750, 2000
# Focus:              ARI and RI scores
# ============================================================
setwd("D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd")

.just_load_functions <- TRUE
source("R/acj/new_library/acj_clustering_newlib_testvalues.R")
.just_load_functions <- FALSE
source("R/acj/new_library/run_scenario_helper.R")

set.seed(123)

scenario     <- "A"
num_replicas <- 2
est_choice   <- "multinomial"
n_values     <- c(100, 500, 1000)
t_values     <- c(300, 750, 2000)

temp_folder <- file.path("outputs", "clustersims",
                         paste0("A_", num_replicas, "_", est_choice, "_neworder"))
if (dir.exists(temp_folder)) unlink(temp_folder, recursive = TRUE)
dir.create(temp_folder, recursive = TRUE)
cat("Output folder:", temp_folder, "\n")

results_A <- RunScenarioSims(scenario, n_values, t_values,
                              num_replicas, est_choice, temp_folder)

if (run_parallel) {
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
