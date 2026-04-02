# ============================================================
# build_table2.R
# Assemble Table 2: ARI and RI for Scenarios A and B | t=2000, n=100/500/1000
#
# Data sources (all are 100-replica summaries in est_values$cluster_table_est):
#
#   Scenario A:
#     n=100  : outputs/clustersims/A_100_binomial_hazel_table2/
#                ClusterSim_100_2000_A_100_binomial_FALSE_neworder.RData
#              (produced by combine_table2_A.R from Table2Files)
#     n=500  : same folder, n=500 variant
#     n=1000 : outputs/clustersims/A_100_binomial_hazel_table1/
#                ClusterSim_1000_2000_A_100_binomial_TRUE_neworder.RData
#              (produced by combine_table1_A.R after array batch jobs)
#
#   Scenario B:
#     n=100  : outputs/PaperTables/Table2Files/B_100_binomial_hazel_table2/
#                ClusterSim_100_2000_B_100_binomial_FALSE_neworder.RData
#     n=500  : same folder, n=500 variant
#     n=1000 : outputs/PaperTables/B_100_binomial_hazel/
#                ClusterSim_1000_2000_B_100_binomial_FALSE_neworder.RData
# ============================================================

proj_root <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd"

# Helper: load est_values from a file; returns NULL (with warning) if missing
load_ev <- function(path) {
  if (!file.exists(path)) {
    warning("File not found (will appear as NA in table): ", path)
    return(NULL)
  }
  e <- new.env()
  load(path, envir = e)
  if (is.null(e$est_values))
    stop("est_values not found in: ", path)
  e$est_values
}

# Helper: extract metric from ev; returns NA list if ev is NULL
extract_safe <- function(ev, metric) {
  if (is.null(ev)) return(list(mean = NA_real_, se = NA_real_))
  list(
    mean = unname(ev$cluster_table_est[metric]),
    se   = unname(ev$cluster_table_est_se[metric])
  )
}

# ---------- File paths ----------
file_A100 <- file.path(proj_root, "outputs", "clustersims",
  "A_100_binomial_hazel_table2",
  "ClusterSim_100_2000_A_100_binomial_FALSE_neworder.RData")

file_A500 <- file.path(proj_root, "outputs", "clustersims",
  "A_100_binomial_hazel_table2",
  "ClusterSim_500_2000_A_100_binomial_FALSE_neworder.RData")

file_A1000 <- file.path(proj_root, "outputs", "clustersims",
  "A_100_binomial_hazel_table1",
  "ClusterSim_1000_2000_A_100_binomial_TRUE_neworder.RData")

file_B100 <- file.path(proj_root, "outputs", "clustersims",
  "RAW_HAZEL_BINOMIAL", "B_100_binomial_hazel_table2",
  "ClusterSim_100_2000_B_100_binomial_FALSE_neworder.RData")

file_B500 <- file.path(proj_root, "outputs", "clustersims",
  "RAW_HAZEL_BINOMIAL", "B_100_binomial_hazel_table2",
  "ClusterSim_500_2000_B_100_binomial_FALSE_neworder.RData")

# n=1000 Scenario B — populated after job_B_n1000t2000 runs on Hazel
file_B1000 <- file.path(proj_root, "outputs", "clustersims",
  "RAW_HAZEL_BINOMIAL", "B_100_binomial_hazel_table2",
  "ClusterSim_1000_2000_B_100_binomial_FALSE_neworder.RData")

# ---------- Load ----------
cat("Loading results...\n")
ev_A100  <- load_ev(file_A100)
ev_A500  <- load_ev(file_A500)
ev_A1000 <- load_ev(file_A1000)
ev_B100  <- load_ev(file_B100)
ev_B500  <- load_ev(file_B500)
ev_B1000 <- load_ev(file_B1000)  # NULL if not yet available

# ---------- Extract ARI and RI (mean and SE) ----------
cells <- list(
  A_n100  = list(ARI = extract_safe(ev_A100,  "dbscan ARI"), RI = extract_safe(ev_A100,  "dbscan RI")),
  A_n500  = list(ARI = extract_safe(ev_A500,  "dbscan ARI"), RI = extract_safe(ev_A500,  "dbscan RI")),
  A_n1000 = list(ARI = extract_safe(ev_A1000, "dbscan ARI"), RI = extract_safe(ev_A1000, "dbscan RI")),
  B_n100  = list(ARI = extract_safe(ev_B100,  "dbscan ARI"), RI = extract_safe(ev_B100,  "dbscan RI")),
  B_n500  = list(ARI = extract_safe(ev_B500,  "dbscan ARI"), RI = extract_safe(ev_B500,  "dbscan RI")),
  B_n1000 = list(ARI = extract_safe(ev_B1000, "dbscan ARI"), RI = extract_safe(ev_B1000, "dbscan RI"))  # NA if pending
)

# ---------- Build display table ----------
# Rows: n=100, n=500, n=1000
# Cols: Scenario A (ARI, SE, RI, SE) | Scenario B (ARI, SE, RI, SE)
fmt <- function(m, s) if (is.na(m)) "pending" else sprintf("%.4f (%.4f)", m, s)

table2 <- data.frame(
  n = c(100, 500, 1000),
  A_ARI = c(fmt(cells$A_n100$ARI$mean,  cells$A_n100$ARI$se),
            fmt(cells$A_n500$ARI$mean,  cells$A_n500$ARI$se),
            fmt(cells$A_n1000$ARI$mean, cells$A_n1000$ARI$se)),
  A_RI  = c(fmt(cells$A_n100$RI$mean,   cells$A_n100$RI$se),
            fmt(cells$A_n500$RI$mean,   cells$A_n500$RI$se),
            fmt(cells$A_n1000$RI$mean,  cells$A_n1000$RI$se)),
  B_ARI = c(fmt(cells$B_n100$ARI$mean,  cells$B_n100$ARI$se),
            fmt(cells$B_n500$ARI$mean,  cells$B_n500$ARI$se),
            fmt(cells$B_n1000$ARI$mean, cells$B_n1000$ARI$se)),
  B_RI  = c(fmt(cells$B_n100$RI$mean,   cells$B_n100$RI$se),
            fmt(cells$B_n500$RI$mean,   cells$B_n500$RI$se),
            fmt(cells$B_n1000$RI$mean,  cells$B_n1000$RI$se)),
  stringsAsFactors = FALSE
)

cat("\n========== Table 2: catFDA-dbscan Clustering (t = 2000) ==========\n")
cat("Format: mean (SE) over 100 replicas\n\n")
print(table2, row.names = FALSE)

# Also save the raw numeric values for downstream use
table2_raw <- data.frame(
  n         = c(100, 500, 1000),
  A_ARI     = c(cells$A_n100$ARI$mean,  cells$A_n500$ARI$mean,  cells$A_n1000$ARI$mean),
  A_ARI_se  = c(cells$A_n100$ARI$se,    cells$A_n500$ARI$se,    cells$A_n1000$ARI$se),
  A_RI      = c(cells$A_n100$RI$mean,   cells$A_n500$RI$mean,   cells$A_n1000$RI$mean),
  A_RI_se   = c(cells$A_n100$RI$se,     cells$A_n500$RI$se,     cells$A_n1000$RI$se),
  B_ARI     = c(cells$B_n100$ARI$mean,  cells$B_n500$ARI$mean,  cells$B_n1000$ARI$mean),
  B_ARI_se  = c(cells$B_n100$ARI$se,    cells$B_n500$ARI$se,    cells$B_n1000$ARI$se),
  B_RI      = c(cells$B_n100$RI$mean,   cells$B_n500$RI$mean,   cells$B_n1000$RI$mean),
  B_RI_se   = c(cells$B_n100$RI$se,     cells$B_n500$RI$se,     cells$B_n1000$RI$se)
)

out_csv <- file.path(proj_root, "outputs", "table2_binomial_raw.csv")
write.csv(table2_raw, out_csv, row.names = FALSE)
cat("\nRaw values saved to:", out_csv, "\n")
