# export_timing_json_binomial.R
# Reads time_elapsed RData files and LSF log files for binomial runs,
# exports timing_summary_binomial.json and timing_summary_binomial.csv
#
# Decisions encoded:
#  - Job 480569 (t=750, A, n=1000): all 20 batches recorded in timing files
#  - Job 481036 (A, t=2000, n=1000): superseded by 481053 (same output folder
#    overwritten); RData unavailable — log data only, RData fields = NULL
#  - Job 480573 (test, A, t=2000, n=1000, 5 reps): log only, RData in a
#    non-RAW_HAZEL_BINOMIAL folder → RData fields = NULL

library(jsonlite)

# ── Paths ──────────────────────────────────────────────────────────────────
raw_hazel  <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/RAW_HAZEL_BINOMIAL"
log_dir    <- file.path(raw_hazel, "logs")
folder_a100_t1 <- file.path(raw_hazel, "A_100_binomial_hazel_table1")  # t=300 (batch_0), t=750 (batch_1..20)
folder_a50_t1  <- file.path(raw_hazel, "A_50_binomial_hazel_table1")   # t=2000, A, n=1000 (batch_1,2)
folder_a100_t2 <- file.path(raw_hazel, "A_100_binomial_hazel_table2")  # t=2000, A, n=100/500 (batch_0)
folder_b100_t2 <- file.path(raw_hazel, "B_100_binomial_hazel_table2")  # t=2000, B, n=100/500/1000 (no batch subdir)
out_json   <- file.path(raw_hazel, "timing_summary_binomial.json")
out_csv    <- file.path(raw_hazel, "timing_summary_binomial.csv")

# ── Log-file parser ────────────────────────────────────────────────────────
parse_log <- function(path) {
  txt  <- readLines(path, warn = FALSE)
  full <- paste(txt, collapse = "\n")

  m_r   <- regmatches(full, regexpr("Time difference of ([\\d\\.]+) (\\w+)", full, perl = TRUE))
  m_lsf <- regmatches(full, regexpr("Run time\\s*:\\s*(\\d+) sec", full, perl = TRUE))
  m_ari <- regmatches(full, regexpr("Mean ARI\\s*:\\s*([\\d\\.]+)", full, perl = TRUE))
  m_tlen <- regmatches(full, regexpr("Timeseries Len:\\s*(\\d+)", full, perl = TRUE))
  m_nrep <- regmatches(full, regexpr("Num Replicas:\\s*(\\d+)", full, perl = TRUE))
  m_nind <- regmatches(full, regexpr("Num Indvs:\\s*(\\d+)", full, perl = TRUE))
  m_out  <- regmatches(full, regexpr("Output folder:\\s*(\\S+)", full, perl = TRUE))
  m_scen <- regmatches(full, regexpr("Scenario:\\s*(\\w+)", full, perl = TRUE))

  list(
    r_time_str       = if (length(m_r))    m_r    else NA_character_,
    lsf_walltime_sec = if (length(m_lsf))  as.integer(sub(".*Run time\\s*:\\s*(\\d+) sec.*",  "\\1", m_lsf,  perl=TRUE)) else NA_integer_,
    mean_ari         = if (length(m_ari))  as.numeric(sub(".*Mean ARI\\s*:\\s*([\\d\\.]+).*", "\\1", m_ari,  perl=TRUE)) else NA_real_,
    t_len            = if (length(m_tlen)) as.integer(sub(".*Timeseries Len:\\s*(\\d+).*",    "\\1", m_tlen, perl=TRUE)) else NA_integer_,
    n_reps_reported  = if (length(m_nrep)) as.integer(sub(".*Num Replicas:\\s*(\\d+).*",      "\\1", m_nrep, perl=TRUE)) else NA_integer_,
    n_indvs          = if (length(m_nind)) as.integer(sub(".*Num Indvs:\\s*(\\d+).*",         "\\1", m_nind, perl=TRUE)) else NA_integer_,
    output_folder    = if (length(m_out))  sub(".*Output folder:\\s*(\\S+).*",                "\\1", m_out,  perl=TRUE)  else NA_character_,
    scenario         = if (length(m_scen)) sub(".*Scenario:\\s*(\\w+).*",                     "\\1", m_scen, perl=TRUE)  else NA_character_
  )
}

# ── time_elapsed RData reader ──────────────────────────────────────────────
read_time_elapsed <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  te    <- env$time_elapsed
  regen <- if (!is.null(env$total_regens)) as.integer(env$total_regens[[1L]]) else NA_integer_
  out   <- list(total_regens = regen)
  for (method in names(te)) out[[method]] <- as.numeric(te[[method]])
  out
}

# ── RData path lookup ─────────────────────────────────────────────────────
get_rdata_path <- function(job_id, batch_num) {
  switch(job_id,
    "480442" = file.path(folder_a100_t2, "batch_0",
                         "time_elapsed_100_2000_A_100_binomial_FALSE_neworder.RData"),
    "480467" = file.path(folder_a100_t2, "batch_0",
                         "time_elapsed_500_2000_A_100_binomial_FALSE_neworder.RData"),
    "480544" = file.path(folder_b100_t2,
                         "time_elapsed_500_2000_B_100_binomial_FALSE_neworder.RData"),
    "480547" = file.path(folder_b100_t2,
                         "time_elapsed_100_2000_B_100_binomial_FALSE_neworder.RData"),
    "480566" = file.path(folder_a100_t1, "batch_0",
                         "time_elapsed_1000_300_A_100_binomial_TRUE_neworder.RData"),
    "480569" = file.path(folder_a100_t1, paste0("batch_", batch_num),
                         "time_elapsed_1000_750_A_100_binomial_TRUE_neworder.RData"),
    "480573" = NULL,   # test job — output not in RAW_HAZEL_BINOMIAL
    "481036" = NULL,   # superseded by 481053 (same folder overwritten)
    "481053" = file.path(folder_a50_t1, paste0("batch_", batch_num),
                         "time_elapsed_1000_2000_A_50_binomial_TRUE_neworder.RData"),
    "482185" = file.path(folder_b100_t2,
                         "time_elapsed_1000_2000_B_100_binomial_FALSE_neworder.RData"),
    NULL
  )
}

# ── Job metadata ──────────────────────────────────────────────────────────
job_meta <- list(
  "480442" = list(scenario="A", n=100L,  t=2000L, cores=64L, reps_per_batch=100L, n_batches=1L,
                  table="table2", job_type="single",
                  note="Table2 A n=100, t=2000"),
  "480467" = list(scenario="A", n=500L,  t=2000L, cores=64L, reps_per_batch=100L, n_batches=1L,
                  table="table2", job_type="single",
                  note="Table2 A n=500, t=2000"),
  "480544" = list(scenario="B", n=500L,  t=2000L, cores=64L, reps_per_batch=100L, n_batches=1L,
                  table="table2", job_type="single",
                  note="Table2 B n=500, t=2000"),
  "480547" = list(scenario="B", n=100L,  t=2000L, cores=64L, reps_per_batch=100L, n_batches=1L,
                  table="table2", job_type="single",
                  note="Table2 B n=100, t=2000"),
  "480566" = list(scenario="A", n=1000L, t=300L,  cores=64L, reps_per_batch=100L, n_batches=1L,
                  table="table1", job_type="single",
                  note="Table1 A n=1000, t=300"),
  "480569" = list(scenario="A", n=1000L, t=750L,  cores=64L, reps_per_batch=100L, n_batches=20L,
                  table="table1", job_type="array",
                  note="Table1 A n=1000, t=750 — 20 batches x 100 reps (only batch_1 used in paper)"),
  "480573" = list(scenario="A", n=1000L, t=2000L, cores=64L, reps_per_batch=5L,  n_batches=1L,
                  table="none",   job_type="test",
                  note="Test job (5 reps); output not in RAW_HAZEL_BINOMIAL — log only"),
  "481036" = list(scenario="A", n=1000L, t=2000L, cores=64L, reps_per_batch=50L, n_batches=2L,
                  table="table1", job_type="array",
                  note="Table1 A n=1000, t=2000 — superseded: 481053 overwrote same output folders; RData unavailable"),
  "481053" = list(scenario="A", n=1000L, t=2000L, cores=64L, reps_per_batch=50L, n_batches=2L,
                  table="table1", job_type="array",
                  note="Table1 A n=1000, t=2000 — final run (2 batches x 50 reps = 100 total)"),
  "482185" = list(scenario="B", n=1000L, t=2000L, cores=64L, reps_per_batch=100L, n_batches=1L,
                  table="table2", job_type="single",
                  note="Table2 B n=1000, t=2000")
)

# ── Collect log files grouped by job_id ───────────────────────────────────
log_files    <- list.files(log_dir, pattern = "^out\\.", full.names = TRUE)
logs_by_job  <- list()

for (lf in sort(log_files)) {
  parts  <- strsplit(basename(lf), "\\.")[[1]]
  job_id <- parts[2]
  batch_n <- if (length(parts) >= 3) as.integer(parts[3]) else 1L
  if (!job_id %in% names(job_meta)) next
  if (is.null(logs_by_job[[job_id]])) logs_by_job[[job_id]] <- list()
  logs_by_job[[job_id]] <- c(logs_by_job[[job_id]], list(list(file = lf, batch = batch_n)))
}

# ── Build records ─────────────────────────────────────────────────────────
records <- list()

for (job_id in names(logs_by_job)) {
  meta        <- job_meta[[job_id]]
  job_batches <- list()

  for (entry in logs_by_job[[job_id]]) {
    log_info  <- parse_log(entry$file)
    batch_num <- entry$batch
    rdata_path <- get_rdata_path(job_id, batch_num)
    te_data    <- read_time_elapsed(rdata_path)

    batch_status <- if (is.na(log_info$mean_ari)) "FAILED_OR_TIMEOUT" else "OK"
    if (meta$job_type == "test")       batch_status <- paste0("TEST_", batch_status)
    if (job_id == "481036")            batch_status <- paste0("SUPERSEDED_", batch_status)

    job_batches <- c(job_batches, list(list(
      batch              = batch_num,
      log_file           = basename(entry$file),
      batch_status       = batch_status,
      rdata_file         = if (!is.null(rdata_path)) basename(rdata_path) else NULL,
      log_parsed = list(
        r_total_time_str  = log_info$r_time_str,
        lsf_walltime_sec  = log_info$lsf_walltime_sec,
        mean_ari          = log_info$mean_ari,
        t_len_from_log    = log_info$t_len,
        n_reps_from_log   = log_info$n_reps_reported,
        n_indvs_from_log  = log_info$n_indvs,
        output_folder     = log_info$output_folder
      ),
      time_elapsed_per_replica = te_data
    )))
  }

  records <- c(records, list(list(
    job_id         = job_id,
    scenario       = meta$scenario,
    n              = meta$n,
    t              = meta$t,
    model          = "binomial",
    cores          = meta$cores,
    reps_per_batch = meta$reps_per_batch,
    n_batches      = meta$n_batches,
    table          = meta$table,
    job_type       = meta$job_type,
    note           = meta$note,
    batches        = job_batches
  )))
}

# ── Write JSON ─────────────────────────────────────────────────────────────
output <- list(
  generated    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  description  = paste(
    "Per-replica timing data (Xiaoxia / univfpca / dbscan methods) and LSF job metrics",
    "for catFDA binomial clustering simulations (Tables 1 & 2, Scenarios A & B).",
    "time_elapsed values are in seconds. lsf_walltime_sec is the LSF-reported wall clock time.",
    "Job 480569 (t=750): all 20 batches recorded; only batch_1 used in the paper.",
    "Job 481036: superseded by 481053 (output folders overwritten); RData unavailable."
  ),
  source_folders = list(
    a100_table1 = folder_a100_t1,
    a50_table1  = folder_a50_t1,
    a100_table2 = folder_a100_t2,
    b100_table2 = folder_b100_t2,
    log_folder  = log_dir
  ),
  jobs = records
)

write(toJSON(output, pretty = TRUE, auto_unbox = TRUE, na = "null"), out_json)
cat("JSON written to:", out_json, "\n")

# ── Write CSV (one row per replica) ───────────────────────────────────────
rows <- list()

for (job in records) {
  for (b in job$batches) {
    te  <- b$time_elapsed_per_replica
    lp  <- b$log_parsed
    n_rep <- if (!is.null(te) && !is.null(te$Xiaoxia)) length(te$Xiaoxia) else 0L

    base_row <- list(
      job_id            = job$job_id,
      scenario          = job$scenario,
      n                 = job$n,
      t                 = job$t,
      model             = job$model,
      cores             = job$cores,
      reps_per_batch    = job$reps_per_batch,
      table             = job$table,
      job_type          = job$job_type,
      batch             = b$batch,
      log_file          = b$log_file,
      batch_status      = b$batch_status,
      lsf_walltime_sec  = lp$lsf_walltime_sec,
      mean_ari_batch    = lp$mean_ari,
      t_len_from_log    = lp$t_len_from_log,
      n_reps_from_log   = lp$n_reps_from_log,
      n_indvs_from_log  = lp$n_indvs_from_log
    )

    if (n_rep == 0L) {
      rows <- c(rows, list(as.data.frame(c(base_row, list(
        total_regens = NA_integer_,
        replica      = NA_integer_,
        Xiaoxia_sec  = NA_real_,
        univfpca_sec = NA_real_,
        dbscan_sec   = NA_real_
      ), stringsAsFactors = FALSE)))
      )
    } else {
      for (i in seq_len(n_rep)) {
        rows <- c(rows, list(as.data.frame(c(base_row, list(
          total_regens = if (!is.null(te$total_regens) && length(te$total_regens) > 0L) as.integer(te$total_regens[[1L]]) else NA_integer_,
          replica      = i,
          Xiaoxia_sec  = te$Xiaoxia[i],
          univfpca_sec = te$univfpca[i],
          dbscan_sec   = te$dbscan[i]
        ), stringsAsFactors = FALSE)))
        )
      }
    }
  }
}

csv_df <- do.call(rbind, rows)
write.csv(csv_df, out_csv, row.names = FALSE)
cat("CSV written to:", out_csv, "\n")
cat("Rows:", nrow(csv_df), "| Columns:", ncol(csv_df), "\n")
cat("Jobs:", length(records), "\n")
