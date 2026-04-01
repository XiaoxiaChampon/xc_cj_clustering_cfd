# export_timing_json.R
# Reads time_elapsed RData files from batch folders and LSF log files,
# then exports a comprehensive JSON with all timing/parameter details.

library(jsonlite)

# ── Paths ──────────────────────────────────────────────────────────────────
raw_hazel  <- "D:/PROJECTS/PAPERS/jasa_paper/xc_cj_clustering_cfd/outputs/clustersims/RAW_HAZEL"
log_dir    <- file.path(raw_hazel, "logs")
folder_5   <- file.path(raw_hazel, "A_5_multinomial_hazel_table1")   # t=750, t=2000 (5 reps/batch, 20 batches)
folder_10  <- file.path(raw_hazel, "A_10_multinomial_hazel_table1")  # t=300 (10 reps/batch, 10 batches)
out_json   <- file.path(raw_hazel, "timing_summary.json")

# ── Log-file parser ────────────────────────────────────────────────────────
parse_log <- function(path) {
  txt <- readLines(path, warn = FALSE)
  full <- paste(txt, collapse = "\n")

  # R total time ("Time difference of X.XX mins/secs")
  m_r <- regmatches(full, regexpr("Time difference of ([\\d\\.]+) (\\w+)", full, perl = TRUE))
  r_time_str <- if (length(m_r) == 1L) m_r else NA_character_

  # LSF wall time in seconds
  m_lsf <- regmatches(full, regexpr("Run time\\s*:\\s*(\\d+) sec", full, perl = TRUE))
  lsf_sec <- if (length(m_lsf) == 1L) {
    as.integer(sub(".*Run time\\s*:\\s*(\\d+) sec.*", "\\1", m_lsf, perl = TRUE))
  } else NA_integer_

  # Mean ARI
  m_ari <- regmatches(full, regexpr("Mean ARI\\s*:\\s*([\\d\\.]+)", full, perl = TRUE))
  mean_ari <- if (length(m_ari) == 1L) {
    as.numeric(sub(".*Mean ARI\\s*:\\s*([\\d\\.]+).*", "\\1", m_ari, perl = TRUE))
  } else NA_real_

  # T length from "Timeseries Len:"
  m_tlen <- regmatches(full, regexpr("Timeseries Len:\\s*(\\d+)", full, perl = TRUE))
  t_len <- if (length(m_tlen) == 1L) {
    as.integer(sub(".*Timeseries Len:\\s*(\\d+).*", "\\1", m_tlen, perl = TRUE))
  } else NA_integer_

  # Num replicas
  m_nrep <- regmatches(full, regexpr("Num Replicas:\\s*(\\d+)", full, perl = TRUE))
  n_reps <- if (length(m_nrep) == 1L) {
    as.integer(sub(".*Num Replicas:\\s*(\\d+).*", "\\1", m_nrep, perl = TRUE))
  } else NA_integer_

  list(
    r_time_str = r_time_str,
    lsf_walltime_sec = lsf_sec,
    mean_ari = mean_ari,
    t_len = t_len,
    n_reps_reported = n_reps
  )
}

# ── time_elapsed RData reader ──────────────────────────────────────────────
read_time_elapsed <- function(path) {
  if (!file.exists(path)) return(NULL)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  te  <- env$time_elapsed          # list: Xiaoxia / univfpca / dbscan
  regen <- env$total_regens        # integer
  out <- list(total_regens = regen)
  for (method in names(te)) {
    out[[method]] <- as.numeric(te[[method]])  # difftime → numeric seconds
  }
  out
}

# ── Job metadata ──────────────────────────────────────────────────────────
# Maps job_id → known parameters
job_meta <- list(
  "469942" = list(t = 2000L, cores = 64L, reps_per_batch = 5L,  n_batches = 20L, status = "array",  note = "original t=2000 job"),
  "469943" = list(t = 300L,  cores = 32L, reps_per_batch = 10L, n_batches = 10L, status = "array",  note = "original t=300 job"),
  "469946" = list(t = 750L,  cores = 32L, reps_per_batch = 5L,  n_batches = 20L, status = "array",  note = "original t=750 job"),
  "476240" = list(t = 300L,  cores = 32L, reps_per_batch = 10L, n_batches = 1L,  status = "rerun",  note = "rerun t=300 batch_6 (NA/NaN crash in original)"),
  "476621" = list(t = 2000L, cores = 64L, reps_per_batch = 5L,  n_batches = 1L,  status = "rerun",  note = "rerun t=2000 batch_1 (wall-time exceeded in original)")
)

# Map job_id → which batch folder directory + t-value file suffix
# A_5 holds t=750 and t=2000 (5 reps/batch); A_10 holds t=300 (10 reps/batch)
get_rdata_path <- function(t_val, batch_num) {
  if (t_val == 300L) {
    folder <- folder_10
    suffix <- sprintf("time_elapsed_1000_300_A_10_multinomial_TRUE_neworder.RData")
  } else if (t_val == 750L) {
    folder <- folder_5
    suffix <- sprintf("time_elapsed_1000_750_A_5_multinomial_TRUE_neworder.RData")
  } else if (t_val == 2000L) {
    folder <- folder_5
    suffix <- sprintf("time_elapsed_1000_2000_A_5_multinomial_TRUE_neworder.RData")
  } else {
    return(NULL)
  }
  file.path(folder, paste0("batch_", batch_num), suffix)
}

# ── Build per-job, per-batch records ──────────────────────────────────────
log_files <- list.files(log_dir, pattern = "^out\\.", full.names = TRUE)

# Group log files by job_id
logs_by_job <- list()
for (lf in sort(log_files)) {
  parts   <- strsplit(basename(lf), "\\.")[[1]]
  job_id  <- parts[2]
  batch_n <- if (length(parts) >= 3) as.integer(parts[3]) else 1L
  if (!job_id %in% names(job_meta)) next
  if (is.null(logs_by_job[[job_id]])) logs_by_job[[job_id]] <- list()
  logs_by_job[[job_id]] <- c(logs_by_job[[job_id]], list(list(file = lf, batch = batch_n)))
}

# Build full records
records <- list()

for (job_id in names(logs_by_job)) {
  meta <- job_meta[[job_id]]
  job_batches <- list()

  for (entry in logs_by_job[[job_id]]) {
    log_info  <- parse_log(entry$file)
    batch_num <- entry$batch

    # For reruns (single-file jobs), the batch number stored in the note
    # 476240 → batch_6 of t=300; 476621 → batch_1 of t=2000
    rdata_batch <- if (job_id == "476240") 6L else if (job_id == "476621") 1L else batch_num

    rdata_path <- get_rdata_path(meta$t, rdata_batch)
    te_data    <- read_time_elapsed(rdata_path)

    # Determine status
    batch_status <- if (is.na(log_info$mean_ari)) "FAILED_OR_TIMEOUT" else "OK"
    if (meta$status == "rerun") batch_status <- paste0("RERUN_", batch_status)

    job_batches <- c(job_batches, list(list(
      batch              = batch_num,
      log_file           = basename(entry$file),
      batch_status       = batch_status,
      rdata_batch_folder = if (!is.null(rdata_path)) paste0("batch_", rdata_batch) else NULL,
      rdata_file         = if (!is.null(rdata_path)) basename(rdata_path) else NULL,
      log_parsed = list(
        r_total_time_str  = log_info$r_time_str,
        lsf_walltime_sec  = log_info$lsf_walltime_sec,
        mean_ari          = log_info$mean_ari,
        t_len_from_log    = log_info$t_len,
        n_reps_from_log   = log_info$n_reps_reported
      ),
      time_elapsed_per_replica = te_data
    )))
  }

  records <- c(records, list(list(
    job_id         = job_id,
    scenario       = "A",
    n              = 1000L,
    t              = meta$t,
    model          = "multinomial",
    cores          = meta$cores,
    reps_per_batch = meta$reps_per_batch,
    n_batches      = meta$n_batches,
    job_status     = meta$status,
    note           = meta$note,
    batches        = job_batches
  )))
}

# ── Write JSON ─────────────────────────────────────────────────────────────
output <- list(
  generated    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  description  = paste(
    "Per-replica timing data (Xiaoxia / univfpca / dbscan methods) and LSF job metrics",
    "for catFDA multinomial clustering simulations (Table 1, Scenario A, n=1000).",
    "time_elapsed values are in seconds. lsf_walltime_sec is the LSF-reported wall clock time."
  ),
  source_folders = list(
    a5_batch_folder  = folder_5,
    a10_batch_folder = folder_10,
    log_folder       = log_dir
  ),
  jobs = records
)

write(toJSON(output, pretty = TRUE, auto_unbox = TRUE, na = "null"), out_json)
cat("Written to:", out_json, "\n")
cat("Jobs:", length(records), "\n")

# ── Write CSV (one row per replica) ───────────────────────────────────────
out_csv <- file.path(raw_hazel, "timing_summary.csv")

rows <- list()
for (job in records) {
  for (b in job$batches) {
    te  <- b$time_elapsed_per_replica
    lp  <- b$log_parsed
    # number of replicas = length of Xiaoxia vector (or 0 if te is NULL)
    n_rep <- if (!is.null(te) && !is.null(te$Xiaoxia)) length(te$Xiaoxia) else 0L
    if (n_rep == 0L) {
      # Batch failed / no RData: one summary row with NA timings
      rows <- c(rows, list(data.frame(
        job_id            = job$job_id,
        scenario          = job$scenario,
        n                 = job$n,
        t                 = job$t,
        model             = job$model,
        cores             = job$cores,
        reps_per_batch    = job$reps_per_batch,
        job_status        = job$job_status,
        batch             = b$batch,
        log_file          = b$log_file,
        batch_status      = b$batch_status,
        lsf_walltime_sec  = lp$lsf_walltime_sec,
        mean_ari_batch    = lp$mean_ari,
        t_len_from_log    = lp$t_len_from_log,
        n_reps_from_log   = lp$n_reps_from_log,
        total_regens      = NA_integer_,
        replica           = NA_integer_,
        Xiaoxia_sec       = NA_real_,
        univfpca_sec      = NA_real_,
        dbscan_sec        = NA_real_,
        stringsAsFactors  = FALSE
      )))
    } else {
      for (i in seq_len(n_rep)) {
        rows <- c(rows, list(data.frame(
          job_id            = job$job_id,
          scenario          = job$scenario,
          n                 = job$n,
          t                 = job$t,
          model             = job$model,
          cores             = job$cores,
          reps_per_batch    = job$reps_per_batch,
          job_status        = job$job_status,
          batch             = b$batch,
          log_file          = b$log_file,
          batch_status      = b$batch_status,
          lsf_walltime_sec  = lp$lsf_walltime_sec,
          mean_ari_batch    = lp$mean_ari,
          t_len_from_log    = lp$t_len_from_log,
          n_reps_from_log   = lp$n_reps_from_log,
          total_regens      = if (!is.null(te$total_regens) && length(te$total_regens) > 0L) as.integer(te$total_regens[[1L]]) else NA_integer_,
          replica           = i,
          Xiaoxia_sec       = te$Xiaoxia[i],
          univfpca_sec      = te$univfpca[i],
          dbscan_sec        = te$dbscan[i],
          stringsAsFactors  = FALSE
        )))
      }
    }
  }
}

csv_df <- do.call(rbind, rows)
write.csv(csv_df, out_csv, row.names = FALSE)
cat("CSV written to:", out_csv, "\n")
cat("Rows:", nrow(csv_df), "| Columns:", ncol(csv_df), "\n")
