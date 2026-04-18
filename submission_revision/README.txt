============================================================
  REPRODUCIBILITY MATERIALS — JASA Revision
  Categorical Functional Data Clustering with Twitter Application
============================================================

DATA AVAILABILITY
-----------------
Large pre-computed data files (*.RData) are NOT included in the code repository
due to file size. Download them from Zenodo before running any scripts:

  Zenodo DOI: https://doi.org/10.5281/zenodo.XXXXXXX   [to be updated on upload]

Place all downloaded .RData files in:
  submission_revision/data/processed/

The CSV files in data/raw/ and data/processed/ are included in the repository.


WORKING DIRECTORY
-----------------
Scripts use paths relative to their own location; no global working directory
needs to be set. Each script's header documents its expected working directory.


FOLDER STRUCTURE
----------------
submission_revision/
  README.txt                    <- this file
  data/
    data_dictionary.csv         <- variable descriptions for all data files
    raw/
      reference_tweet_data.txt  <- tweet IDs and company labels (9,802 tweets)
      users_profile_data.csv    <- anonymized profile metadata (4,776 users)
    processed/
      mentions_dataset.csv           <- binary activity matrix: @mention tweets
      no_mentions_dataset.csv        <- binary activity matrix: non-mention tweets
      no_tweet_dataset.csv           <- binary activity matrix: no tweet activity
      twitter_3875_clusterlabel.csv  <- cluster assignments for 3,875 users
      Twiiter_figure_logit_mul_final_Dec.RData  <- pre-computed estimates for figures
      W_matrix_final.RData           <- raw binary state matrix (optional; see Step 3)
  code/
    catfda_package/             <- R package source for catFDA method
    simulation/                 <- scripts to reproduce Tables 1 and 2
    twitter_application/        <- script for the Twitter real-data application
  outputs/
    figures/                    <- generated PNG figures (produced by Step 3)


============================================================
  STEP 0 — Install the catfda package
============================================================

The catFDA method is implemented as a local R package. Install it with:

  install.packages("devtools")
  devtools::install_local("code/catfda_package")

Required CRAN packages (install once):
  install.packages(c("mgcv", "pracma", "refund", "FADPclust",
                     "dbscan", "elbow", "fossil", "clickstream",
                     "ClickClust"))


============================================================
  STEP 1 — Reproduce Table 1 (Simulation MSE, Scenario A)
============================================================

Table 1 requires running 100 simulation replicas across three time series
lengths (t = 300, 750, 2000) with n = 1000 subjects. These were originally
run on an HPC cluster. To reproduce locally or on a cluster:

  1a. Submit the LSF/SLURM job scripts (or run the R scripts directly):
        code/simulation/1_run_jobs/job_A_n1000t300.sh    -> runs run_hazel_A_n1000t300.R
        code/simulation/1_run_jobs/job_A_n1000t750.sh    -> runs run_hazel_A_n1000t750.R
        code/simulation/1_run_jobs/job_A_n1000t2000.sh   -> runs run_hazel_A_n1000t2000.R

      Each .sh file requests 64 cores for 2 hours (LSF/bsub syntax). Adjust
      #BSUB headers and the conda activate path to match your HPC environment.
      Output goes to: outputs/clustersims/  (created automatically)

  1b. Combine batch outputs:
        source("code/simulation/2_combine/combine_table1_A.R")

  1c. Build the final Table 1:
        source("code/simulation/3_tables/table1_binomial.R")


============================================================
  STEP 2 — Reproduce Table 2 (Clustering ARI/RI, Scenarios A & B)
============================================================

Table 2 requires 100 replicas for each of n = 100, 500, 1000 under two
scenarios (A and B) at t = 2000.

  2a. Submit the job scripts (or run R scripts directly):
        code/simulation/1_run_jobs/job_A_n100t2000.sh
        code/simulation/1_run_jobs/job_A_n500t2000.sh
        code/simulation/1_run_jobs/job_A_n1000t2000.sh   (shared with Table 1)
        code/simulation/1_run_jobs/job_B_n100t2000.sh
        code/simulation/1_run_jobs/job_B_n500t2000.sh
        code/simulation/1_run_jobs/job_B_n1000t2000.sh

  2b. Combine batch outputs:
        source("code/simulation/2_combine/combine_table2_A.R")

  2c. Build the final Table 2:
        source("code/simulation/3_tables/table2_binomial.R")


============================================================
  STEP 3 — Reproduce Twitter Real-Data Application (Paper Figures 2 & 3)
============================================================

Run the following script from its own directory:

  cd code/twitter_application/
  Rscript Twitter_DataAnalysis_Final.R

Or from R:

  setwd("code/twitter_application")
  source("Twitter_DataAnalysis_Final.R")

Required data (already in data/processed/):
  - Twiiter_figure_logit_mul_final_Dec.RData   <- pre-computed estimates

Outputs written to outputs/figures/:
  - twitter_cluster_scatter.png   <- Paper Figure 2 (right panel): FPC score scatter
  - twitter_cluster1_curves.png   <- Paper Figure 3 (top row):    Cluster 1 curves
  - twitter_cluster2_curves.png   <- Paper Figure 3 (middle row): Cluster 2 curves
  - twitter_cluster0_curves.png   <- Paper Figure 3 (bottom row): Cluster 0 curves

Optional — individual trajectory plots:
  The script also contains a commented-out block that loads W_matrix_final.RData
  and plots example individual-level state trajectories using the cfda package.
  Uncomment that block and install cfda to use it.


============================================================
  DATA NOTES
============================================================

- All Twitter usernames have been replaced with anonymous pseudonymous
  identifiers (U_1 through U_8619). The U_xxx codes are consistent
  across reference_tweet_data.txt and users_profile_data.csv.

- Full tweet text is not included per X/Twitter Terms of Service.
  Tweet IDs in reference_tweet_data.txt can be used to retrieve tweet
  content via the Twitter/X API ("rehydration").

- Raw user timeline data (tweets_of_4776users) is not included:
  tweet IDs were not collected during data acquisition, and the data
  contains PII and full tweet text prohibited from redistribution.
  The processed binary matrices (mentions_dataset.csv, no_mentions_dataset.csv,
  no_tweet_dataset.csv) are the direct inputs to all statistical analyses
  and are provided in full.

- See data/data_dictionary.csv for a complete description of all variables
  in each data file.


============================================================
  SOFTWARE
============================================================

  R version:    4.x (developed and tested on R 4.3+)
  Key packages: catfda (local), mgcv, pracma, refund, FADPclust,
                dbscan, elbow, fossil, ggplot2



