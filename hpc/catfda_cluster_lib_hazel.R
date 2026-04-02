# ============================================================
# catfda_cluster_lib_hazel.R
# Library file: all function definitions for catFDA clustering.
# Source this file from a runner script; do NOT run directly.
# ============================================================

# For: profiling and visualization of profiling
#library(profvis)

# For: gam
library(mgcv)

# For: cubicspline
library(pracma)

# For: fpca.face
library(refund)

# For: FADPclust
library(FADPclust)

# For: kNNdist
library(dbscan)

# For: elbow
# devtools::install_github("ahasverus/elbow")
library(elbow)

# For: rand.index
library(fossil)

# For: cfda method
library(cfda)

# For: gather method
library(tidyverse)

# For: set_names (used in new_order.R functions)
library(rlang)

# Load the our clustering library
#devtools::load_all("D:/PROJECTS/PAPERS/jasa_paper/catfda")
library(catfda)

# Source required functions
source("../trapzfnum_function.R")
# Source the improved estimation functions
#source("R/XCATFDA/R/new_order.R")

# ---- For: parallelization ----
# For: foreach loop
library(foreach)

run_parallel <- TRUE
time_elapsed <- list()
if (run_parallel) {
  print("RUNNING PARALLEL")

  # For: makeCluster
  library(doParallel)

  # For: %dorng% or registerDoRNG for reproducable parallel random number generation
  library(doRNG)

  if (exists("initialized_parallel") && initialized_parallel == TRUE) {
    parallel::stopCluster(cl = my.cluster)
  }
  
  n.cores <- suppressWarnings(as.numeric(Sys.getenv("LSB_DJOB_NUMPROC")))
  if (is.na(n.cores)) {
    detected <- parallel::detectCores()
    n.cores <- if (!is.na(detected) && detected > 1L) detected - 1L else 1L
  }
  print(paste("Using", n.cores, "cores"))
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  cat("Parellel Registered: ", foreach::getDoParRegistered(), " [cores = ", n.cores, "]\n")
  initialized_parallel <- TRUE

  # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
}

#' Create directories
if (!dir.exists("outputs")) {
  dir.create("outputs")
}
if (!dir.exists("outputs/clustersims")) {
  dir.create("outputs/clustersims")
}

# profvis({


#**** evaluate cluster accuracy using RI and ARI
# INPUT
# VECTORS with true memberships (true_cluster)
#        with estimated memberships (coming from cluster_method$label)
# noise: scalar, 0 or 3, if dbscan clustering, use 0, else, use 3.
# OUTPUT: RI and ARI, CPN
#
evaluate_cluster <- function(true_cluster, new_cluster, noise) {
  ri <- rand.index(true_cluster, new_cluster)
  ari <- adj.rand.index(true_cluster, new_cluster)

  c3 <- which(true_cluster == noise)
  cc3 <- which(new_cluster == noise)
  cpn <- sum(cc3 %in% c3) / length(c3)
  return(list(ri = ri, ari = ari, cpn = cpn))
}

#**** clustering the scores using KMEANS
# INPUT
# n by K scores - or features extracted from the multivariate FD
# OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
kmeans_cluster <- function(data = scores_K) {
  out_kmeans <- NbClust::NbClust(
    data = data, diss = NULL,
    distance = "euclidean", min.nc = 2, max.nc = 5,
    method = "kmeans", index = "silhouette"
  )

  return(list(nclust = as.numeric(out_kmeans$Best.nc[1]), label = out_kmeans$Best.partition))
}


#**** clustering the scores using KMEANS
# INPUT
# n by K scores - or features extracted from the multivariate FD
# scale_eps, scalar, the scale of the epsilon
# OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
dbscan_cluster <- function(data = scores_K, scale_eps) {
  dimz <- dim(data)[2]
  if (dimz <= 2) {
    minPts <- 4
  }
  if (dimz > 2) {
    minPts <- 2 * dimz + 1
  }
  dist <- kNNdist(data, k = minPts - 1)
  ######### change to max increase
  distdataelbow <- data.frame(sort(dist))
  distdataelbow$index <- 1:(dim(data)[1])
  ipoint <- elbow(data = distdataelbow, plot = FALSE)
  epsoptimal <- (ipoint$sort.dist._selected) * scale_eps

  out_dbscan <- dbscan(data, eps = epsoptimal, minPts = minPts)
  return(list(nclust = dim(table(out_dbscan$cluster)), label = out_dbscan$cluster))
}


#*** clustering the curves directly using FADP
# INPUT
# Z1,Z2 - m by n latent processes for the categ FD with 3 categories
# tt - grid of points
# PVE - percentage of explained variance (default value 0.95)
# OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
fadp_cluster <- function(mZ1, mZ2, tt = tt, PVE = 0.95) {
  fdabasis <- fda::create.bspline.basis(c(0, 1), 25, 4)
  fdatime <- tt
  fdafd1 <- fda::smooth.basis(tt, mZ1, fdabasis)$fd
  fdafd2 <- fda::smooth.basis(tt, mZ2, fdabasis)$fd

  FADlist <- list(fdafd1, fdafd2)

  FADP <- FADPclust(
    fdata = FADlist, cluster = 2:5, method = "FADP1",
    proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
    pve = PVE, stats = "Avg.silhouette"
  )

  return(list(nclust = FADP$nclust, label = FADP$clust))
}



#*** extract scores using UNIVARIAFE fpca, FAST COV estimation.  SUPER FAST
# INPUT
# Z1,Z2 - m by n latent processes for the categ FD with 3 categories
# tt - grid of points
# PVE - percentage of explained variance (default value 0.95)
# OUTPUT
# list with estimated PHI (m by K) with function norm 1
#           projection of data onto PHI (called scores - nby K matrix)
#

extract_scores_UNIVFPCA <- function(mZ1, mZ2, tt, PVE = 0.95) {
  m <- nrow(mZ1)
  n <- ncol(mZ1)

  out1 <- fpca.face(Y = t(mZ1), argvals = tt, pve = 0.99)
  out2 <- fpca.face(Y = t(mZ2), argvals = tt, pve = 0.99)
  O_m1 <- matrix(0, nrow = m, ncol = out1$npc)
  O_m2 <- matrix(0, nrow = m, ncol = out2$npc)

  # construct PHI
  Phi_est0 <- rbind(
    cbind(out1$efunctions * sqrt(m), O_m2),
    cbind(O_m1, out2$efunctions * sqrt(m))
  )
  Scores0 <- cbind(out1$scores, out2$scores)
  ScoresCov_0 <- cov(Scores0)
  oute <- eigen(ScoresCov_0)

  K <- which(cumsum(oute$values) / sum(oute$values) >= PVE)[1]
  count_iter <- 0
  delta <- 0.01
  while (K < 2 && count_iter < 100) {
    count_iter <- count_iter + 1
    cat("count_iter: ", count_iter, "\n")
    K <- which(cumsum(oute$values) / sum(oute$values) >= (PVE + delta))[1]
    delta <- delta + 0.01
  }

  Phi_est <- Phi_est0 %*% oute$vectors[, 1:K] # correct eigenfns

  mZ <- rbind(mZ1, mZ2)
  Scores_est <- t(mZ) %*% Phi_est / sqrt(m) # they are not demeaned

  return(list(scores = Scores_est, Phi = Phi_est))
}

# Function to use trapzfnum function and find L2 distance for 2D array, n of them
#' @param  truecurve- 2D array (true curve)
#' @param   estcurve-2D array (estimated curve)
#' @return : scalar-L2 distance

mse_bw_matrix <- function(truecurve, estcurve, timestamps01) {
  n <- dim(truecurve)[2]
  # datapoints=dim(truecurve)[1]
  # mseall=c(0)
  ###### could probably use apply function here it's also subject level
  mseall <- foreach(i = 1:n, .combine = c, .packages = c("pracma"), .export=c("trapzfnum")) %dorng% {
    #source("../trapzfnum_function.R")
    return(rbind(trapzfnum(truecurve[, i], estcurve[, i], timestamps01)))
  }

  return(mseall)
}


# Function to use trapzfnump function and find Hellinger distance for 2D array, n of them
mse_bw_matrixp <- function(truecurve, estcurve, timestamps01) {
  n <- dim(truecurve)[2]
  # datapoints= dim(truecurve)[1]
  # mseall=c(0)
  # ######could probably use apply function here it's also subject level
  # for (i in 1:n)
  # {
  #   mseall[i]=trapzfnump(truecurve[,i],estcurve[,i])/sqrt(2)
  # }

  sqrt_2 <- sqrt(2)
  mseall <- foreach(i = 1:n, .combine = c, .packages = c("pracma"), .export=c("trapzfnump")) %dorng% {
    #source("../trapzfnum_function.R")
    return(rbind(trapzfnump(truecurve[, i], estcurve[, i], timestamps01) / sqrt_2))
  }

  return(mseall)
}


# cfda
####
# write a function to produce the scores
# input one  X_nt matrix n*t: 1 of  Nth simulation, n*t* N
# output scores matrix for that specific Nth simulation   n*M     M is the column number of scores
cfda_score_function <- function(cfda_data, nCores, timestamps01, basis_num) {
  cfda_data <- t(cfda_data)
  timeseries_length <- dim(cfda_data)[2]
  times <- seq(1, timeseries_length, by = 1)
  colnames(cfda_data) <- paste0("time_", times)
  time <- rep(timestamps01, times = dim(cfda_data)[1])
  subjects <- seq(1, dim(cfda_data)[1], by = 1)
  cfda_new <- data.frame(cbind(cfda_data, subjects))
  data_long <- gather(cfda_new, time_t, state, paste0("time_", times[1]):paste0("time_", times[timeseries_length]), factor_key = TRUE)
  cfda_d <- data_long %>%
    arrange(subjects) %>%
    dplyr::rename(id = subjects) %>%
    dplyr::select(id, state)
  cfda_final <- data.frame(cbind(cfda_d, time))
  Tmax <- max(cfda_final$time)
  cfda_cut <- cut_data(cfda_final, Tmax = Tmax)
  b <- create.bspline.basis(c(0, Tmax), nbasis = basis_num, norder = 4)
  fmca <- compute_optimal_encoding(cfda_cut, b, nCores = nCores, verbose = FALSE)
  delta <- 0.01
  pve <- 0.95
  nPc90 <- which(cumsum(prop.table(fmca$eigenvalues)) > pve)[1]
  while (nPc90 < 2 && pve < 1) {
    pve <- pve + delta
    nPc90 <- which(cumsum(prop.table(fmca$eigenvalues)) > pve)[1]
  }
  cfda_score <- fmca$pc[, 1:nPc90]
  return(cfda_score)
}


ClusterSimulation <- function(num_indvs, timeseries_length,
                              scenario, num_replicas, est_choice_list, run_hellinger, temp_folder,
                              run_univfpca = TRUE, run_kmeans = TRUE, run_fadp = TRUE,
                              run_dbscan = TRUE, run_cfda = TRUE,
                              save_curves = TRUE) {
  cat(
    "Cluster Simulation\nNum Indvs:\t", num_indvs,
    "\nTimeseries Len:\t", timeseries_length,
    "\nScenario:\t", scenario,
    "\nNum Replicas:\t", num_replicas
  )

  occur_fraction <- GetOccurrenceFractions(scenario)
  cat("\nOccurrence Fractions: ", occur_fraction)

  cluster_allocation <- occur_fraction * num_indvs
  cat("\nIndividuals in each cluster (cluster alloc.): ", cluster_allocation)

  true_cluster <- rep(c(1, 2, 3), cluster_allocation)
  cat("\nTrue Cluster:\n", true_cluster)

  true_cluster_db <- rep(c(1, 2, 0), cluster_allocation)
  cat("\nTrue Cluster DB:\n", true_cluster_db)

  rmse <- array(0, c(num_replicas, 2, 3)) # 2 components and 3 settings
  hellinger <- array(0, c(num_replicas, 3, 3)) # 3 events probab and 3 settings
  true_kmeans <- est_kmeans <- NULL # records clustering membership on the TRUE SCORES/ESTIMATED SCORES
  true_fadp <- est_fadp <- NULL # records clustering membership on the TRUE Z/ESTIMATED Z
  true_dbscan <- est_dbscan <- NULL
  true_dbscan_cfda <- est_dbscan_cfda <- NULL

  time_elapsed <<- list()
  # "Xiaoxia"=NULL, "univfpca"=NULL, "kmeans"=NULL, "fadp"=NULL, "dbscan"=NULL, "cfd"=NULL)
  last_time <- 0
  row_name <- NULL
  timeKeeperStart <- function(rn) {
    row_name <<- rn
    if (FALSE == row_name %in% names(time_elapsed)) {
      time_elapsed[[row_name]] <<- NULL
    }
    last_time <<- Sys.time()
  }
  timeKeeperNext <- function() {
    this_time <- Sys.time()
    this_section_time <- this_time - last_time
    cat(row_name, "calc time taken:", capture.output(this_section_time), "\n")
    time_elapsed[[row_name]] <<- append(time_elapsed[[row_name]], this_section_time)
    last_time <<- this_time
  }

  ######
  Z_true_curve <- Z_est_curve <- list()
  p_true_curve <- p_est_curve <- list()
  W_cfd <- list()
  total_regens <- 0
  for (replica_idx in 1:num_replicas)
  {
    cat(
      "\nNum Indvs:", num_indvs,
      "\tTimeseries Len:", timeseries_length,
      "\tScenario:", scenario,
      "\tNum Replicas:", num_replicas
    )
    cat("\nReplica: ", replica_idx, "\n")

    # generate clusters
    # set.seed(seed_cluster + 100 * replica_idx)
    # cat("\nCluster", replica_idx, " --> seed: ", seed_cluster + 100 * replica_idx, "\n")
    cat("Cluster", replica_idx, "\n")

    cluster_f1 <- GenerateClusterData(1, scenario, 3, cluster_allocation[1], timeseries_length)
    cluster_f2 <- GenerateClusterData(2, scenario, 3, cluster_allocation[2], timeseries_length)
    cluster_f3 <- GenerateClusterData(3, scenario, 3, cluster_allocation[3], timeseries_length)

    # Z1, Z2: the two latent Gaussian processes. Always exactly 2 because the
    # multinomial logit model with 3 categories requires K-1=2 latent processes
    # (category 3 is the reference with implicit log-odds = 0).
    Z1 <- cbind(cluster_f1$Z1, cluster_f2$Z1, cluster_f3$Z1)
    Z2 <- cbind(cluster_f1$Z2, cluster_f2$Z2, cluster_f3$Z2)

    # p1, p2, p3: true category probability curves. Always 3 to match the
    # 3-category simulation setup: p1=exp(Z1)/d, p2=exp(Z2)/d, p3=1/d
    # where d = 1+exp(Z1)+exp(Z2).
    p1 <- cbind(cluster_f1$p1, cluster_f2$p1, cluster_f3$p1)
    p2 <- cbind(cluster_f1$p2, cluster_f2$p2, cluster_f3$p2)
    p3 <- cbind(cluster_f1$p3, cluster_f2$p3, cluster_f3$p3)
    prob_curves <- list(p1 = p1, p2 = p2, p3 = p3)
    ######### 9/11/2023
    if (save_curves) {
      Z_true_curve[[replica_idx]] <- array(c(Z1, Z2), dim = c(timeseries_length, num_indvs, 2))
      p_true_curve[[replica_idx]] <- array(c(p1, p2, p3), dim = c(timeseries_length, num_indvs, 3))
    }
    ############


    # generate categFuncData
    # set.seed(seed_cfd + 100 * replica_idx)
    # cat("\nCategFD", replica_idx, " --> seed: ", seed_cfd + 100 * replica_idx, "\n")
    cat("CategFD", replica_idx, "\n")

    categ_func_data_list <- generate_categ_func_data(prob_curves)

    # category_labels: the sorted set of distinct category values present in the
    # data (e.g. {1, 2, 3}). Sorted so one-hot encoding indices are consistent.
    category_labels <- unique(c(categ_func_data_list$w_mat))
    if (is.numeric(category_labels)) {
      category_labels <- sort(category_labels)
    }

    # Data quality check: for each individual, verify their categorical time
    # series has all categories present and each appears often enough for reliable
    # estimation. If not, regenerate just that individual (within their cluster).
    for (indv in 1:num_indvs)
    {
      if (indv %in% 1:cluster_allocation[1]) {
        setting_choice <- 1
      }
      if (indv %in% (cluster_allocation[1] + 1):(cluster_allocation[1] + cluster_allocation[2])) {
        setting_choice <- 2
      }
      if (indv %in% (cluster_allocation[1] + cluster_allocation[2] + 1):num_indvs) {
        setting_choice <- 3
      }

      # cat_counts: frequency table of each category for this individual.
      # n_observed_cats: number of distinct categories actually seen.
      # Regeneration thresholds ensure every category is observed enough times
      # for robust estimation (roughly >=1.3% of time points for each category):
      #   t=300  -> min 4 occurrences per category
      #   t=750  -> min 10 occurrences per category
      cat_counts <- table(categ_func_data_list$w_mat[, indv])
      n_observed_cats <- length(cat_counts)
      count_iter <- 0
      while (count_iter < 100 &&
        ((n_observed_cats < length(category_labels)) ||
          (timeseries_length == 300 && min(as.numeric(cat_counts)) < 4) ||
          (timeseries_length == 750 && min(as.numeric(cat_counts)) < 10)
        )
      ) {
        count_iter <- count_iter + 1

        new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)

        new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
        new_categ_func_data_list <- generate_categ_func_data(new_prob_curves)

        # GenerateClusterData produces 5 candidate individuals; take the 3rd
        # (middle candidate) to avoid boundary effects from the first/last draws.
        categ_func_data_list$w_mat[, indv] <- new_categ_func_data_list$w_mat[, 3]
        Z1[, indv] <- new_cluster_data$Z1[, 3]
        categ_func_data_list$x_array[indv, , ] <- 0
        Z2[, indv] <- new_cluster_data$Z2[, 3]

        for (this_time in 1:timeseries_length)
        {
          categ_func_data_list$x_array[indv, this_time, which(category_labels == categ_func_data_list$w_mat[, indv][this_time])] <- 1
        }

        cat_counts <- table(categ_func_data_list$w_mat[, indv])
        n_observed_cats <- length(cat_counts)
      } # end while
      total_regens <- total_regens + count_iter
    } # end for(indv in 1:num_indvs)

    return_vals <- list()

    # Estimation
    for (est_choice in est_choice_list)
    {
      cat("estimate_categ_func_data", replica_idx, "\n")
      timestamps01 <- seq(from = 0.0001, to = 1, length = timeseries_length)
      timeKeeperStart("Xiaoxia")
      categFD_est <- estimate_categ_func_data(est_choice, timestamps01, categ_func_data_list$w_mat)
      ##################### 9/11/2023
      if (save_curves) {
        Z_est_curve[[replica_idx]] <- array(c(categFD_est$z1_est, categFD_est$z2_est), dim = c(timeseries_length, num_indvs, 2))
        p_est_curve[[replica_idx]] <- array(c(categFD_est$p1_est, categFD_est$p2_est, categFD_est$p3_est), dim = c(timeseries_length, num_indvs, 3))
        W_cfd[[replica_idx]] <- categ_func_data_list$w_mat
      }
      ##################### 9/11/2023
      timeKeeperNext()

      if (run_hellinger) {
        # evaluate performance Z and P
        rmse1_temp <- c(by(mse_bw_matrix(Z1, categFD_est$z1_est, timestamps01), true_cluster, mean))
        rmse2_temp <- c(by(mse_bw_matrix(Z2, categFD_est$z2_est, timestamps01), true_cluster, mean))
        rmse[replica_idx, , ] <- rbind(rmse1_temp, rmse2_temp)

        error.p1 <- mse_bw_matrixp(p1, categFD_est$p1_est, timestamps01)
        error.p2 <- mse_bw_matrixp(p2, categFD_est$p2_est, timestamps01)
        error.p3 <- mse_bw_matrixp(p3, categFD_est$p3_est, timestamps01)


        hellinger[replica_idx, , ] <- rbind(
          c(by(error.p1, true_cluster, mean)),
          c(by(error.p2, true_cluster, mean)),
          c(by(error.p3, true_cluster, mean))
        )
      }

      # Record clustering performance. USE only Z
      #**** extract_scores_UNIVFPCA is the FASTEST

      if (run_univfpca || run_kmeans || run_dbscan) {
        mfpca_true <- extract_scores_UNIVFPCA(mZ1 = Z1, mZ2 = Z2, tt = timestamps01, PVE = 0.95)
        # plot(  scores_true$scores[, 1:2])
        timeKeeperStart("univfpca")
        mfpca_est <- extract_scores_UNIVFPCA(mZ1 = categFD_est$z1_est, mZ2 = categFD_est$z2_est, tt = timestamps01, PVE = 0.95)
        timeKeeperNext()
      }

      # KMEANS
      if (run_kmeans) {
        true_kmeans_temp <- kmeans_cluster(data = mfpca_true$scores)$label
        timeKeeperStart("kmeans")
        est_kmeans_temp <- kmeans_cluster(data = mfpca_est$scores)$label
        timeKeeperNext()
      } else {
        true_kmeans_temp <- rep(NA, num_indvs)
        est_kmeans_temp <- rep(NA, num_indvs)
      }

      # FADP
      if (run_fadp) {
        true_fadp_temp <- fadp_cluster(mZ1 = Z1, mZ2 = Z2, tt = timestamps01)$label
        timeKeeperStart("fadp")
        est_fadp_temp <- fadp_cluster(mZ1 = categFD_est$z1_est, mZ2 = categFD_est$z2_est, tt = timestamps01)$label
        timeKeeperNext()
      } else {
        true_fadp_temp <- rep(NA, num_indvs)
        est_fadp_temp <- rep(NA, num_indvs)
      }
      # dbscan
      if (run_dbscan) {
        true_dbscan_temp <- dbscan_cluster(data = mfpca_true$scores, 1)$label
        timeKeeperStart("dbscan")
        est_dbscan_temp <- dbscan_cluster(data = mfpca_est$scores, 1)$label
        timeKeeperNext()
      } else {
        true_dbscan_temp <- rep(NA, num_indvs)
        est_dbscan_temp <- rep(NA, num_indvs)
      }
      ## dbscan cfda
      if (run_cfda) {
        basis_num <- 10
        nCores <- n.cores
        true_dbscan_temp_cfda <- true_cluster_db
        parallel::stopCluster(cl = my.cluster)
        timeKeeperStart("cfd")
        cfd_scores <- cfda_score_function(categ_func_data_list$w_mat, nCores, timestamps01, basis_num)
        my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
        doParallel::registerDoParallel(cl = my.cluster)
        cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
        est_dbscan_temp_cfda <- dbscan_cluster(data = cfd_scores, 1)$label
      } else {
        true_dbscan_temp_cfda <- rep(NA, length(true_cluster_db))
        est_dbscan_temp_cfda <- rep(NA, length(true_cluster_db))
      }
      # record results
      true_dbscan <- cbind(true_dbscan, true_dbscan_temp)
      est_dbscan <- cbind(est_dbscan, est_dbscan_temp)

      true_dbscan_cfda <- cbind(true_dbscan_cfda, true_dbscan_temp_cfda)
      est_dbscan_cfda <- cbind(est_dbscan_cfda, est_dbscan_temp_cfda)

      true_kmeans <- cbind(true_kmeans, true_kmeans_temp)
      est_kmeans <- cbind(est_kmeans, est_kmeans_temp)

      true_fadp <- cbind(true_fadp, true_fadp_temp)
      est_fadp <- cbind(est_fadp, est_fadp_temp)

      cat("Done replica:", replica_idx, "\n")
    } # End for est_choice_list
  } # END of "for(replica_idx in 1:num_replicas)'

  cat("\n replicas done \n")

  if (run_hellinger) {
    # Assess accuracy in the simulation study
    mse_sim <- apply(rmse, c(2, 3), mean)
    hellinger_sim <- apply(hellinger, c(2, 3), mean)
  }

  # Evaluate clustering accuracy (RI, ARI, CPN) for each method across replicas.
  ### truth
  # dbscan
  if (run_dbscan) {
    true_dbscan_ri <- apply(true_dbscan, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ri
    })
    true_dbscan_ari <- apply(true_dbscan, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ari
    })
    true_dbscan_cpn <- apply(true_dbscan, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$cpn
    })
  } else {
    true_dbscan_ri <- rep(NA, num_replicas)
    true_dbscan_ari <- rep(NA, num_replicas)
    true_dbscan_cpn <- rep(NA, num_replicas)
  }
  # dbscan cfda
  if (run_cfda) {
    true_dbscan_ri_cfda <- apply(true_dbscan_cfda, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ri
    })
    true_dbscan_ari_cfda <- apply(true_dbscan_cfda, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ari
    })
    true_dbscan_cpn_cfda <- apply(true_dbscan_cfda, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$cpn
    })
  } else {
    true_dbscan_ri_cfda <- rep(NA, num_replicas)
    true_dbscan_ari_cfda <- rep(NA, num_replicas)
    true_dbscan_cpn_cfda <- rep(NA, num_replicas)
  }

  ### kmeans
  if (run_kmeans) {
    true_kmeans_ri <- apply(true_kmeans, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ri
    })
    true_kmeans_ari <- apply(true_kmeans, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ari
    })
    true_kmeans_cpn <- apply(true_kmeans, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$cpn
    })
  } else {
    true_kmeans_ri <- rep(NA, num_replicas)
    true_kmeans_ari <- rep(NA, num_replicas)
    true_kmeans_cpn <- rep(NA, num_replicas)
  }

  # fadp
  if (run_fadp) {
    true_fadp_ri <- apply(true_fadp, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ri
    })
    true_fadp_ari <- apply(true_fadp, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ari
    })
    true_fadp_cpn <- apply(true_fadp, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$cpn
    })
  } else {
    true_fadp_ri <- rep(NA, num_replicas)
    true_fadp_ari <- rep(NA, num_replicas)
    true_fadp_cpn <- rep(NA, num_replicas)
  }



  ##### estimate results
  # dbscan
  if (run_dbscan) {
    est_dbscan_ri <- apply(est_dbscan, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ri
    })
    est_dbscan_ari <- apply(est_dbscan, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ari
    })
    est_dbscan_cpn <- apply(est_dbscan, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$cpn
    })
  } else {
    est_dbscan_ri <- rep(NA, num_replicas)
    est_dbscan_ari <- rep(NA, num_replicas)
    est_dbscan_cpn <- rep(NA, num_replicas)
  }
  # dbscan cfda
  if (run_cfda) {
    est_dbscan_ri_cfda <- apply(est_dbscan_cfda, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ri
    })
    est_dbscan_ari_cfda <- apply(est_dbscan_cfda, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ari
    })
    est_dbscan_cpn_cfda <- apply(est_dbscan_cfda, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$cpn
    })
  } else {
    est_dbscan_ri_cfda <- rep(NA, num_replicas)
    est_dbscan_ari_cfda <- rep(NA, num_replicas)
    est_dbscan_cpn_cfda <- rep(NA, num_replicas)
  }
  ### kmeans
  if (run_kmeans) {
    est_kmeans_ri <- apply(est_kmeans, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ri
    })
    est_kmeans_ari <- apply(est_kmeans, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ari
    })
    est_kmeans_cpn <- apply(est_kmeans, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$cpn
    })
  } else {
    est_kmeans_ri <- rep(NA, num_replicas)
    est_kmeans_ari <- rep(NA, num_replicas)
    est_kmeans_cpn <- rep(NA, num_replicas)
  }

  # fadp
  if (run_fadp) {
    est_fadp_ri <- apply(est_fadp, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ri
    })
    est_fadp_ari <- apply(est_fadp, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ari
    })
    est_fadp_cpn <- apply(est_fadp, 2, function(cluster) {
      evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$cpn
    })
  } else {
    est_fadp_ri <- rep(NA, num_replicas)
    est_fadp_ari <- rep(NA, num_replicas)
    est_fadp_cpn <- rep(NA, num_replicas)
  }





  print("cluster_table")
  cluster_table_true <- c(
    mean(true_dbscan_ri_cfda),
    mean(true_dbscan_ari_cfda),
    mean(true_dbscan_cpn_cfda),
    mean(true_fadp_ri),
    mean(true_fadp_ari),
    mean(true_fadp_cpn),
    mean(true_kmeans_ri),
    mean(true_kmeans_ari),
    mean(true_kmeans_cpn),
    mean(true_dbscan_ri),
    mean(true_dbscan_ari),
    mean(true_dbscan_cpn)
  )
  names(cluster_table_true) <- c(
    "cfda-db RI", "cfda-db ARI", "cfda-db cpn",
    "fadp RI", "fadp ARI", "fadp cpn",
    "kmeans RI", "kmeans ARI", "kmeans cpn",
    "dbscan RI", "dbscan ARI", "dbscan cpn"
  )

  cluster_table_est <- c(
    mean(est_dbscan_ri_cfda),
    mean(est_dbscan_ari_cfda),
    mean(est_dbscan_cpn_cfda),
    mean(est_fadp_ri),
    mean(est_fadp_ari),
    mean(est_fadp_cpn),
    mean(est_kmeans_ri),
    mean(est_kmeans_ari),
    mean(est_kmeans_cpn),
    mean(est_dbscan_ri),
    mean(est_dbscan_ari),
    mean(est_dbscan_cpn)
  )
  names(cluster_table_est) <- c(
    "cfda-db RI", "cfda-db ARI", "cfda-db cpn",
    "fadp RI", "fadp ARI", "fadp cpn",
    "kmeans RI", "kmeans ARI", "kmeans cpn",
    "dbscan RI", "dbscan ARI", "dbscan cpn"
  )


  cluster_table_est_se <- c(
    sd(est_dbscan_ri_cfda) / sqrt(num_replicas),
    sd(est_dbscan_ari_cfda) / sqrt(num_replicas),
    sd(est_dbscan_cpn_cfda) / sqrt(num_replicas),
    sd(est_fadp_ri) / sqrt(num_replicas),
    sd(est_fadp_ari) / sqrt(num_replicas),
    sd(est_fadp_cpn) / sqrt(num_replicas),
    sd(est_kmeans_ri) / sqrt(num_replicas),
    sd(est_kmeans_ari) / sqrt(num_replicas),
    sd(est_kmeans_cpn) / sqrt(num_replicas),
    sd(est_dbscan_ri) / sqrt(num_replicas),
    sd(est_dbscan_ari) / sqrt(num_replicas),
    sd(est_dbscan_cpn) / sqrt(num_replicas)
  )
  names(cluster_table_est_se) <- c(
    "cfda-db RI", "cfda-db ARI", "cfda-db cpn",
    "fadp RI", "fadp ARI", "fadp cpn",
    "kmeans RI", "kmeans ARI", "kmeans cpn",
    "dbscan RI", "dbscan ARI", "dbscan cpn"
  )
  print("returning")

  save(time_elapsed, total_regens, file = file.path(temp_folder, paste("time_elapsed_", num_indvs, "_", timeseries_length, "_",
    scenario, "_", num_replicas, "_", est_choice, "_", run_hellinger, "_neworder.RData",
    sep = ""
  )))

  est_values <- list(
    "cluster_table_est" = cluster_table_est,
    "cluster_table_est_se" = cluster_table_est_se
  )
  if (save_curves) {
    est_values$Z_est_curves <- Z_est_curve
    est_values$p_est_curves <- p_est_curve
  }

  if (run_hellinger) {
    est_values$mse <- mse_sim
    est_values$hellinger <- hellinger_sim
  }

  # Per-replica data for batch combining
  est_values$est_dbscan_ari_reps <- est_dbscan_ari
  est_values$est_dbscan_ri_reps  <- est_dbscan_ri
  if (run_hellinger) {
    est_values$rmse_reps      <- rmse
    est_values$hellinger_reps <- hellinger
  }

  save(est_values, file = file.path(temp_folder, paste("ClusterSim_", num_indvs, "_", timeseries_length, "_",
    scenario, "_", num_replicas, "_", est_choice, "_", run_hellinger, "_neworder.RData",
    sep = ""
  )))

  return_vals[[est_choice]] <- est_values
  return_vals[["cluster_table_true"]] <- cluster_table_true
  return_vals[["Z_true_curves"]] <- Z_true_curve
  return_vals[["p_true_curves"]] <- p_true_curve
  return_vals[["W_cfd"]] <- W_cfd

  return(return_vals)
}

# ---- Legacy version kept for reference and validation --------------------------------
# The active code now calls generate_categ_func_data() from the catfda package.
# The package version adds input validation and returns $w_mat / $x_array instead of
# $W / $X.  The legacy version is preserved below for cross-checking.
#
GenerateCategFuncData_legacy <- function(prob_curves) {
  curve_count <- length(prob_curves)
  num_indvs <- ncol(prob_curves$p1)
  timeseries_length <- nrow(prob_curves$p1)

  W <- matrix(0, ncol = num_indvs, nrow = timeseries_length)
  X_array <- array(0, c(num_indvs, timeseries_length, curve_count))

  for (indv in c(1:num_indvs)) {
    X <- sapply(
      c(1:timeseries_length),
      function(this_time) {
        rmultinom(
          n = 1, size = 1,
          prob = c(
            prob_curves$p1[this_time, indv],
            prob_curves$p2[this_time, indv],
            prob_curves$p3[this_time, indv]
          )
        )
      }
    )
    W[, indv] <- apply(X, 2, which.max)
    X_array[indv, , ] <- t(X)
  }
  return(list(X = X_array, W = W))
}

#' Validate that generate_categ_func_data() (package) matches GenerateCategFuncData_legacy()
#'
#' Runs both functions with the same seed on synthetic probability curves and
#' checks that the categorical output (W / w_mat) and one-hot array (X / x_array)
#' are element-wise identical.
#'
#' @param seed Integer seed for reproducibility (default 42)
#' @param n_time Number of time points (default 30)
#' @param n_indvs Number of individuals (default 10)
#' @return Invisible TRUE if all checks pass; stops with an error message otherwise.
#'
validate_generate_categ_func_data <- function(seed = 42, n_time = 30, n_indvs = 10) {
  # Build simple probability curves that sum to 1
  set.seed(seed)
  raw <- matrix(runif(n_time * n_indvs * 3), nrow = n_time * n_indvs)
  raw <- raw / rowSums(raw)
  p1 <- matrix(raw[, 1], nrow = n_time, ncol = n_indvs)
  p2 <- matrix(raw[, 2], nrow = n_time, ncol = n_indvs)
  p3 <- matrix(raw[, 3], nrow = n_time, ncol = n_indvs)
  prob_curves <- list(p1 = p1, p2 = p2, p3 = p3)

  # Run legacy
  set.seed(seed + 1)
  old <- GenerateCategFuncData_legacy(prob_curves)

  # Run package version (needs matching prob_curves structure)
  set.seed(seed + 1)
  new <- generate_categ_func_data(prob_curves)

  # Compare W vs w_mat
  if (!identical(old$W, new$w_mat)) {
    stop("VALIDATION FAILED: W (legacy) != w_mat (package). Output categories differ.")
  }

  # Compare X vs x_array
  if (!identical(old$X, new$x_array)) {
    stop("VALIDATION FAILED: X (legacy) != x_array (package). One-hot arrays differ.")
  }

  cat("validate_generate_categ_func_data: PASSED",
      "-- W/w_mat and X/x_array are identical across", n_indvs,
      "individuals x", n_time, "time points.\n")
  invisible(TRUE)
}
# ---- End legacy / validation section ------------------------------------------------

#' Get clustered data
#'
#'
GenerateClusterData <- function(setting, scenario, k, num_indvs, timeseries_length) {
  setting_object <- GetMuAndScore(setting, scenario, k)
  cluster_f <- GenerateClusterDataScenario(num_indvs,
    timeseries_length,
    k,
    mu_1 = setting_object$mu_1,
    mu_2 = setting_object$mu_2,
    score_vals = setting_object$score_vals
  )
  return(cluster_f)
}

#' Get fraction of occurrence of each class for a given scenario
#' @param scenario scenario name as a string "A", "B", "C"
#' @return a vector containing the fractions
#'
GetOccurrenceFractions <- function(scenario) {
  occur_fraction <- switch(scenario,
    "A" = c(0.75, 0.22, 0.03),
    "B" = c(0.5, 0.3, 0.2),
    "C" = c(0.1, 0.6, 0.3)
  )

  return(occur_fraction)
}

#' Get mu_1, mu_2 functions, and score_vals objects for a given context.
#' @param setting setting identified as an integer 1,2,3
#' @param scenario scenario name as a string "A", "B", "C"
#' @param k number of points along the score decay axis
#' @return A list that contains mu_1, mu_2, score_vals
#'
GetMuAndScore <- function(setting, scenario, k) {
  all_score_values <- rep(0, k)

  if (1 == setting) {
    mu_1 <- function(t) -1 + 2 * t + 2 * t^2

    mu_2 <- switch(scenario,
      "A" = function(t) -2.5 + exp(t * 2),
      "B" = function(t) -0.5 + exp(t * 2),
      "C" = function(t) -2.5 + exp(t * 2)
    )
    score_front <- switch(scenario,
      "A" = c(1, 1 / 2, 1 / 4),
      "B" = c(1, 1 / 2, 1 / 4),
      "C" = c(50, 25, 5)
    )
  } else if (2 == setting) {
    mu_1 <- function(t) 4 * t^2 - 1.2

    mu_2 <- function(t) 4 * t^2 - 3.5

    score_front <- c(1, 1 / 2, 1 / 4)
  } else if (3 == setting) {
    mu_1 <- function(t) -2.2 + 4 * t^2

    mu_2 <- function(t) -7 + 6 * t^2

    score_front <- c(1, 1 / 4, 1 / 16)
  }

  for (idx in 1:length(score_front))
  {
    all_score_values[idx] <- score_front[idx]
  }

  return(list("mu_1" = mu_1, "mu_2" = mu_2, "score_vals" = all_score_values))
}


#' Generate cluster data for a given scenario
#' @param num_indvs number of individuals
#' @param timeseries_length length of time-series as an integer
#' @param k  number of eigen(psi) functions
#' @param mu_1 mean function for the first latent curve
#' @param mu_2 mean function for the second latent curve
#' @param score_vals the variance of the principal component scores
#'
GenerateClusterDataScenario <- function(num_indvs,
                                        timeseries_length,
                                        k = 3,
                                        mu_1,
                                        mu_2,
                                        score_vals) {
  timestamps01 <- seq(from = 0.0001, to = 1, length = timeseries_length)

  # noise octaves
  # cat("octave", num_indvs, k, num_indvs * k, "\n")
  scores_standard <- matrix(rnorm(num_indvs * k), ncol = k)
  scores <- scores_standard %*% diag(sqrt(score_vals))

  #
  BIG_mu <- c(mu_1(timestamps01), mu_2(timestamps01))
  BIG_phi <- PsiFunc(k, timestamps01)

  Z <- BIG_phi %*% t(scores) + BIG_mu
  Z1 <- Z[1:timeseries_length, ]
  Z2 <- Z[1:timeseries_length + timeseries_length, ]
  expZ1 <- exp(Z1)
  expZ2 <- exp(Z2)
  denom <- 1 + expZ1 + expZ2
  p1 <- expZ1 / denom
  p2 <- expZ2 / denom
  p3 <- 1 / denom

  # vectorize for future work!!!
  return(list(
    Z1 = Z1, Z2 = Z2,
    p1 = p1, p2 = p2, p3 = p3,
    MEAN = BIG_mu, PHI = BIG_phi, MFPC = scores
  ))
}

#' Psi function
#'
PsiFunc <- function(klen, timestamps01) {
  psi_k1 <- sapply(c(1:klen), function(i) sin((2 * i + 1) * pi * timestamps01))
  psi_k2 <- sapply(c(1:klen), function(i) cos(2 * i * pi * timestamps01))
  return(rbind(psi_k1, psi_k2))
}

