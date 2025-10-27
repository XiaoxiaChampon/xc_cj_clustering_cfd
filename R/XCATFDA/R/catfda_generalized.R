library(refund)



cluster_simulation <- function(categ_func_data_list_W) {
  # what is Q vals ? better name???
  Q_vals <- unique(c(categ_func_data_list_W))
  if (is.numeric(Q_vals)) {
    Q_vals <- sort(Q_vals)
  }

  # I need to know what this loop is meant to do !??? maybe there is a better way
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

    # better names for the following variables ???
    # 1. check weather one category only appears 1 time and is it in the end of the timeseries
    # 2. OR is it appearing only one time in the begining
    # 3. OR if the category is less thant Q total category
    # In general, one category only occurs 2 times
    # If timepoints=300, one category only occurs less than 4 times 3/300=0.01
    # If timepoints=750, one category only occurs less than 6 times 5/750=0.0067
    tolcat <- table(categ_func_data_list_W[, indv])
    catorder <- order(tolcat, decreasing = TRUE)
    numcat <- length(catorder)
    refcat <- catorder[numcat]
    count_iter <- 0
    while (count_iter < 100 &&
      ((numcat < length(Q_vals)) ||
        (timeseries_length == 300 && min(as.numeric(tolcat)) < 4) ||
        (timeseries_length == 750 && min(as.numeric(tolcat)) < 10)
      )
    ) {
      count_iter <- count_iter + 1

      new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)

      new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
      new_categ_func_data_list <- GenerateCategFuncData(new_prob_curves)

      # what is this 3 ?? arbitrarily chosen?
      categ_func_data_list_W[, indv] <- new_categ_func_data_list_W[, 3]
      Z1[, indv] <- new_cluster_data$Z1[, 3] # latent curves Z1 and Z2
      categ_func_data_list$X[indv, , ] <- 0
      Z2[, indv] <- new_cluster_data$Z2[, 3]

      for (this_time in 1:timeseries_length)
      {
        categ_func_data_list$X[indv, this_time, which(Q_vals == categ_func_data_list_W[, indv][this_time])] <- 1
      }

      tolcat <- table(categ_func_data_list_W[, indv])
      catorder <- order(tolcat, decreasing = TRUE)
      numcat <- length(catorder)
      refcat <- catorder[numcat]
    } # end while
    total_regens <- total_regens + count_iter
  } # end for(indv in 1:num_indvs)

  return_vals <- list()
}





ClusterSimulation <- function(num_indvs, timeseries_length,
                              scenario, num_replicas, est_choice_list, run_hellinger, temp_folder) {
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

    # Recover the latent Gaussian process --> is this always 2 ???
    Z1 <- cbind(cluster_f1$Z1, cluster_f2$Z1, cluster_f3$Z1)
    Z2 <- cbind(cluster_f1$Z2, cluster_f2$Z2, cluster_f3$Z2)

    # Recover the true probability curves --> could there be more than 3 ???
    p1 <- cbind(cluster_f1$p1, cluster_f2$p1, cluster_f3$p1)
    p2 <- cbind(cluster_f1$p2, cluster_f2$p2, cluster_f3$p2)
    p3 <- cbind(cluster_f1$p3, cluster_f2$p3, cluster_f3$p3)
    prob_curves <- list(p1 = p1, p2 = p2, p3 = p3)
    ######### 9/11/2023
    Z_true_curve[[replica_idx]] <- array(c(Z1, Z2), dim = c(timeseries_length, num_indvs, 2))
    p_true_curve[[replica_idx]] <- array(c(p1, p2, p3), dim = c(timeseries_length, num_indvs, 3))
    ############


    # generate categFuncData
    # set.seed(seed_cfd + 100 * replica_idx)
    # cat("\nCategFD", replica_idx, " --> seed: ", seed_cfd + 100 * replica_idx, "\n")
    cat("CategFD", replica_idx, "\n")

    categ_func_data_list <- GenerateCategFuncData(prob_curves)

    # what is Q vals ? better name???
    Q_vals <- unique(c(categ_func_data_list$W))
    if (is.numeric(Q_vals)) {
      Q_vals <- sort(Q_vals)
    }

    # I need to know what this loop is meant to do !??? maybe there is a better way
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

      # better names for the following variables ???
      # 1. check weather one category only appears 1 time and is it in the end of the timeseries
      # 2. OR is it appearing only one time in the begining
      # 3. OR if the category is less thant Q total category
      # In general, one category only occurs 2 times
      # If timepoints=300, one category only occurs less than 4 times 3/300=0.01
      # If timepoints=750, one category only occurs less than 6 times 5/750=0.0067
      tolcat <- table(categ_func_data_list$W[, indv])
      catorder <- order(tolcat, decreasing = TRUE)
      numcat <- length(catorder)
      refcat <- catorder[numcat]
      count_iter <- 0
      while (count_iter < 100 &&
        ((numcat < length(Q_vals)) ||
          (timeseries_length == 300 && min(as.numeric(tolcat)) < 4) ||
          (timeseries_length == 750 && min(as.numeric(tolcat)) < 10)
        )
      ) {
        count_iter <- count_iter + 1

        new_cluster_data <- GenerateClusterData(setting_choice, scenario, 3, 5, timeseries_length)

        new_prob_curves <- list(p1 = new_cluster_data$p1, p2 = new_cluster_data$p2, p3 = new_cluster_data$p3)
        new_categ_func_data_list <- GenerateCategFuncData(new_prob_curves)

        # what is this 3 ?? arbitrarily chosen?
        categ_func_data_list$W[, indv] <- new_categ_func_data_list$W[, 3]
        Z1[, indv] <- new_cluster_data$Z1[, 3] # latent curves Z1 and Z2
        categ_func_data_list$X[indv, , ] <- 0
        Z2[, indv] <- new_cluster_data$Z2[, 3]

        for (this_time in 1:timeseries_length)
        {
          categ_func_data_list$X[indv, this_time, which(Q_vals == categ_func_data_list$W[, indv][this_time])] <- 1
        }

        tolcat <- table(categ_func_data_list$W[, indv])
        catorder <- order(tolcat, decreasing = TRUE)
        numcat <- length(catorder)
        refcat <- catorder[numcat]
      } # end while
      total_regens <- total_regens + count_iter
    } # end for(indv in 1:num_indvs)

    return_vals <- list()

    # Estimation
    for (est_choice in est_choice_list)
    {
      cat("EstimateCategFuncData", replica_idx, "\n")
      timestamps01 <- seq(from = 0.0001, to = 1, length = timeseries_length)
      timeKeeperStart("Xiaoxia")
      categFD_est <- EstimateCategFuncData(est_choice, timestamps01, categ_func_data_list$W)
      ##################### 9/11/2023
      Z_est_curve[[replica_idx]] <- array(c(categFD_est$Z1_est, categFD_est$Z2_est), dim = c(timeseries_length, num_indvs, 2))
      p_est_curve[[replica_idx]] <- array(c(categFD_est$p1_est, categFD_est$p2_est, categFD_est$p3_est), dim = c(timeseries_length, num_indvs, 3))
      W_cfd[[replica_idx]] <- categ_func_data_list$W
      ##################### 9/11/2023
      timeKeeperNext()

      if (run_hellinger) {
        # evaluate performance Z and P
        rmse1_temp <- c(by(mse_bw_matrix(Z1, categFD_est$Z1_est, timestamps01), true_cluster, mean))
        rmse2_temp <- c(by(mse_bw_matrix(Z2, categFD_est$Z2_est, timestamps01), true_cluster, mean))
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

      mfpca_true <- extract_scores_UNIVFPCA(mZ1 = Z1, mZ2 = Z2, tt = timestamps01, PVE = 0.95)
      # plot(  scores_true$scores[, 1:2])
      timeKeeperStart("univfpca")
      mfpca_est <- extract_scores_UNIVFPCA(mZ1 = categFD_est$Z1_est, mZ2 = categFD_est$Z2_est, tt = timestamps01, PVE = 0.95)
      timeKeeperNext()

      # KMEANS
      true_kmeans_temp <- kmeans_cluster(data = mfpca_true$scores)$label
      timeKeeperStart("kmeans")
      est_kmeans_temp <- kmeans_cluster(data = mfpca_est$scores)$label
      timeKeeperNext()

      # FADP
      true_fadp_temp <- fadp_cluster(mZ1 = Z1, mZ2 = Z2, tt = timestamps01)$label
      timeKeeperStart("fadp")
      est_fadp_temp <- fadp_cluster(mZ1 = categFD_est$Z1_est, mZ2 = categFD_est$Z2_est, tt = timestamps01)$label
      timeKeeperNext()
      # dbscan
      true_dbscan_temp <- dbscan_cluster(data = mfpca_true$scores, 1)$label
      timeKeeperStart("dbscan")
      est_dbscan_temp <- dbscan_cluster(data = mfpca_est$scores, 1)$label
      timeKeeperNext()
      ## dbscan cfda
      ########## no cfda
      basis_num <- 10
      nCores <- n.cores
      true_dbscan_temp_cfda <- true_cluster_db
      parallel::stopCluster(cl = my.cluster)
      timeKeeperStart("cfd")
      cfd_scores <- cfda_score_function(categ_func_data_list$W, nCores, timestamps01, basis_num)
      my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
      doParallel::registerDoParallel(cl = my.cluster)
      cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
      est_dbscan_temp_cfda <- dbscan_cluster(data = cfd_scores, 1)$label
      ########## no cfda
      timeKeeperNext()
      # record results
      true_dbscan <- cbind(true_dbscan, true_dbscan_temp)
      est_dbscan <- cbind(est_dbscan, est_dbscan_temp)

      true_dbscan_cfda <- cbind(true_dbscan, true_dbscan_temp_cfda)
      est_dbscan_cfda <- cbind(est_dbscan, est_dbscan_temp_cfda)

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

  ### Code could be simplified into a loop???
  ### truth
  # dbscan
  true_dbscan_ri <- apply(true_dbscan, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ri
  })
  true_dbscan_ari <- apply(true_dbscan, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ari
  })
  true_dbscan_cpn <- apply(true_dbscan, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$cpn
  })
  # dbscan cfda
  true_dbscan_ri_cfda <- apply(true_dbscan_cfda, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ri
  })
  true_dbscan_ari_cfda <- apply(true_dbscan_cfda, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ari
  })
  true_dbscan_cpn_cfda <- apply(true_dbscan_cfda, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$cpn
  })

  ### kmeans
  true_kmeans_ri <- apply(true_kmeans, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ri
  })
  true_kmeans_ari <- apply(true_kmeans, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ari
  })
  true_kmeans_cpn <- apply(true_kmeans, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$cpn
  })

  # fadp
  true_fadp_ri <- apply(true_fadp, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ri
  })
  true_fadp_ari <- apply(true_fadp, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ari
  })
  true_fadp_cpn <- apply(true_fadp, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$cpn
  })



  ##### estimate results
  # dbscan
  est_dbscan_ri <- apply(est_dbscan, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ri
  })
  est_dbscan_ari <- apply(est_dbscan, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ari
  })
  est_dbscan_cpn <- apply(est_dbscan, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$cpn
  })
  # dbscan cfda
  est_dbscan_ri_cfda <- apply(est_dbscan_cfda, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ri
  })
  est_dbscan_ari_cfda <- apply(est_dbscan_cfda, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$ari
  })
  est_dbscan_cpn_cfda <- apply(est_dbscan_cfda, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster_db, new_cluster = cluster, 0)$cpn
  })
  ### kmeans
  est_kmeans_ri <- apply(est_kmeans, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ri
  })
  est_kmeans_ari <- apply(est_kmeans, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ari
  })
  est_kmeans_cpn <- apply(est_kmeans, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$cpn
  })

  # fadp
  est_fadp_ri <- apply(est_fadp, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ri
  })
  est_fadp_ari <- apply(est_fadp, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$ari
  })
  est_fadp_cpn <- apply(est_fadp, 2, function(cluster) {
    evaluate_cluster(true_cluster = true_cluster, new_cluster = cluster, 3)$cpn
  })





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
    scenario, "_", num_replicas, "_", est_choice, "_", run_hellinger, ".RData",
    sep = ""
  )))

  est_values <- list(
    "cluster_table_est" = cluster_table_est,
    "cluster_table_est_se" = cluster_table_est_se,
    "Z_est_curves" = Z_est_curve,
    "p_est_curves" = p_est_curve
  )

  if (run_hellinger) {
    est_values$mse <- mse_sim
    est_values$hellinger <- hellinger_sim

    return_vals$mse <- mse_sim
    return_vals$hellinger <- hellinger_sim
  }

  save(est_values, file = file.path(temp_folder, paste("ClusterSim_", num_indvs, "_", timeseries_length, "_",
    scenario, "_", num_replicas, "_", est_choice, "_", run_hellinger, ".RData",
    sep = ""
  )))

  return_vals[est_choice] <- est_values
  return_vals["cluster_table_true"] <- cluster_table_true
  return_vals["cluster_table_est"] <- cluster_table_est
  return_vals["cluster_table_est_se"] <- cluster_table_est_se
  return_vals["Z_est_curve"] <- Z_est_curve
  return_vals["p_est_curve"] <- p_est_curve
  return_vals["Z_true_curves"] <- Z_true_curve
  return_vals["p_true_curves"] <- p_true_curve
  return_vals["W_cfd"] <- W_cfd

  return(return_vals)
}





#' Generate cluster data for a given scenario
#' @param num_indvs number of individuals
#' @param timeseries_length length of time-series as an integer
#' @param k  number of eigen(psi) functions
#' @param mu_1 mean function for the first latent curve. A function of time t.
#' @param mu_2 mean function for the second latent curve. A function of time t.
#' @param score_vals the variance of the principal component scores. An vector like c(1.0, 0.25, 0.625)
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


generate_cluster_data_scenario <- function(num_indvs,
                                           timeseries_length,
                                           k = 3,
                                           mu_list,
                                           score_vals) {
  num_cats <- length(mu_list)

  if (length(score_vals) != k) {
    stop("Length of score_vals must match number of components k.")
  }

  timestamps01 <- seq(from = 0.0001, to = 1, length.out = timeseries_length)

  # Generate multivariate PC scores
  scores_standard <- matrix(rnorm(num_indvs * k), ncol = k)
  scores <- scores_standard %*% diag(sqrt(score_vals))

  # Construct mean vector (stacked)
  big_mu <- unlist(lapply(mu_list, function(mu_fn) mu_fn(timestamps01)))

  # Generate phi matrix for all categories
}






#*** extract scores using UNIVARIAFE fpca, FAST COV estimation.  SUPER FAST
# INPUT
# latentProcesses: a list of matrices, each of shape m x n (latent processes
#                  for each category)
# tt: grid of length m (same across all processes)
# PVE: proportion of variance explained (default = 0.95)
# OUTPUT
# list with estimated PHI (m by K) with function norm 1
#           projection of data onto PHI (called scores - nby K matrix)
#

extract_scores_univfpca <- function(latent_proc_list, tt, pve_threshold = 0.95) {
  num_latent_procs <- length(latent_proc_list) # L
  num_timesteps <- nrow(latent_proc_list[[1]]) # m
  num_values <- ncol(latent_proc_list[[1]]) # n
  sqrt_num_timesteps <- sqrt(num_timesteps) # sqrt(m)

  # Pass 1: Calculate fpca for each latent process and
  #           store the values to be used in the pass 2.
  scores_list <- vector("list", num_latent_procs)
  eigen_funcs_list <- vector("list", num_latent_procs)
  npc_vec <- integer(num_latent_procs)
  total_npcs <- 0
  for (lpidx in 1:num_latent_procs) {
    fpca_i <- fpca.face(Y = t(latent_proc_list[[lpidx]]), argvals = tt, pve = 0.99)
    scores_list[[lpidx]] <- fpca_i$scores
    eigen_funcs_list[[lpidx]] <- fpca_i$efunctions
    npc_vec[lpidx] <- fpca_i$npc
    total_npcs <- total_npcs + fpca_i$npc
  }

  # Pass 2: Create the actual scores output and estimate Phi
  phi_est0 <- matrix(0, nrow = num_timesteps * num_latent_procs, ncol = total_npcs)

  ef_col_start <- 1
  for (lpidx in 1:num_latent_procs) {
    npci <- npc_vec[lpidx]

    ef_col_end <- ef_col_start + npci - 1
    ef_row_start <- 1 + num_timesteps * (lpidx - 1)
    ef_row_end <- num_timesteps * lpidx

    phi_est0[ef_row_start:ef_row_end, ef_col_start:ef_col_end] <- eigen_funcs_list[[lpidx]] * sqrt_num_timesteps

    # calc col start for next iteration
    ef_col_start <- ef_col_end + 1
  }
  scores0 <- do.call(cbind, scores_list)

  out_eigen0 <- eigen(cov(scores0)) # ,symmetric = TRUE, only.values = FALSE)
  cum_var <- cumsum(out_eigen0$values) / sum(out_eigen0$values)

  last_eigen_idx <- which(cum_var >= pve_threshold)[1]
  count_iter <- 0
  delta <- 0.01
  while (last_eigen_idx < 2 && count_iter < 100) {
    count_iter <- count_iter + 1
    cat("count_iter: ", count_iter, "\n")
    last_eigen_idx <- which(cum_var >= (pve_threshold + delta))[1]
    delta <- delta + 0.01
  }

  phi_est <- phi_est0 %*% out_eigen0$vectors[, 1:last_eigen_idx]

  mZ <- do.call(rbind, latent_proc_list) # shape: (L * m) x n

  scores_est <- crossprod(mZ, phi_est) / sqrt_num_timesteps

  return(list(scores = scores_est, Phi = phi_est))
}
