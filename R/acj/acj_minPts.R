library(profvis)

library(elbow)
library(foreach)
library(doParallel)


MOP <- function(count_iter, max_iter, var_alpha=0.5)
{
  pow <- 1 / var_alpha
  return(1 - (count_iter^pow) / (max_iter^pow))
}

MOA <- function(count_iter, max_iter, var_gamma=1, var_delta=0.2)
{
  return(var_delta + count_iter * (var_gamma - var_delta) / max_iter)
}

sum_cohesion <- function(colA, colB)
{
  n <- length(colA)

  center <- c(mean(colA), mean(colB))
  t <- foreach(ind = 1:length(colA), .combine = "+") %do%
  {
    sum((c(colA[ind], colB[ind]) - center)^2)
  }

  # multiplier <- 1 / (n * (n - 1) * 0.5)
  # t <- 0
  # for(ind in 1:n)
  # {
  #   # cat("\nA:", ind, "\n")
  #   for(otr in ind:n)
  #   {
  #     if(ind != otr)
  #     {
  #       # cat("\tB:", otr)
  #       t <- t + multiplier * sqrt((colA[ind] - colA[otr])^2 + (colB[ind] - colB[otr])^2)
  #     }
  #   }
  # }
  # # print(t)
  return(t)
}

FitnessFunction <- function(i, eps, minpts, datapoints)
{
  # cat("FF: ", i, eps, minpts, "\n")
  dbscan_results <- dbscan(datapoints, eps = eps, MinPts = floor(minpts))
  # print(dbscan_results$cluster)
  df <- data.frame(cbind(datapoints, X3=dbscan_results$cluster))
  # print("--DF--")
  # print(df)
  # print("===")
  cluster_value_counts <- table(df$X3)
  # print(cluster_value_counts)
  centers <- df |>
    summarize(c1 = mean(x), c2 = mean(y), ss=sum_cohesion(x,y), .by = X3)
  # print(centers)
  cohesion <- sum(centers$ss)
  # print(cohesion)
  core <- c(mean(centers$c1), mean(centers$c2))
  # cat("Core: ", core, "\n")
  separation <- 0
  for(cluster_id in unique(df$X3))
  {
    cluster_id_str <- toString(cluster_id)
    # cat(cluster_id, cluster_value_counts[[cluster_id_str]], "\n")
    # print(centers[centers$X3==cluster_id_str,])
    cluster_size <- cluster_value_counts[[cluster_id_str]]
    center_vals <- centers[centers$X3==cluster_id_str,]
    this_separation <- cluster_size * sum((c(center_vals$c1, center_vals$c2) - core)^2)
    separation <- separation + this_separation
  }
  # print(separation)
  # cat("Answer: ", separation + cohesion, "\n")
  return(cohesion - separation)
}

OBLAOA <- function(X, lb, ub, datapoints, max_count_iter=100, var_mu=0.5)
{
  eps = 1
  count_iter <- 0

  X_previous <- X
  print(X_previous)
  print(dim(X)[1])
  print(X_previous[1,1])
  print(X_previous[1,2])
  print(datapoints)

  fitness_X_previous <- foreach(i = 1:dim(X_previous)[1], .combine=rbind) %do% {
    FitnessFunction(i, X_previous[i,1], X_previous[i,2], datapoints)
  }
  # print(fitness_X_previous)
  # print("done")
  best_fitness <- min(fitness_X_previous)
  idx_x_best <- which(fitness_X_previous == best_fitness)[1]
  x_best <- X_previous[idx_x_best,]
  cat("best", x_best, best_fitness, idx_x_best,"\n")
  print(X_previous)

  X_current <- X_previous  # current is a copy of previous at the beginning

  while(count_iter < max_count_iter)
  {
    moa <- MOA(count_iter, max_count_iter)
    mop <- MOP(count_iter, max_count_iter)

    # cat("top: ", X_current, "\n")

    for(i in 1:dim(X)[1])
    {
      for(j in 1:dim(X)[2])
      {
        r_vals = runif(3)
        if(moa < r_vals[1])
        {
          if(r_vals[2] < 0.5)
          {
            # division
            X_current[i,j] <- x_best[j] / (mop + eps) * ((ub[j] - lb[j]) * var_mu + lb[j])
          }else{
            # multiplication
            X_current[i,j] <- x_best[j] * mop * ((ub[j] - lb[j]) * var_mu + lb[j])
          }
        }else{
          if(r_vals[3] < 0.5)
          {
            # subtraction
            X_current[i,j] <- x_best[j] - mop * ((ub[j] - lb[j]) * var_mu + lb[j])
          }else{
            # addition
            X_current[i,j] <- x_best[j] + mop * ((ub[j] - lb[j]) * var_mu + lb[j])
          }# end if r3
        }# end if r1
      }# end for j
    }# end for i

    # print(X_current)

    # Fitness of all candidates
    fitness_X_current <- foreach(i = 1:dim(X_current)[1], .combine=rbind) %do% {
      FitnessFunction(i, X_current[i,1], X_current[i,2], datapoints)
    }
    # best of the current candidates
    best_fitness <- min(fitness_X_current)
    idx_x_best <- which(fitness_X_current == best_fitness)[1]
    x_best <- X_current[idx_x_best,]

    # Fitness of opposition candidate
    x_opposition <- ub + lb - x_best
    oppose_fitness <- FitnessFunction(-1, x_opposition[1], x_opposition[2], datapoints)
    if(oppose_fitness < best_fitness)
    {
      x_best <- x_opposition
      best_fitness <- oppose_fitness
    }

    cat(x_best, " --> ", best_fitness, "\n")

    X_previous <- X_current # now the current is used and old

    count_iter <- count_iter + 1
  }# end while

  return(x_best)
}# end function



# dist <- kNNdist(data, k=3)
# distlbow <- data.frame(sort(dist))
# distlbow$index = 1:dim(data)[2]

set.seed(123)
lower_bounds <- c(0,0)
upper_bounds <- c(40, 50)
candidates <- matrix(runif(20), ncol=2)
candidates <- candidates %*% diag(upper_bounds - lower_bounds)
# data <- matrix(c(1,2,94,95,3,1,2,3,94,95,14,15,78,88,16,15,20,14,77,87), ncol=2)
data <- multishapes[,1:2]
result <- OBLAOA(candidates, lower_bounds, upper_bounds, data, max_count_iter=10)
answer <- dbscan(data, eps = result[1], MinPts = result[2])
print(answer)
print(answer$cluster)
plot(data, col=answer$cluster+1)
