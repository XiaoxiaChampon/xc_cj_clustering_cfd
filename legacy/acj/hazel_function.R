############################################################
# Copyright 2024 Xiaoxia Champon

# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the “Software”), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
############################################################
# Purpose: hazel functions for real data, including elbow, cfda package
# Author:  Xiaoxia Champon
# Date: 11/26/2024
##############################################################

#################################
elbow <- function(data, plot = FALSE) {
  
  
  ## Argument checks ----
  
  if (missing(data)) {
    stop("Please provide a two-columns data frame.")
  }
  
  if (!is.list(data)) {
    stop("`data` must be a two-columns data frame.")
  }
  
  if (ncol(data) > 2) {
    warning("Only the first two columns will be considered.")
  }
  
  if (!is.numeric(data[ , 1]) || !is.numeric(data[ , 2])) {
    stop("Non-numeric data detected.")
  }
  
  if (sum(is.na(data[ , 1])) + sum(is.na(data[ , 2]))) {
    stop("Missing values detected.")
  }
  
  if (!is.logical(plot)) {
    stop("`plot` must be a boolean.")
  }
  
  
  ## Data transformation ----
  
  data <- data[ , 1:2]
  data <- data[order(data[ , 1]), ]
  
  
  ## Get constant increase/decrease in y ----
  
  constant <- data[c(1, nrow(data)), ]
  colnames(constant) <- c("x", "y")
  
  mod <- stats::lm(y ~ x, data = constant)
  
  data[ , "constant"] <- round(mod$coef[[1]] + mod$coef[[2]] * data[ , 1], 3)
  
  
  ## Detect inflection point ----
  
  pos <- round(nrow(data) / 2)
  
  if (data[pos, "constant"] < data[pos, 2]) { # Concave Down
    
    ymin <- min(data[ , 2])
    data[ , "benefits"] <- ymin + round(data[ , 2] - data[ , "constant"], 3)
    maxi <- data[which.max(data[ , "benefits"]), ]
    
  } else { # Concave Up
    
    ymax <- max(data[ , 2])
    data[ , "benefits"] <- ymax - round(data[ , "constant"] - data[ , 2], 3)
    maxi <- data[which.min(data[ , "benefits"]), ]
  }
  
  
  ## Store results ----
  
  xxx <- list()
  xxx[[1]] <- maxi[1, 1]
  xxx[[2]] <- data
  names(xxx) <- c(paste(colnames(data)[1], "selected", sep = "_"), "data")
  
  
  ## Plot ----
  
  if (plot) {
    
    xlims <- range(data[ , 1])
    ylims <- c(min(c(data[ , 2], data[ , 3], data[ , 4])), max(data[ , 2]))
    
    
    ## Graphical parameters ----
    
    graphics::par(
      mar      = c(2.5, 2.5, 1.5, 1.5),
      family   = "serif",
      cex.axis = 0.85,
      mgp      = c(2, .15, 0),
      tcl      = -0.25,
      fg       = "#666666",
      col      = "#666666",
      col.axis = "#666666"
    )
    
    
    ## Background plot ----
    
    graphics::plot(
      x    = data[ , 1],
      y    = data[ , 2],
      xlim = xlims,
      ylim = ylims,
      ann  = FALSE,
      axes = FALSE,
      type = "n"
    )
    
    graphics::grid()
    graphics::box()
    
    
    ## Add axes ----
    
    graphics::par(mgp = c(2, 0.00, 0))
    graphics::axis(1, lwd = 0)
    graphics::axis(1, maxi[ , 1], lwd = 0, font = 2, col.axis = "black")
    
    graphics::par(mgp = c(2, 0.25, 0))
    graphics::axis(2, lwd = 0, las = 1)
    at <- round(maxi[ , 2], 3)
    graphics::axis(2, at, lwd = 0, font = 2, col.axis = "black", las = 1)
    
    graphics::mtext(side = 1, cex = 1, line = 1.25, text = expression("x"))
    graphics::mtext(side = 2, cex = 1, line = 1.45, text = expression("f(x)"))
    
    
    ## Real gains/losses in y while x increases ----
    
    graphics::polygon(
      x      = c(data[ , 1], data[1, 1]),
      y      = c(data[ , "benefits"], data[1, "benefits"]),
      col    = "#aaaaaa66",
      border = "#aaaaaa"
    )
    
    
    ## Data serie ----
    
    graphics::points(
      x    = data[ , 1],
      y    = data[ , 2],
      type = "b",
      pch  = 19,
      col  = "black"
    )
    
    
    ## Inflection point informations ----
    
    graphics::lines(
      x   = rep(maxi[1, 1], 2),
      y   = c(par()$usr[3], maxi[1, 2]),
      col = "black",
      lwd = 0.5
    )
    
    graphics::lines(
      x   = c(par()$usr[1], maxi[1, 1]),
      y   = rep(maxi[1, 2], 2),
      col = "black",
      lwd = 0.5
    )
    
    graphics::points(
      x    = maxi[ , 1],
      y    = maxi[ , 2],
      type = "b",
      pch  = 19,
      cex  = 1.5,
      col  = "black"
    )
    
    graphics::points(
      x    = maxi[ , 1],
      y    = maxi[ , 4],
      type = "b",
      pch  = 19,
      cex  = 1,
      col  = "#666666"
    )
  }
  
  return(xxx)
}

#############
#generate new Z
extract_scores_UNIVFPCA <- function (mZ1,mZ2, tt , PVE=0.95)
{
  m<- nrow(mZ1)
  n<-ncol(mZ1)
  
  out1 <- refund::fpca.face(Y=t(mZ1), argvals =tt, pve = 0.99)
  out2 <- refund::fpca.face(Y=t(mZ2), argvals =tt, pve = 0.99)
  O_m1 <- matrix(0, nrow=m, ncol=out1$npc)
  O_m2  <- matrix(0, nrow=m, ncol=out2$npc)
  
  #construct PHI
  Phi_est0 <-  rbind(cbind(out1$efunctions*sqrt(m),O_m2 ),
                     cbind(O_m1,out2$efunctions*sqrt(m)))
  Scores0 <- cbind(out1$scores, out2$scores)
  ScoresCov_0 <- cov(Scores0 )
  oute <- eigen(ScoresCov_0)
  
  K<- which(cumsum( oute$values)/sum(oute$values)>=PVE)[1]
  count_iter = 0
  delta=0.01
  while (K<2 && count_iter<100) {
    count_iter = count_iter + 1
    cat("count_iter: ", count_iter, "\n")
    K<- which(cumsum( oute$values)/sum(oute$values)>=(PVE+delta))[1]
    delta=delta+0.01
  }
  
  Phi_est <-  Phi_est0%*% oute$vectors[,1:K] # correct eigenfns
  
  mZ <- rbind(mZ1, mZ2)
  Scores_est <- t(mZ) %*%Phi_est/sqrt(m)  # they are not demeaned
  
  return (list(scores=Scores_est, Phi= Phi_est))
}

#**** clustering the scores using KMEANS
#INPUT
# n by K scores - or features extracted from the multivariate FD
#OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
kmeans_cluster <- function(data=scores_K){
  out_kmeans = NbClust::NbClust(data = data, diss = NULL,
                                distance = "euclidean", min.nc = 2, max.nc = 5,
                                method = "kmeans",index="silhouette")
  
  return(list(nclust=as.numeric(out_kmeans$Best.nc[1]), label = out_kmeans$Best.partition))
}

#**** clustering the scores using KMEANS
#INPUT
# n by K scores - or features extracted from the multivariate FD
# scale_eps, scalar, the scale of the epsilon
#OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
dbscan_cluster <- function(data=scores_K, scale_eps){
  
  dimz=dim(data)[2]
  if (dimz<=2){minPts=4}
  if (dimz>2){minPts=2*dimz+1}
  dist=kNNdist(data, k = minPts-1)
  #########change to max increase
  distdataelbow=data.frame(sort(dist))
  distdataelbow$index=1:(dim(data)[1])
  ipoint <- elbow(data = distdataelbow)
  epsoptimal=(ipoint$sort.dist._selected)*scale_eps
  
  out_dbscan <- dbscan(data, eps =epsoptimal , minPts = minPts)
  return(list(nclust=dim(table(out_dbscan$cluster)), label = out_dbscan$cluster))
}


#*** clustering the curves directly using FADP
#INPUT
# Z1,Z2 - m by n latent processes for the categ FD with 3 categories
# tt - grid of points
# PVE - percentage of explained variance (default value 0.95)
# OUTPUT
# list with 2 elements: nclust (number of clusters), label (vector with cluster membership)
#
fadp_cluster <- function(mZ1, mZ2, tt=tt, PVE=0.95){
  # t_length <- 2000
  # mZ1 <- Z_after_final[1:t_length,]
  # mZ2 <- Z_after_final[2001:(2000+t_length),]
  # tt <- t
  fdabasis = fda::create.bspline.basis(c(0,1), 25, 4)
  fdatime = tt
  fdafd1 = fda::smooth.basis(tt, mZ1, fdabasis)$fd
  fdafd2 = fda::smooth.basis(tt, mZ2, fdabasis)$fd
  
  FADlist=list(fdafd1,fdafd2)
  
  FADP <- FADPclust(fdata = FADlist, cluster = 2:5, method = "FADP1",
                    proportion = seq(0.02, 0.2, 0.02), f.cut = 0.15,
                    pve = PVE, stats = "Avg.silhouette")
  
  return(list(nclust=FADP$nclust, label = FADP$clust))
}
##############################################
#FADPclust package
FADPclust <- function (fdata, cluster = 2:10, method = "FADP1", proportion = NULL, 
                       f.cut = 0.15, pve = 0.9, stats = "Avg.silhouette") 
{
  if (!inherits(fdata, "fd") & !inherits(fdata, "list")) {
    stop("Error in fdata! For univariate FADP: fdata should be a functional data object produced by fd() function of fda package, for multivariate FADP: a list of functional data objects.", 
         sep = "\n")
  }
  if (is.null(proportion) & method == "FADP1") {
    proportion <- seq(0.1, 1, 0.1)
  }
  if (is.null(proportion) & method == "FADP2") {
    proportion <- seq(0.1, 1, 0.1)
  }
  knn_density1 <- function(distance, k) {
    distmat <- as.matrix(distance)
    n <- nrow(distmat)
    v.d <- pi^(1/2)/gamma(1/2 + 1)
    r.k <- apply(distmat, 1, sort)[k + 1, ]
    den <- k/(n * v.d * r.k)
    return(den)
  }
  knn_density2 <- function(score, k) {
    if (is.null(ncol(score)) == FALSE) {
      n <- nrow(score)
      p <- ncol(score)
    }
    else {
      n <- length(score)
      p <- 1
    }
    if (k >= n) {
      stop("k is not a reasonable and valid number!")
    }
    distance <- dist(score, method = "euclidean", upper = TRUE)
    distmat <- as.matrix(distance)
    v.d <- pi^(p/2)/gamma(p/2 + 1)
    r.k <- apply(distmat, 1, sort)[k + 1, ]
    f_hat <- k/(n * v.d * r.k)
    eta.k <- apply(distmat, 1, sort)[k + 1, ]
    mu.k <- mean(eta.k)
    h <- mu.k + sqrt(sum((eta.k - mu.k)^2)/(n - 1))
    den <- c()
    k <- c()
    for (i in 1:nrow(distmat)) {
      k[i] <- sum(as.vector(distmat[i, ]) <= h) - 1
    }
    den <- k/(n * f_hat * (h)^p)
    return(den)
  }
  assignment <- function(center, distance, den) {
    n.curve <- nrow(distance)
    assign <- rep(0, n.curve)
    assign[center] <- 1:length(center)
    assigned.index <- center
    temp.index <- center
    den.max <- which.max(den)
    while (length(which(assign == 0)) != 0) {
      loc <- rep(0, n.curve)
      for (j in setdiff(1:n.curve, assigned.index)) {
        if (j %in% den.max) {
          nearest.neighbor <- center[which.min(distance[center, 
                                                        j])]
        }
        else {
          neighbors.index <- which(den[j] < den)
          nearest.neighbor <- neighbors.index[which.min(distance[j, 
                                                                 neighbors.index])]
        }
        loc[j] <- nearest.neighbor
      }
      for (l in 1:length(temp.index)) {
        loc.index <- which(loc == temp.index[l])
        assign[loc.index] = rep(assign[temp.index[l]], 
                                length(loc.index))
      }
      assigned.index <- which(assign != 0)
      temp.index <- setdiff(assigned.index, temp.index)
    }
    return(assign)
  }
  adp1 <- function(distance, clusters) {
    distance <- as.matrix(distance)
    n.curve <- nrow(distance)
    k.list <- unique(ceiling(proportion * n.curve^(4/5)))
    lapply.adp1 <- function(para.index) {
      k <- paralist[para.index, 1]
      m <- paralist[para.index, 2]
      den <- knn_density1(distance = distance, k = k)
      del <- c()
      den.max <- which(den == max(den), arr.ind = TRUE)
      del <- mapply(function(j) if (j %in% den.max) {
        del[j] <- max(distance[j, ])
      }
      else {
        del[j] <- min(distance[j, which(den[j] < den)])
      }, 1:n.curve)
      if (round(n.curve * f.cut) < m) {
        alpha <- m/n.curve
      }
      else {
        alpha = f.cut
      }
      den.index <- which(den %in% sort(den, decreasing = T)[1:round(n.curve * 
                                                                      alpha)])
      center.temp <- which(del %in% sort(del[den.index], 
                                         decreasing = T)[1:m])
      result.temp <- assignment(center = center.temp, distance = distance, 
                                den = den)
      stats.temp <- cluster.stats(d = distance, clustering = result.temp)
      if (stats == "Avg.silhouette") {
        evaluation <- stats.temp$avg.silwidth
      }
      if (stats == "Dunn") {
        evaluation <- stats.temp$dunn2
      }
      if (stats == "CH") {
        evaluation <- stats.temp$ch
      }
      return(list(evaluation = evaluation, clustering = result.temp, 
                  density = den, delta = del, center = center.temp))
    }
    paralist <- expand.grid(k.list, clusters)
    task.len <- nrow(paralist)
    result.temp <- lapply(1:task.len, lapply.adp1)
    evaluation.list <- unlist(lapply(1:task.len, function(l) result.temp[[l]]$evaluation))
    selected.index <- which.max(evaluation.list)
    k.selected <- paralist[selected.index, 1]
    nclust.selected <- paralist[selected.index, 2]
    clustering.selected <- result.temp[[selected.index]]$clustering
    density.selected <- result.temp[[selected.index]]$density
    delta.selected <- result.temp[[selected.index]]$delta
    center.selected <- result.temp[[selected.index]]$center
    result <- list(nclust.selected, proportion[which(k.list == 
                                                       k.selected)], method, clustering.selected, density.selected, 
                   delta.selected, center.selected, max(evaluation.list))
    names(result) <- c("nclust", "para", "method", "clust", 
                       "density", "delta", "center", stats)
    return(result)
  }
  adp2 <- function(distance, clusters, score) {
    distance <- as.matrix(distance)
    n.curve <- nrow(distance)
    k.list <- unique(ceiling(proportion * n.curve^(4/5)))
    lapply.adp2 <- function(para.index) {
      k <- paralist[para.index, 1]
      m <- paralist[para.index, 2]
      den <- knn_density2(score = score, k = k)
      del <- c()
      den.max <- which(den == max(den), arr.ind = TRUE)
      del <- mapply(function(j) if (j %in% den.max) {
        del[j] <- max(distance[j, ])
      }
      else {
        del[j] <- min(distance[j, which(den[j] < den)])
      }, 1:n.curve)
      if (round(n.curve * f.cut) < m) {
        alpha <- m/n.curve
      }
      else {
        alpha = f.cut
      }
      den.index <- which(den %in% sort(den, decreasing = T)[1:round(n.curve * 
                                                                      alpha)])
      center.temp <- which(del %in% sort(del[den.index], 
                                         decreasing = T)[1:m])
      result.temp <- assignment(center = center.temp, distance = distance, 
                                den = den)
      stats.temp <- cluster.stats(d = distance, clustering = result.temp)
      if (stats == "Avg.silhouette") {
        evaluation <- stats.temp$avg.silwidth
      }
      if (stats == "Dunn") {
        evaluation <- stats.temp$dunn2
      }
      if (stats == "CH") {
        evaluation <- stats.temp$ch
      }
      return(list(evaluation = evaluation, clustering = result.temp, 
                  density = den, delta = del, center = center.temp))
    }
    paralist <- expand.grid(k.list, clusters)
    task.len <- nrow(paralist)
    result.temp <- lapply(1:task.len, lapply.adp2)
    evaluation.list <- unlist(lapply(1:task.len, function(l) result.temp[[l]]$evaluation))
    selected.index <- which.max(evaluation.list)
    k.selected <- paralist[selected.index, 1]
    nclust.selected <- paralist[selected.index, 2]
    clustering.selected <- result.temp[[selected.index]]$clustering
    density.selected <- result.temp[[selected.index]]$density
    delta.selected <- result.temp[[selected.index]]$delta
    center.selected <- result.temp[[selected.index]]$center
    eta.k.opt <- apply(distance, 1, sort)[k.selected + 1, 
    ]
    mu.k.opt <- mean(eta.k.opt)
    h.opt <- mu.k.opt + sqrt(sum((eta.k.opt - mu.k.opt)^2)/(n.curve - 
                                                              1))
    result <- list(nclust.selected, proportion[which(k.list == 
                                                       k.selected)], method, clustering.selected, density.selected, 
                   delta.selected, center.selected, max(evaluation.list))
    names(result) <- c("nclust", "para", "method", "clust", 
                       "density", "delta", "center", stats)
    return(result)
  }
  if (inherits(fdata, "fd") & method == "FADP1") {
    nbasis <- fdata$basis$nbasis
    type <- fdata$basis$type
    distance_L2 <- semimetric.basis(fdata1 = fdata, fdata2 = fdata, 
                                    nderiv = 0, type.basis1 = type, nbasis1 = nbasis, 
                                    type.basis2 = type, nbasis2 = nbasis)
    result <- adp1(distance = distance_L2, clusters = cluster)
  }
  if (inherits(fdata, "list") & method == "FADP1") {
    basis <- fdata[[1]]$basis
    nbasis <- fdata[[1]]$basis$nbasis
    rangeval <- fdata[[1]]$basis$rangeval
    p <- length(fdata)
    n.curve <- ncol(fdata[[1]]$coefs)
    dot <- nbasis * 10
    t <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - 
                                               rangeval[1])/(dot - 1))
    data <- array(0, dim = c(n.curve, dot, p))
    for (i in 1:p) {
      data[, , i] <- t(eval.fd(evalarg = t, fdobj = fdata[[i]]))
    }
    for (i in 1:p) {
      m1 <- min(data[, , i])
      m2 <- max(data[, , i])
      data.temp <- (data[, , i] - m1)/(m2 - m1)
      data[, , i] <- data.temp
    }
    lam <- matrix(0, p, dot)
    a <- array(data = 0, dim = c(p, dot, p))
    for (i in 1:dot) {
      xir <- matrix(0, p, n.curve)
      for (r in 1:p) {
        xir[r, ] <- data[, i, r]
      }
      sigma <- xir %*% t(xir)
      out <- eigen(sigma)
      lam_i <- out[["values"]]
      a_i <- out[["vectors"]]
      lam[, i] <- lam_i
      for (r in 1:p) {
        a[, i, r] <- a_i[, r]
      }
    }
    w <- 5 * (t[2] - t[1])
    l <- c()
    l[1] <- 1
    for (k in 2:dot) {
      if ((t[k] - w) < 1) {
        l[k] <- 1
      }
      else {
        l[k] <- k - 5
      }
    }
    judge <- function(r, a) {
      adj <- c()
      for (k in 2:dot) {
        sum1 <- 0
        sum2 <- 0
        for (j in l[k]:k - 1) {
          sum1 <- sum1 + sqrt(sum((a[, j, r] - a[, k, 
                                                 r])^2))
          sum2 <- sum2 + sqrt(sum((a[, j, r] + a[, k, 
                                                 r])^2))
        }
        if (sum1 > sum2) {
          adj[k] <- 1
        }
      }
      return(adj)
    }
    for (r in 1:p) {
      loc <- which(judge(r, a) == 1)
      a[, loc, r] <- -a[, loc, r]
    }
    v <- apply(lam, 2, sum)
    phi <- c()
    for (k in 1:p) {
      phi_temp <- sum(lam[k, ] * 0.1)/sum(v * 0.1)
      phi[k] <- phi_temp
    }
    m <- 1
    s <- 0
    repeat {
      s <- s + phi[m]
      m <- m + 1
      if (s >= pve) {
        num <- m - 1
        break
      }
    }
    component <- function(r) {
      bind.matrix <- data[, , 1]
      for (m in 2:p) {
        bind.matrix <- cbind(bind.matrix, data[, , m])
      }
      z <- c()
      zir <- matrix(0, n.curve, dot)
      for (i in 1:n.curve) {
        xi <- matrix(bind.matrix[i, ], p, dot, byrow = TRUE)
        for (j in 1:dot) {
          z[j] <- t(a[, j, r]) %*% xi[, j]
        }
        zir[i, ] <- z
      }
      return(zir)
    }
    comp_list <- list()
    beta <- phi[1:num]/sum(phi[1:num])
    distance_list2 <- list()
    distance_L2 <- matrix(0, n.curve, n.curve)
    for (r in 1:num) {
      comp <- component(r)
      compfd <- Data2fd(argvals = t, y = t(comp), basisobj = basis)
      L2_distance <- semimetric.basis(fdata1 = compfd, 
                                      fdata2 = compfd, nderiv = 0, type.basis1 = "bspline", 
                                      nbasis1 = nbasis, type.basis2 = "bspline", nbasis2 = nbasis)
      comp_list[[r]] <- compfd
      distance_list2[[r]] <- L2_distance
      distance_L2 <- distance_L2 + distance_list2[[r]]^2 * 
        beta
    }
    distance <- sqrt(distance_L2)
    result <- adp1(distance = distance, clusters = cluster)
  }
  if (inherits(fdata, "fd") & method == "FADP2") {
    nbasis <- fdata$basis$nbasis
    rangeval <- fdata$basis$rangeval
    dot <- nbasis * 10
    t <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - 
                                               rangeval[1])/(dot - 1))
    fundata <- funData(argvals = t, X = t(eval.fd(evalarg = t, 
                                                  fdobj = fdata)))
    pca <- PACE(fundata, pve = pve)
    score <- pca$scores
    distance_L2 <- dist(score, method = "euclidean", upper = TRUE)
    result <- adp2(distance = distance_L2, clusters = cluster, 
                   score = score)
  }
  if (inherits(fdata, "list") & method == "FADP2") {
    nbasis <- fdata[[1]]$basis$nbasis
    rangeval <- fdata[[1]]$basis$rangeval
    dot <- nbasis * 10
    t <- seq(rangeval[1], rangeval[2], by = (rangeval[2] - 
                                               rangeval[1])/(dot - 1))
    fundata <- list()
    type <- list()
    comp_num <- c()
    for (i in 1:length(fdata)) {
      fundata[[i]] <- funData(argvals = t, X = t(eval.fd(evalarg = t, 
                                                         fdobj = fdata[[i]])))
      type[[i]] <- list(type = "uFPCA")
      comp_num[i] <- PACE(fundata[[i]], pve = pve)$npc
    }
    fundata <- multiFunData(fundata)
    pca <- MFPCA(mFData = fundata, M = max(comp_num), uniExpansions = type)
    score <- pca$scores
    distance_L2 <- dist(score, method = "euclidean", upper = TRUE)
    result <- adp2(distance = distance_L2, clusters = cluster, 
                   score = score)
  }
  return(result)
}

##############################################




#####
#cfda package
cut_data <- function(
    data,
    Tmax,
    prolongLastState = "all",
    NAstate = "Not observed",
    warning = FALSE)
{
  checkData(data)
  checkLogical(warning, "warning")
  if (length(NAstate) > 1) {
    stop("NAstate must have a length of 1")
  }
  if (any(is.na(Tmax)) || !is.numeric(Tmax) || (length(Tmax) != 
                                                1)) {
    stop("Tmax must be a real.")
  }
  d <- do.call(rbind, by(data, data$id, function(x) {
    cut_cfd(x, Tmax, prolongLastState, NAstate, warning)
  }))
  rownames(d) <- NULL
  return(d)
}



estimate_pt <- function (data, NAafterTmax = FALSE) 
{
  checkData(data)
  t_jumps <- sort(unique(data$time))
  uniqueId <- unique(data$id)
  if (is.factor(data$state)) {
    states <- levels(data$state)
  }
  else {
    states <- sort(unique(data$state))
  }
  res <- matrix(0, nrow = length(states), ncol = length(t_jumps), 
                dimnames = list(states, round(t_jumps, 3)))
  for (id in uniqueId) {
    x <- data[data$id == id, ]
    for (i in seq_along(t_jumps)) {
      aux <- id_get_state(x, t_jumps[i], NAafterTmax)
      res[match(aux, states), i] <- res[match(aux, states), 
                                        i] + 1
    }
  }
  res <- prop.table(res, margin = 2)
  out <- list(pt = res, t = t_jumps)
  class(out) <- "pt"
  return(out)
}


compute_optimal_encoding <- function(
    data, basisobj, computeCI = TRUE, nBootstrap = 50, propBootstrap = 1, method = c("precompute", "parallel"),
    verbose = TRUE, nCores = max(1, ceiling(detectCores() / 2))) {
  t1 <- proc.time()
  
  ## check parameters
  check_compute_optimal_encoding_parameters(data, basisobj, nCores, verbose, computeCI, nBootstrap, propBootstrap)
  method <- match.arg(method)
  nCores <- min(max(1, nCores), detectCores() - 1)
  ## end check
  
  # used to determine the moments where the probability is 0 in plot.fmca
  pt <- estimate_pt(data)
  
  if (verbose) {
    cat("######### Compute encoding #########\n")
  }
  
  
  # change state as integer
  out <- stateToInteger(data$state)
  data$state <- out$state
  label <- out$label
  rm(out)
  
  uniqueId <- as.character(unique(data$id))
  nId <- length(uniqueId)
  K <- length(label$label)
  nBasis <- basisobj$nbasis
  
  if (verbose) {
    cat(paste0("Number of individuals: ", nId, "\n"))
    cat(paste0("Number of states: ", K, "\n"))
    cat(paste0("Basis type: ", basisobj$type, "\n"))
    cat(paste0("Number of basis functions: ", nBasis, "\n"))
    cat(paste0("Number of cores: ", nCores, "\n"))
  }
  
  if (method == "precompute") {
    uniqueTime <- sort(unique(data$time))
    
    V <- computeVmatrix2(data, basisobj, K, uniqueId, uniqueTime, nCores, verbose)
    
    Uval <- computeUmatrix2(data, basisobj, K, uniqueId, uniqueTime, nCores, verbose)
  } else {
    V <- computeVmatrix(data, basisobj, K, uniqueId, nCores, verbose)
    
    Uval <- computeUmatrix(data, basisobj, K, uniqueId, nCores, verbose)
  }
  
  fullEncoding <- computeEncoding(Uval, V, K, nBasis, uniqueId, label, verbose, manage0 = TRUE)
  
  
  if (computeCI) {
    signRef <- getSignReference(fullEncoding$alpha)
    
    bootEncoding <- computeBootStrapEncoding(Uval, V, K, nBasis, label, nId, propBootstrap, nBootstrap, signRef, verbose)
    if (length(bootEncoding) == 0) {
      warning("All bootstrap samples return an error. Try to change the basis.")
      out <- c(fullEncoding, list(V = V, basisobj = basisobj, label = label, pt = pt))
    } else {
      varAlpha <- computeVarianceAlpha(bootEncoding, K, nBasis)
      out <- c(fullEncoding, list(V = V, basisobj = basisobj, label = label, pt = pt,
                                  bootstrap = bootEncoding, varAlpha = varAlpha))
    }
  } else {
    out <- c(fullEncoding, list(V = V, basisobj = basisobj, label = label, pt = pt))
  }
  
  class(out) <- "fmca"
  t2 <- proc.time()
  
  out$runTime <- as.numeric((t2 - t1)[3])
  if (verbose) {
    cat(paste0("Run Time: ", round(out$runTime, 2), "s\n"))
  }
  
  return(out)
}

check_compute_optimal_encoding_parameters <- function(data, basisobj, nCores, verbose, computeCI, nBootstrap, propBootstrap) {
  checkData(data)
  checkDataBeginTime(data)
  checkDataEndTmax(data)
  checkDataNoDuplicatedTimes(data)
  
  if (!is.basis(basisobj)) {
    stop("basisobj is not a basis object.")
  }
  if (any(is.na(nCores)) || !is.whole.number(nCores) || (nCores < 1)) {
    stop("nCores must be an integer > 0.")
  }
  checkLogical(verbose, "verbose")
  checkLogical(computeCI, "computeCI")
  if (computeCI) {
    if (any(is.na(nBootstrap)) || (length(nBootstrap) != 1) || !is.whole.number(nBootstrap) || (nBootstrap < 1)) {
      stop("nBootstrap must be an integer > 0.")
    }
    if (any(is.na(propBootstrap)) || !is.numeric(propBootstrap) || (length(propBootstrap) != 1) || (propBootstrap > 1) || (
      propBootstrap <= 0)) {
      stop("propBootstrap must be a real between 0 and 1.")
    }
  }
}

# return a matrix with nId rows and nBasis * nState columns
computeVmatrix <- function(data, basisobj, K, uniqueId, nCores, verbose) {
  nBasis <- basisobj$nbasis
  
  phi <- fd(diag(nBasis), basisobj) # basis function as functional data
  
  # declare parallelization
  if (nCores > 1) {
    cl <- makeCluster(nCores)
  } else {
    cl <- NULL
  }
  
  
  if (verbose) {
    cat("---- Compute V matrix:\n")
    pbo <- pboptions(type = "timer", char = "=")
  } else {
    pboptions(type = "none")
  }
  t2 <- proc.time()
  
  
  # we compute V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  V <- do.call(rbind, pblapply(cl = cl, split(data, data$id), compute_Vxi, phi = phi, K = K)[uniqueId])
  rownames(V) <- NULL
  
  t3 <- proc.time()
  
  if (verbose) {
    cat(paste0("\nDONE in ", round((t3 - t2)[3], 2), "s\n"))
  }
  
  # stop parallelization
  if (nCores > 1) {
    stopCluster(cl)
  }
  
  return(V)
}



# compute_Vxi  (Vxi = \int_0^Tphi_i(t)X_t=x dt)
#
# @param x one individual (id, time, state)
# @param phi basis functions (e.g. output of \code{\link{fd}} on a \code{\link{create.bspline.basis}} output)
# @param K number of state
# @param ... parameters for integrate function
#
# @return vector of size K*nBasis: V[(x=1,i=1)],
# V[(x=1,i=2)],..., V[(x=1,i=nBasis)], V[(x=2,i=1)], V[(x=2,i=2)]
#
# @examples
# K <- 4
# Tmax <- 6
# PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_PJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
# d_JK2 <- cut_data(d_JK, Tmax)
#
# m <- 6
# b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
# I <- diag(rep(1, m))
# phi <- fd(I, b)
# compute_Vxi(d_JK2[d_JK2$id == 1, ], phi, K)
#
# @author Cristian Preda
compute_Vxi <- function(x, phi, K) {
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis) # V11, V12,...V1m, V21, V22, V2m etc VK1... VKm
  
  for (u in seq_len(nrow(x) - 1)) {
    state <- x$state[u]
    
    for (j in seq_len(nBasis)) { # j = la base
      ind <- (state - 1) * nBasis + j
      aux[ind] <- aux[ind] +
        integrate(
          function(t) {
            eval.fd(t, phi[j])
          },
          lower = x$time[u], upper = x$time[u + 1],
          stop.on.error = FALSE
        )$value
    }
  }
  
  return(aux)
}


# return a matrix with nId rows and nBasis * nState columns
computeVmatrix2 <- function(data, basisobj, K, uniqueId, uniqueTime, nCores, verbose) {
  nBasis <- basisobj$nbasis
  
  phi <- fd(diag(nBasis), basisobj) # basis function as functional data
  
  t2 <- proc.time()
  
  index <- data.frame(seq_along(uniqueTime), row.names = uniqueTime)
  
  # V_ij = int(0,T){phi_j(t)*1_X(t)=i} dt
  integrals <- compute_integral_V(phi, uniqueTime, verbose)
  # V <- do.call(
  #   rbind,
  #   pblapply(cl = cl, split(data, data$id), fill_V, integrals = integrals, index = index, K = K, nBasis = nBasis)[uniqueId]
  # )
  V <- do.call(
    rbind,
    lapply(split(data, data$id), fill_V, integrals = integrals, index = index, K = K, nBasis = nBasis)[uniqueId]
  )
  rownames(V) <- NULL
  
  # V <- do.call(rbind, pblapply(cl = cl, split(data, data$id), compute_Vxi, phi = phi, K = K)[uniqueId])
  # rownames(V) <- NULL
  
  t3 <- proc.time()
  
  if (verbose) {
    cat(paste0("\nDONE in ", round((t3 - t2)[3], 2), "s\n"))
  }
  
  return(V)
}


compute_integral_V <- function(phi, uniqueTime, verbose) {
  nBasis <- phi$basis$nbasis
  if (verbose) {
    pb <- timerProgressBar(min = 0, max = nBasis, width = 50)
    on.exit(close(pb))
    jj <- 1
  }
  integrals <- list()
  for (i in seq_len(nBasis)) { # TODO parallel
    integrals[[i]] <- rep(0., length(uniqueTime))
    for (ii in seq_len(length(uniqueTime) - 1)) {
      integrals[[i]][ii] <- integrate(
        function(t) {
          eval.fd(t, phi[i])
        },
        lower = uniqueTime[ii], upper = uniqueTime[ii + 1],
        stop.on.error = FALSE
      )$value
    }
    if (verbose) {
      setTimerProgressBar(pb, jj)
      jj <- jj + 1
    }
  }
  
  return(integrals)
}


fill_V <- function(x, integrals, index, K, nBasis) {
  aux <- rep(0, K * nBasis)
  for (u in seq_len(nrow(x) - 1)) {
    state <- x$state[u]
    s <- as.character(x$time[u])
    e <- as.character(x$time[u + 1])
    for (i in seq_len(nBasis)) {
      integral <- sum(integrals[[i]][index[s, ]:(index[e, ] - 1)])
      ind <- (state - 1) * nBasis + i
      aux[ind] <- aux[ind] + integral
    }
  }
  return(aux)
}


# return a matrix with nId rows and nBasis * nState columns
computeUmatrix <- function(data, basisobj, K, uniqueId, nCores, verbose) {
  nBasis <- basisobj$nbasis
  
  phi <- fd(diag(nBasis), basisobj) # basis functions as functional data
  
  # declare parallelization
  if (nCores > 1) {
    cl <- makeCluster(nCores)
  } else {
    cl <- NULL
  }
  
  t3 <- proc.time()
  if (verbose) {
    cat(paste0("---- Compute U matrix:\n"))
    pbo <- pboptions(type = "timer", char = "=")
  } else {
    pbo <- pboptions(type = "none")
  }
  
  Uval <- do.call(rbind, pblapply(cl = cl, split(data, data$id), compute_Uxij, phi = phi, K = K)[uniqueId])
  
  # stop parallelization
  if (nCores > 1) {
    stopCluster(cl)
  }
  
  
  t4 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t4 - t3)[3], 2), "s\n"))
  }
  
  
  return(Uval)
}



# compute Uxij
#
# compute int_0^T phi_i(t) phi_j(t) 1_{X_{t}=x}
#
# @param x one individual (id, time, state)
# @param phi basis functions (e.g. output of \code{\link{fd}} on a \code{\link{create.bspline.basis}} output)
# @param K number of state
# @param ... parameters for integrate function
#
# @return vector of size K*nBasis*nBasis: U[(x=1,i=1),(x=1,j=1)],
# U[(x=1,i=1), (x=1,j=2)],..., U[(x=1,i=1), (x=1,j=nBasis)], U[(x=2,i=1), (x=2,j=1)], U[(x=2,i=1), (x=2,j=2)]
#
# @examples
# K <- 4
# Tmax <- 6
# PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_PJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
# d_JK2 <- cut_data(d_JK, Tmax)
#
# m <- 6
# b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
# I <- diag(rep(1, m))
# phi <- fd(I, b)
# compute_Uxij(d_JK2[d_JK2$id == 1, ], phi, K)
#
# @author Cristian Preda, Quentin Grimonprez
compute_Uxij <- function(x, phi, K) {
  nBasis <- phi$basis$nbasis
  aux <- rep(0, K * nBasis * nBasis)
  
  for (u in seq_len(nrow(x) - 1)) {
    state <- x$state[u]
    for (i in seq_len(nBasis)) {
      for (j in i:nBasis) { # symmetry between i and j
        integral <- integrate(
          function(t) {
            eval.fd(t, phi[i]) * eval.fd(t, phi[j])
          },
          lower = x$time[u], upper = x$time[u + 1],
          stop.on.error = FALSE
        )$value
        
        ind <- (state - 1) * nBasis * nBasis + (i - 1) * nBasis + j
        aux[ind] <- aux[ind] + integral
        
        # when i == j, we are on the diagonal of the matrix, no symmetry to apply
        if (i != j) {
          ind <- (state - 1) * nBasis * nBasis + (j - 1) * nBasis + i
          aux[ind] <- aux[ind] + integral
        }
      }
    }
  }
  
  return(aux)
}


# return a matrix with nId rows and nBasis * nState columns
computeUmatrix2 <- function(data, basisobj, K, uniqueId, uniqueTime, nCores, verbose) {
  nBasis <- basisobj$nbasis
  
  phi <- fd(diag(nBasis), basisobj) # basis function as functional data
  # TODO https://r-coder.com/progress-bar-r/
  
  t3 <- proc.time()
  if (verbose) {
    cat(paste0("---- Compute U matrix:\n"))
  }
  index <- data.frame(seq_along(uniqueTime), row.names = uniqueTime)
  
  integrals <- compute_integral_U(phi, uniqueTime, verbose)
  
  # Uval <- do.call(
  #   rbind,
  #   pblapply(cl = cl, split(data, data$id), fill_U, integrals = integrals, index = index, K = K, nBasis = nBasis)[uniqueId]
  # )
  
  Uval <- do.call(
    rbind,
    lapply(split(data, data$id), fill_U, integrals = integrals, index = index, K = K, nBasis = nBasis)[uniqueId]
  )
  
  t4 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t4 - t3)[3], 2), "s\n"))
  }
  
  return(Uval)
}


compute_integral_U <- function(phi, uniqueTime, verbose) {
  nBasis <- phi$basis$nbasis
  if (verbose) {
    nIter <- nBasis * (1 + nBasis) / 2
    pb <- timerProgressBar(min = 0, max = nIter, width = 50)
    on.exit(close(pb))
    jj <- 1
  }
  
  integrals <- list()
  for (i in seq_len(nBasis)) { # TODO parallel
    integrals[[i]] <- list()
    for (j in i:nBasis) {
      integrals[[i]][[j]] <- rep(0, length(uniqueTime))
      for (ii in seq_len(length(uniqueTime) - 1)) {
        integrals[[i]][[j]][ii] <- integrate(
          function(t) {
            eval.fd(t, phi[i]) * eval.fd(t, phi[j])
          },
          lower = uniqueTime[ii], upper = uniqueTime[ii + 1],
          stop.on.error = FALSE
        )$value
      }
      if (verbose) {
        setTimerProgressBar(pb, jj)
        jj <- jj + 1
      }
    }
  }
  return(integrals)
}

fill_U <- function(x, integrals, index, K, nBasis) {
  aux <- rep(0, K * nBasis * nBasis)
  for (u in seq_len(nrow(x) - 1)) {
    state <- x$state[u]
    s <- as.character(x$time[u])
    e <- as.character(x$time[u + 1])
    for (i in seq_len(nBasis)) {
      for (j in i:nBasis) { # symmetry between i and j
        integral <- sum(integrals[[i]][[j]][index[s, ]:(index[e, ] - 1)])
        ind <- (state - 1) * nBasis * nBasis + (i - 1) * nBasis + j
        aux[ind] <- aux[ind] + integral
        
        # when i == j, we are on the diagonal of the matrix, no symmetry to apply
        if (i != j) {
          ind <- (state - 1) * nBasis * nBasis + (j - 1) * nBasis + i
          aux[ind] <- aux[ind] + integral
        }
      }
    }
  }
  return(aux)
}

# @author Cristian Preda, Quentin Grimonprez
computeEncoding <- function(Uval, V, K, nBasis, uniqueId, label, verbose, manage0 = FALSE) {
  t4 <- proc.time()
  if (verbose) {
    cat(paste0("---- Compute encoding: "))
  }
  
  G <- cov(V)
  
  # create F matrix
  Fval <- colMeans(Uval)
  Fmat <- matrix(0, ncol = K * nBasis, nrow = K * nBasis) # diagonal-block matrix with K blocks of size nBasis*nBasis
  for (i in seq_len(K)) {
    Fmat[((i - 1) * nBasis + 1):(i * nBasis), ((i - 1) * nBasis + 1):(i * nBasis)] <-
      matrix(Fval[((i - 1) * nBasis * nBasis + 1):(i * nBasis * nBasis)], ncol = nBasis, byrow = TRUE)
  }
  
  # save matrices before modifying them
  outMat <- list(F = Fmat, G = G)
  
  
  # manage the case where there is a column full of 0 (non invertible matrix)
  # if TRUE, we remove the 0-columns (and rows) and throw a warning
  # otherwise we throw an error and the process stops
  if (manage0) {
    # column full of 0
    ind0 <- (colSums(Fmat == 0) == nrow(Fmat))
    if (sum(ind0) > 0) {
      warning(paste(
        "The F matrix contains at least one column of 0s.",
        "At least one state is not present in the support of one basis function.",
        "Corresponding coefficients in the alpha output will have a 0 value."
      ))
    }
    
    F05 <- t(mroot(Fmat[!ind0, !ind0])) # F  = t(F05)%*%F05
    G <- G[!ind0, !ind0]
    V <- V[, !ind0]
  } else {
    ind0 <- (colSums(Fmat == 0) == nrow(Fmat))
    
    # res = eigen(solve(F)%*%G)
    F05 <- t(mroot(Fmat)) # F  = t(F05)%*%F05
    
    if (any(dim(F05) != rep(K * nBasis, 2))) {
      cat("\n")
      if (any(colSums(Fmat) == 0)) {
        stop(paste("F matrix is not invertible. In the support of each basis function,",
                   "each state must be present at least once (p(x_t) != 0 for t in the support).",
                   "You can try to change the basis."
        ))
      }
      
      stop("F matrix is not invertible. You can try to change the basis.")
    }
  }
  
  
  invF05 <- solve(F05)
  # res <- eigen(F05 %*% solve(F) %*% G %*% solve(F05))
  res <- eigen(t(invF05) %*% G %*% invF05)
  
  # eigenvectors (they give the coefficients of the m=nBasis encoding, for each eigenvalue) as a list of matrices of size m x K
  # The first matrix is the first eigenvector. This vector (1rst column in res$vectors) contains the coefficients of the
  # encoding for the state 1 in the first m (=nBasis) positions, then for the state 2 in the m next positions and so on
  # until the k-th state.
  # I put this first column as and matrix and the coefficientsare column of length m
  
  # transform the matrix of eigenvectors as a list
  
  # aux1 = split(res$vectors, rep(1:ncol(res$vectors), each = nrow(res$vectors)))
  invF05vec <- invF05 %*% res$vectors
  aux1 <- split(invF05vec, rep(seq_len(ncol(res$vectors)), each = nrow(res$vectors)))
  
  # we create the matrices m x K for each eigenvalue: 1rst column = coefficients for state 1,
  # 2nd col = coefficients for state 2, etc
  
  alpha <- lapply(aux1, function(w) {
    wb <- rep(NA, nBasis * K)
    wb[!ind0] <- w
    
    return(matrix(wb, ncol = K, dimnames = list(NULL, label$label)))
  })
  
  pc <- V %*% invF05vec
  rownames(pc) <- uniqueId
  
  t5 <- proc.time()
  
  if (verbose) {
    cat(paste0("\nDONE in ", round((t5 - t4)[3], 2), "s\n"))
  }
  
  return(list(eigenvalues = res$values, alpha = alpha, pc = pc, F = outMat$F, G = outMat$G))
}


compute_time_spent <- function(data) {
  ## check parameters
  checkData(data)
  ## end check
  
  if (is.factor(data$state)) {
    labels <- levels(data$state)
  } else {
    labels <- sort(unique(data$state))
  }
  
  res <- by(data, data$id, function(x) {
    compute_time_spent_intern(x, labels)
  })
  out <- do.call(rbind, res)
  colnames(out) <- labels
  class(out) <- "timeSpent"
  
  return(out)
}


# How long an individual stays in each state
# @author Cristian Preda
compute_time_spent_intern <- function(data, labels) {
  aux <- rep(0, length(labels))
  for (i in seq_along(labels)) {
    idx <- which(data$state == labels[i])
    for (u in idx) {
      if (u < nrow(data)) {
        aux[i] <- aux[i] + data$time[u + 1] - data$time[u]
      }
    }
  }
  return(aux)
}

#' Boxplot of time spent in each state
#'
#' @param x output of \code{\link{compute_time_spent}} function
#' @param col a vector containing color for each state
#' @param ... extra parameters for \code{geom_boxplot}
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # cut at Tmax = 8
#' d_JK2 <- cut_data(d_JK, Tmax = 8)
#'
#' # compute time spent by each id in each state
#' timeSpent <- compute_time_spent(d_JK2)
#'
#' # plot the result
#' boxplot(timeSpent, col = c("#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"))
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' boxplot(timeSpent, notch = TRUE, outlier.colour = "black") +
#'   coord_flip() +
#'   labs(title = "Time spent in each state")
#' @author Quentin Grimonprez
#' @seealso \link{compute_time_spent}
#' @family Descriptive statistics
#'
#' @export
boxplot.timeSpent <- function(x, col = NULL) {
  df <- data.frame(timeSpent = as.vector(x), state = factor(rep(colnames(x), each = nrow(x)), levels = colnames(x)))
  p <- ggplot(df, aes_string(x = "state", y = "timeSpent", fill = "state")) +
    geom_boxplot(...) +
    labs(x = "State", y = "Time Spent", fill = "State")
  
  if (!is.null(col)) {
    p <- p + scale_fill_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_fill_hue(drop = FALSE)
  } # keep the same color order as plotData
  
  return(p)
}


#' Compute duration of individuals
#'
#' For each individual, compute the duration
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#'
#' @return a vector containing the duration of each trajectories
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#'
#' # compute duration of each individual
#' duration <- compute_duration(d_JK)
#'
#' hist(duration)
#' @seealso \link{hist.duration}
#' @family Descriptive statistics
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_duration <- function(data) {
  ## check parameters
  checkData(data)
  ## end check
  
  out <- tapply(data$time, as.factor(data$id), function(x) diff(range(x)))
  class(out) <- "duration"
  
  return(out)
}

#' Plot the duration
#'
#'
#' @param x output of \code{\link{compute_duration}} function
#' @param breaks number of breaks. If not given, use the Sturges rule
#' @param ... parameters for \code{geom_histogram}
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#'
#' # compute duration of each individual
#' duration <- compute_duration(d_JK)
#'
#' hist(duration)
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' hist(duration) +
#'   labs(title = "Distribution of the duration")
#' @author Quentin Grimonprez
#' @seealso \link{compute_duration}
#' @family Descriptive statistics
#'
#' @export
hist.duration <- function(x, breaks = NULL) {
  # choose the number of breaks using Sturges rule
  if (is.null(breaks)) {
    breaks <- floor(1 + log2(length(x)))
  }
  
  extraParam <- list(...)
  defaultParam <- list(fill = "lightblue", color = "black", bins = breaks)
  param <- c(extraParam, defaultParam[which(!(names(defaultParam) %in% names(extraParam)))])
  
  ggplot(data.frame(duration = as.vector(x)), aes_string(x = "duration")) +
    do.call(geom_histogram, param) +
    labs(x = "Duration", y = "Frequency")
}


#' Extract the state of each individual at a given time
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state.
#' @param t time at which extract the state
#' @param NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax
#' (useful when individuals has different lengths)
#'
#' @return a vector containing the state of each individual at time t
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # get the state of each individual at time t = 6
#' get_state(d_JK, 6)
#'
#'
#' # get the state of each individual at time t = 12 (> Tmax)
#' get_state(d_JK, 12)
#' # if NAafterTmax = TRUE, it will return NA for t > Tmax
#' get_state(d_JK, 12, NAafterTmax = TRUE)
#'
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
get_state <- function(data, t, NAafterTmax = FALSE) {
  ## check parameters
  checkData(data)
  if (any(is.na(t)) || !is.numeric(t) || (length(t) != 1)) {
    stop("t must be a real.")
  }
  ## end check
  
  out <- by(data, data$id, function(x) {
    id_get_state(x, t, NAafterTmax)
  })
  out2 <- as.vector(out)
  names(out2) <- names(out)
  
  return(out2)
}

# return the state at time t
#
# x cfda dataframe
# t time value
# NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax
# @author Cristian Preda, Quentin Grimonprez
id_get_state <- function(x, t, NAafterTmax = FALSE) {
  if (NAafterTmax && (t > x$time[length(x$time)])) {
    return(NA)
  }
  
  aux <- which(x$time <= t)
  return(x$state[aux[length(aux)]])
}



#' Estimate probabilities to be in each state
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#' @param NAafterTmax if TRUE, return NA if t > Tmax otherwise return the state associated with Tmax
#' (useful when individuals has different lengths)
#'
#' @return A list of two elements:
#' \itemize{
#'   \item{t: vector of time}
#'   \item{pt: a matrix with K (= number of states) rows and with \code{length(t)} columns containing the
#' probabilities to be in each state at each time.}
#' }
#'
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' d_JK2 <- cut_data(d_JK, 10)
#'
#' # estimate probabilities
#' estimate_pt(d_JK2)
#' @author Cristian Preda, Quentin Grimonprez
#' @seealso \link{plot.pt}
#' @family Descriptive statistics
#'
#' @export
estimate_pt <- function(data, NAafterTmax = FALSE) {
  ## check parameters
  checkData(data)
  ## end check
  
  t_jumps <- sort(unique(data$time))
  uniqueId <- unique(data$id)
  
  if (is.factor(data$state)) {
    states <- levels(data$state)
  } else {
    states <- sort(unique(data$state))
  }
  
  res <- matrix(0, nrow = length(states), ncol = length(t_jumps), dimnames = list(states, round(t_jumps, 3)))
  
  for (id in uniqueId) {
    x <- data[data$id == id, ]
    for (i in seq_along(t_jumps)) {
      aux <- id_get_state(x, t_jumps[i], NAafterTmax)
      
      res[match(aux, states), i] <- res[match(aux, states), i] + 1
    }
  }
  
  res <- prop.table(res, margin = 2)
  
  out <- list(pt = res, t = t_jumps)
  class(out) <- "pt"
  
  return(out)
}

# Extract probability to be in each state at a given time
#
# @param pt output of \link{estimate_pt} function
# @param t time value at which the probability is required
#
# @return  probability to be in each state at time t
#
#
# @examples
# # Simulate the Jukes-Cantor model of nucleotide replacement
# K <- 4
# PJK <- matrix(1/3, nrow = K, ncol = K) - diag(rep(1/3, K))
# lambda_PJK <- c(1, 1, 1, 1)
# d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#
# d_JK2 <- cut_data(d_JK, 10)
#
# # estimate probabilities
# pt <- estimate_pt(d_JK2)
#
# get_proba(pt, 1.5)
#
#
# @seealso \link{estimate_pt}
#
# @author Quentin Grimonprez
#
# @export
get_proba <- function(pt, t) {
  # if we do not export this function, there is no need to do theses checks
  # if(!inherits(x, "pt"))
  #   stop("pt must be an object of class pt.")
  # if(any(is.na(t)) || !is.numeric(t) || length(t) != 1)
  #   stop("t must be a real.")
  
  # find the index containing the first time value greater than the given time t
  i <- sum(t >= pt$t)
  
  # if i == 0, the given time is lower than any time in pt, we can't estimate probabilities, we return NA
  if (i == 0) {
    p <- rep(NA, nrow(pt$pt))
    names(p) <- rownames(pt$pt)
    
    return(p)
  }
  
  return(pt$pt[, i])
}


#' Plot probabilities
#'
#' Plot the probabilities of each state at each given time
#'
#' @param x output of \code{\link{estimate_pt}}
#' @param col a vector containing color for each state
#' @param ribbon if TRUE, use ribbon to plot probabilities
#' @param ... only if \code{ribbon = TRUE}, parameter \code{addBorder}, if TRUE, add black border to the ribbons.
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' d_JK2 <- cut_data(d_JK, 10)
#'
#' pt <- estimate_pt(d_JK2)
#'
#' plot(pt, ribbon = TRUE)
#' @author Quentin Grimonprez
#' @method plot pt
#' @seealso \link{estimate_pt}
#' @family Descriptive statistics
#'
#' @export
plot.pt <- function(x, col = NULL, ribbon = FALSE) {
  ## check parameters
  checkLogical(ribbon, "ribbon")
  ## end check
  
  if (ribbon) {
    p <- plot_pt_ribbon(x, col)
  } else {
    p <- plot_pt_classic(x, col)
  }
  
  return(p)
}


# plot line
# @author Quentin Grimonprez
plot_pt_classic <- function(pt, col = NULL) {
  plot_data <- data.frame(
    State = factor(rep(rownames(pt$pt), each = ncol(pt$pt)), levels = rownames(pt$pt)),
    proba = as.vector(t(pt$pt)),
    time = rep(pt$t, nrow(pt$pt))
  )
  
  p <- ggplot(plot_data, aes_string(x = "time", y = "proba", group = "State", colour = "State")) +
    geom_line() +
    ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)")
  
  if (!is.null(col)) {
    p <- p + scale_colour_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_fill_hue(drop = FALSE)
  } # keep the same color order as plotData
  
  return(p)
}



# plot probabilities using ribbon
# @author Quentin Grimonprez
plot_pt_ribbon <- function(pt, col = NULL, addBorder = TRUE) {
  ## check parameters
  checkLogical(addBorder, "addBorder")
  ## end check
  
  plot_data <- as.data.frame(t(apply(pt$pt, 2, cumsum)))
  nState <- ncol(plot_data)
  labels <- paste0("state", names(plot_data))
  shortLabels <- factor(names(plot_data), levels = names(plot_data))
  names(plot_data) <- labels
  plot_data$time <- pt$t
  plot_data$state0 <- rep(0, nrow(plot_data))
  labels <- c("state0", labels)
  
  p <- ggplot(plot_data)
  for (i in seq_len(nState)) {
    p <- p + geom_ribbon(aes_string(
      ymin = paste0("`", labels[i], "`"),
      ymax = paste0("`", labels[i + 1], "`"), x = "time",
      fill = shortLabels[i]
    ),
    colour = ifelse(addBorder, "black", NA), alpha = 0.8
    )
  }
  
  if (!is.null(col)) {
    p <- p + scale_fill_manual(values = col, drop = FALSE)
  } else {
    p <- p + scale_fill_hue(drop = FALSE)
  } # keep the same color order as plotData
  
  p <- p + ylim(0, 1) +
    labs(x = "Time", y = "p(t)", title = "P(X(t) = x)", fill = "State")
  
  return(p)
}


#' Compute the number of jumps
#'
#' For each individual, compute the number of jumps performed
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
#' \code{state}, associated state.
#' @param countDuplicated if \code{TRUE}, jumps in the same state are counted as jump
#'
#' @return A vector containing the number of jumps for each individual
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # compute the number of jumps
#' nJump <- compute_number_jumps(d_JK)
#' @seealso \link{hist.njump}
#' @family Descriptive statistics
#' @author Cristian Preda, Quentin Grimonprez
#'
#' @export
compute_number_jumps <- function(data, countDuplicated = FALSE) {
  ## check parameters
  checkData(data)
  checkLogical(countDuplicated, "countDuplicated")
  ## end check
  
  out <- by(data, data$id, function(x) {
    compute_number_jumpsIntern(x, countDuplicated)
  })
  nom <- names(out)
  out <- as.vector(out)
  names(out) <- nom
  class(out) <- "njump"
  
  return(out)
}

# @param state vector with state, ordered by time
# @param countDuplicated if TRUE jump in the same state are counted
# @author Quentin Grimonprez
compute_number_jumpsIntern <- function(x, countDuplicated = TRUE) {
  if (countDuplicated) {
    return(length(x$state) - 1)
  } else {
    out <- rle(as.character(x$state[order(x$time)]))$values # rle does not manage factor, as.character allows it
    return(length(out) - 1)
  }
}


#' Plot the number of jumps
#'
#'
#' @param x output of \code{\link{compute_number_jumps}} function
#' @param breaks number of breaks. If not given, use the Sturges rule
#' @param ... parameters for \code{geom_histogram}
#'
#' @return a \code{ggplot} object that can be modified using \code{ggplot2} package.
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' nJump <- compute_number_jumps(d_JK)
#'
#' hist(nJump)
#'
#' # modify the plot using ggplot2
#' library(ggplot2)
#' hist(nJump, fill = "#984EA3") +
#'   labs(title = "Distribution of the number of jumps")
#' @author Quentin Grimonprez
#' @seealso \link{compute_number_jumps}
#' @family Descriptive statistics
#'
#' @export
hist.njump <- function(x, breaks = NULL) {
  # choose the number of breaks using Sturges rule
  if (is.null(breaks)) {
    breaks <- min(floor(1 + log2(length(x))), max(x) + 1)
  }
  
  extraParam <- list(...)
  defaultParam <- list(fill = "lightblue", color = "black", bins = breaks, center = 0)
  param <- c(extraParam, defaultParam[which(!(names(defaultParam) %in% names(extraParam)))])
  
  ggplot(data.frame(njump = as.vector(x)), aes_string(x = "njump")) +
    do.call(geom_histogram, param) +
    labs(x = "Number of jumps", y = "Frequency") +
    scale_x_continuous(breaks = function(x) pretty(seq(ceiling(x[1]), floor(x[2]), by = 1)))
}


#' Table of transitions
#'
#' Calculates a frequency table counting the number of times each pair of states were observed in successive observation times.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#' @param removeDiagonal if TRUE, does not count transition from a state i to i
#'
#' @return a matrix of size \code{K*K} containing the number of transition for each pair
#'
#' @examples
#' # Simulate the Jukes-Cantor model of nucleotide replacement
#' K <- 4
#' PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
#' lambda_PJK <- c(1, 1, 1, 1)
#' d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)
#'
#' # table of transitions
#' statetable(d_JK)
#' @author Quentin Grimonprez
#' @family Descriptive statistics
#'
#' @export
statetable <- function(data, removeDiagonal = FALSE) {
  ## check parameters
  checkData(data)
  ## end check
  
  newState <- stateToInteger(data$state)
  
  out <- statetable.msm(newState$state, data$id)
  
  # If there is at least 1 absorbing state, the matrix is not a square matrix
  out <- completeStatetable(out)
  
  colnames(out) <- newState$label$label[match(colnames(out), newState$label$code)]
  rownames(out) <- newState$label$label[match(rownames(out), newState$label$code)]
  
  if (removeDiagonal) {
    diag(out) <- 0
  }
  
  return(out)
}


compute_optimal_encoding <- function (data, basisobj, computeCI = TRUE, nBootstrap = 50, 
                                      propBootstrap = 1, method = c("precompute", "parallel"), 
                                      verbose = TRUE, nCores = max(1, ceiling(detectCores()/2))
) 
{
  t1 <- proc.time()
  check_compute_optimal_encoding_parameters(data, basisobj, 
                                            nCores, verbose, computeCI, nBootstrap, propBootstrap)
  method <- match.arg(method)
  nCores <- min(max(1, nCores), detectCores() - 1)
  pt <- estimate_pt(data)
  if (verbose) {
    cat("######### Compute encoding #########\n")
  }
  out <- stateToInteger(data$state)
  data$state <- out$state
  label <- out$label
  rm(out)
  uniqueId <- as.character(unique(data$id))
  nId <- length(uniqueId)
  K <- length(label$label)
  nBasis <- basisobj$nbasis
  if (verbose) {
    cat(paste0("Number of individuals: ", nId, "\n"))
    cat(paste0("Number of states: ", K, "\n"))
    cat(paste0("Basis type: ", basisobj$type, "\n"))
    cat(paste0("Number of basis functions: ", nBasis, "\n"))
    cat(paste0("Number of cores: ", nCores, "\n"))
  }
  if (method == "precompute") {
    uniqueTime <- sort(unique(data$time))
    V <- computeVmatrix2(data, basisobj, K, uniqueId, uniqueTime, 
                         nCores, verbose)
    Uval <- computeUmatrix2(data, basisobj, K, uniqueId, 
                            uniqueTime, nCores, verbose)
  }
  else {
    V <- computeVmatrix(data, basisobj, K, uniqueId, nCores, 
                        verbose)
    Uval <- computeUmatrix(data, basisobj, K, uniqueId, nCores, 
                           verbose)
  }
  fullEncoding <- computeEncoding(Uval, V, K, nBasis, uniqueId, 
                                  label, verbose, manage0 = TRUE)
  if (computeCI) {
    signRef <- getSignReference(fullEncoding$alpha)
    bootEncoding <- computeBootStrapEncoding(Uval, V, K, 
                                             nBasis, label, nId, propBootstrap, nBootstrap, signRef, 
                                             verbose)
    if (length(bootEncoding) == 0) {
      warning("All bootstrap samples return an error. Try to change the basis.")
      out <- c(fullEncoding, list(V = V, basisobj = basisobj, 
                                  label = label, pt = pt))
    }
    else {
      varAlpha <- computeVarianceAlpha(bootEncoding, K, 
                                       nBasis)
      out <- c(fullEncoding, list(V = V, basisobj = basisobj, 
                                  label = label, pt = pt, bootstrap = bootEncoding, 
                                  varAlpha = varAlpha))
    }
  }
  else {
    out <- c(fullEncoding, list(V = V, basisobj = basisobj, 
                                label = label, pt = pt))
  }
  class(out) <- "fmca"
  t2 <- proc.time()
  out$runTime <- as.numeric((t2 - t1)[3])
  if (verbose) {
    cat(paste0("Run Time: ", round(out$runTime, 2), "s\n"))
  }
  return(out)
}

# Check if the data.frame has the required format
# @author Quentin Grimonprez
checkData <- function(data) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }
  
  requiredColNames <- c("id", "time", "state")
  missColNames <- !(requiredColNames %in% colnames(data))
  if (any(missColNames)) {
    stop(paste0("Missing columns in data: ", paste(requiredColNames[missColNames], collapse = ", "), "."))
  }
  
  if (nrow(data) <= 1) {
    stop("There is only one row or less.")
  }
  
  if (any(is.na(data))) {
    stop("There is some missing values.")
  }
  
  invisible(return(NULL))
}


# Check if all individuals end with the same time value
# @author Quentin Grimonprez
checkDataEndTmax <- function(data) {
  lastTime <- tapply(data$time, data$id, function(x) x[length(x)])
  
  nLastTime <- length(unique(lastTime))
  
  if (nLastTime != 1) {
    stop("Each individual must end with the same time value.")
  }
  
  invisible(NULL)
}


# Check if all individuals start with the same time value
# @author Quentin Grimonprez
checkDataBeginTime <- function(data) {
  firstTime <- tapply(data$time, data$id, function(x) x[1])
  
  nFirstTime <- length(unique(firstTime))
  
  if (nFirstTime != 1) {
    stop("Each individual must begin with the same time value.")
  }
  
  invisible(NULL)
}

# Check if all individual has different time values
# @author Quentin Grimonprez
checkDataNoDuplicatedTimes <- function(data) {
  duplicatedTimes <- any(tapply(data$time, data$id, function(x) any(duplicated(x))))
  
  if (duplicatedTimes) {
    warning("Some ids contain duplicated time values.")
  }
  
  invisible(NULL)
}

# Check if the given parameter is a single boolean
# @author Quentin Grimonprez
checkLogical <- function(x, paramName) {
  if (length(x) != 1) {
    stop(paste0(paramName, " must be either TRUE or FALSE."))
  }
  
  if (is.na(x)) {
    stop(paste0(paramName, " must be either TRUE or FALSE."))
  }
  
  if (!is.logical(x)) {
    stop(paste0(paramName, " must be either TRUE or FALSE."))
  }
  
  invisible(NULL)
}

# Check if it is an integer (or vector of integer)
# @author Quentin Grimonprez
is.whole.number <- function(x) {
  x == as.integer(x)
}

######
cut_data <- function(data, Tmax, prolongLastState = "all", NAstate = "Not observed", warning = FALSE) {
  ## check parameters
  checkData(data)
  checkLogical(warning, "warning")
  if (length(NAstate) > 1) {
    stop("NAstate must have a length of 1")
  }
  if (any(is.na(Tmax)) || !is.numeric(Tmax) || (length(Tmax) != 1)) {
    stop("Tmax must be a real.")
  }
  ## end check
  
  d <- do.call(rbind, by(data, data$id, function(x) {
    cut_cfd(x, Tmax, prolongLastState, NAstate, warning)
  }))
  rownames(d) <- NULL
  
  return(d)
}

# @author Cristian Preda
cut_cfd <- function(data, Tmax, prolongLastState = "all", NAstate = NA, warning = FALSE) {
  l <- nrow(data)
  currTmax <- max(data$time)
  
  if (Tmax > currTmax) {
    if (((length(prolongLastState) > 0) && all(prolongLastState == "all")) || (data$state[l] %in% prolongLastState)) {
      return(rbind(data, data.frame(id = data$id[1], state = data$state[l], time = Tmax)))
    } else {
      if (warning) {
        warning(paste0("id ", data$id[1], " does not end with an absorbing state. Cannot impute the state until time ",
                       Tmax, ". Please, add more records or change the Tmax value."))
      }
      d <- data
      if (is.factor(d$state)) {
        levels(d$state) = c(levels(d$state), NAstate)
      }
      d$state[l] <- NAstate
      return(rbind(d, data.frame(id = data$id[1], state = NAstate, time = Tmax)))
    }
  } else {
    if (currTmax == Tmax) {
      return(data)
    } else {
      if (Tmax %in% data$time) {
        k <- which(data$time == Tmax)
        return(data[seq_len(k), ])
      } else {
        k <- max(which(data$time <= Tmax))
        return(rbind(data[seq_len(k), ], data.frame(state = data$state[k], time = Tmax, id = data$id[1])))
      }
    }
  }
}

# change the labels into integer
#
# @param state vector with labels
# @return a list with state containing the new formatted state and label,
# a data.frame containing the labels and corresponding integers
#
# @author Quentin Grimonprez
stateToInteger <- function(state) {
  lab <- data.frame(label = sort(unique(state)), code = seq_along(unique(state)))
  
  newstate <- refactorCategorical(state, lab$label, lab$code)
  
  return(list(state = newstate, label = lab))
}


# Rename a categorical value
#
# @param data matrix/data.frame/vector containing the data
# @param oldCateg vector containing categories to change
# @param newCateg vector containing new categorical values
#
# @return Data with new categorical values
#
# @examples
# dat <- c("single", "married", "married", "divorced", "single")
# refactorCategorical(dat, c("single", "married", "divorced"), 1:3)
#
# @author Quentin Grimonprez
refactorCategorical <- function(data, oldCateg = unique(data), newCateg = seq_along(oldCateg)) {
  ind <- match(data, oldCateg)
  
  if (any(is.na(ind[!is.na(data)]))) {
    warning("NA produced.")
  }
  
  return(newCateg[ind])
}


#' Remove duplicated states
#'
#' Remove duplicated consecutive states from data.
#' If for an individual there is two or more consecutive states that are identical, only the first is kept.
#' Only time when the state changes are kept.
#'
#' @param data data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs
#' and \code{state}, associated state.
#' @param keep.last if TRUE, keep the last state for every individual even if it is a duplicated state.
#'
#' @return \code{data} without duplicated consecutive states
#'
#' @examples
#' data <- data.frame(
#'   id = rep(1:3, c(10, 3, 8)), time = c(1:10, 1:3, 1:8),
#'   state = c(rep(1:5, each = 2), 1:3, rep(1:3, c(1, 6, 1)))
#' )
#' out <- remove_duplicated_states(data)
#' @author Quentin Grimonprez
#' @family format
#' @export
remove_duplicated_states <- function(data, keep.last = TRUE) {
  out <- do.call(rbind, by(data, data$id, remove_duplicated_states.intern, keep.last))
  row.names(out) <- NULL
  
  out
}

# @author Quentin Grimonprez
remove_duplicated_states.intern <- function(data, keep.last = TRUE) {
  data <- data[order(data$time), ]
  
  outRle <- rle(as.character(data$state))
  indToKeep <- 1 + c(0, cumsum(outRle$lengths[-length(outRle$lengths)]))
  
  if (keep.last && indToKeep[length(indToKeep)] != nrow(data)) {
    indToKeep <- c(indToKeep, nrow(data))
  }
  
  return(data[indToKeep, ])
}


#' Convert a matrix to a cfda data.frame
#'
#' @param X matrix containing the states
#' @param times time values. If \code{NULL}, it uses a sequence of integers starting with 1
#' @param labels id labels. If \code{NULL}, it uses the matrix colnames
#' @param byrow if \code{FALSE}, one column = one trajectory
#'
#' @return a data.frame in the cfda format
#'
#' @examples
#' x <- matrix(c("a", "b", "c", "c",
#'               "c", "a", "a", "a",
#'               "b", "c", "a", "b"), ncol = 4, byrow = TRUE,
#'               dimnames = list(NULL, paste0("ind", 1:4)))
#' matrixToCfd(x)
#' @family format
#' @export
matrixToCfd <- function(X, times = NULL, labels = NULL, byrow = FALSE) {
  checkLogical(byrow, "byrow")
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or a data.frame")
  }
  nTimes <- ifelse(byrow, ncol(X), nrow(X))
  nInd <- ifelse(byrow, nrow(X), ncol(X))
  
  if (is.null(times)) {
    times <- seq_len(nTimes)
  } else {
    if (!is.numeric(times) || !((is.vector(times) && (length(times) == nTimes)) ||
                                (is.matrix(times) && (length(times) == nTimes)) ||
                                (is.matrix(times) && (length(times) == nTimes * nInd)))) {
      stop(paste0("times must be a numeric vector of length ", nTimes, " or a matrix of length ", nTimes, "x", nInd))
    }
  }
  
  if (byrow) {
    X <- t(X)
    if (is.vector(times)) {
      times <- matrix(times, ncol = 1)
    } else {
      times <- t(times)
    }
  }
  
  if (is.vector(times)) {
    times <- matrix(times, ncol = 1)
  }
  
  timesPerInd <- ncol(times) > 1
  
  if (is.null(labels)) {
    if (!is.null(colnames(X))) {
      labels <- colnames(X)
    } else {
      labels <- seq_len(ncol(X))
    }
  } else if (!is.vector(labels) || (length(labels) != nInd)) {
    stop(paste0("labels must be a vector of length ", nInd))
  }
  
  outData <- data.frame(id = c(), time = c(), state = c())
  
  for (ind in seq_len(nInd)) {
    indT <- ifelse(timesPerInd, ind, 1)
    outData <- rbind(outData, data.frame(id = labels[ind], time = times[1, indT], state = X[1, ind]))
    if (nTimes > 2) {
      # do not copy 2 consecutive time values with the same state
      for (time in 2:(nTimes - 1)) {
        if (X[time, ind] != X[time - 1, ind]) {
          outData <- rbind(outData, data.frame(id = labels[ind], time = times[time, indT], state = X[time, ind]))
        }
      }
    }
    outData <- rbind(outData, data.frame(id = labels[ind], time = times[nTimes, indT], state = X[nTimes, ind]))
  }
  
  rownames(outData) <- NULL
  return(outData)
}

# convert a fd object to a categorical functional data frame (see convertToCfd)
fdToCfd <- function(fd, breaks, labels = NULL, include.lowest = FALSE, right = TRUE, times = NULL, idLabels = NULL, nx = 200) {
  if (!inherits(fd, "fd")) {
    stop("fd is not a fd object")
  }
  if (is.null(times)) {
    times <- seq(fd$basis$rangeval[1], fd$basis$rangeval[2], length = nx)
  }
  if (is.null(idLabels)) {
    idLabels <- fd$fdnames$reps
  }
  X <- eval.fd(times, fd)
  return(quantiMatrixToCfd(X, breaks, labels = labels, include.lowest = include.lowest, right = right,
                           idLabels = idLabels, times = times, byrow = FALSE))
}

# convert a qualitative matrix to a categorical functional data frame (see convertToCfd)
quantiMatrixToCfd <- function(X, breaks, labels = NULL, include.lowest = FALSE, right = TRUE,
                              times = NULL, idLabels = NULL, byrow = FALSE) {
  X <- matrix(cut(X, breaks = breaks, labels = labels, right = right, include.lowest = include.lowest),
              nrow = nrow(X), dimnames = dimnames(X))
  if (any(is.na(X))) {
    stop("The conversion has generated NA. Please, correct your breaks.")
  }
  return(matrixToCfd(X, times, idLabels, byrow))
}

#' Convert data to categorical functional data
#'
#' @param x matrix or fd object
#' @param breaks either a numeric vector of two or more unique cut points or a single number (greater than or equal to 2)
#' giving the number of intervals into which x is to be cut.
#' @param labels labels for the levels of the resulting category. By default, labels are constructed using "(a,b]"
#' interval notation. If labels = FALSE, simple integer codes are returned instead of a factor.
#' @param include.lowest logical, indicating if an ‘x[i]’ equal to the lowest (or highest, for right = FALSE) ‘breaks’
#' value should be included.
#' @param right logical, indicating if the intervals should be closed on the right (and open on the left) or vice versa.
#' @param times vector containing values at which \code{fd} is to be evaluated
#' @param idLabels vector containing id labels. If NULL it use the names found in the matrix or fd object
#' @param nx Only if \code{x} is a fd object. Number of points to evaluate \code{fd}
#' @param byrow Only if \code{x} is a matrix. If \code{FALSE}, one column = one trajectory
#'
#' @return a data.frame in the cfda format
#'
#' @examples
#' # fd object
#' data("CanadianWeather")
#' temp <- CanadianWeather$dailyAv[,, "Temperature.C"]
#' basis <- create.bspline.basis(c(1, 365), nbasis = 8, norder = 4)
#' fd <- smooth.basis(1:365, temp, basis)$fd
#'
#' # "Very Cold" = [-50:-10), "Cold" = [-10:0)
#' out <- convertToCfd(fd, breaks = c(-50, -10, 0, 10, 20, 50),
#'                     labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"),
#'                     times = 1:365)
#'
#' # matrix
#' out2 <- convertToCfd(temp, breaks = c(-50, -10, 0, 10, 20, 50),
#'                      labels = c("Very Cold", "Cold", "Fresh", "OK", "Hot"),
#'                      times = 1:365, byrow = FALSE)
#' @seealso \link{flours}
#' @family format
#' @export
convertToCfd <- function(x, breaks, labels = NULL, include.lowest = FALSE, right = TRUE, times = NULL,
                         idLabels = NULL, nx = 200, byrow = FALSE) {
  if (inherits(x, "fd")) {
    return(fdToCfd(x, breaks, labels, include.lowest, right, times, idLabels = NULL, nx))
  } else if (is.matrix(x) || is.data.frame(x)) {
    return(quantiMatrixToCfd(x, breaks, labels, include.lowest, right, times, idLabels = NULL, byrow))
  }
}
#######
computeBootStrapEncoding <- function(Uval, V, K, nBasis, label, nId, propBootstrap, nBootstrap, signReference, verbose) {
  outEnc <- list()
  t3 <- proc.time()
  if (verbose) {
    cat("---- Compute Bootstrap Encoding:\n")
  }
  
  for (i in seq_len(nBootstrap)) {
    if (verbose) {
      cat("*")
    }
    idToKeep <- sample(nId, floor(propBootstrap * nId), replace = TRUE)
    
    try({
      outEnc[[i]] <- computeEncoding(Uval[idToKeep, ], V[idToKeep, ], K, nBasis, idToKeep,
                                     label, verbose = FALSE, manage0 = TRUE)
    })
    
    # outEnc[[i]] = c(outEnc[[i]] , list(basisobj = basisobj))
    # class(outEnc[[i]]) = "fmca"
  }
  
  # reorder alpha, pc such that representation have the same sign for each bootstrap sample
  outEnc <- unifySign(outEnc, signReference)
  
  
  t4 <- proc.time()
  if (verbose) {
    cat(paste0("\nDONE in ", round((t4 - t3)[3], 2), "s\n"))
  }
  
  return(outEnc)
}

# get the sign of the alpha of the full encoding to later ensure that bootstrap samples have the same sign
# @author Quentin Grimonprez
getSignReference <- function(alpha) {
  pos <- rep(0, length(alpha))
  isNeg <- rep(FALSE, length(alpha))
  for (i in seq_along(alpha)) {
    pos[i] <- which.max(abs(alpha[[i]]))
    isNeg[i] <- Re(alpha[[i]][pos[i]]) < 0
  }
  
  return(list(position = pos, isNegative = isNeg, allNegative = lapply(alpha, function(x) Re(x) < 0)))
}


# eigenvectors are equivalent up to the sign
# we try to have the same sign for each bootstrap sample
#
# @param out a list of computeEncoding function output
#
# @author Quentin Grimonprez
unifySign <- function(out, signReference) {
  for (i in seq_along(out)) {
    if (!is.null(out[[i]])) { # an output can be NULL due to inversion problem
      # we look if the element at the given position is negative
      for (j in seq_along(out[[i]]$alpha)) {
        signNeg <- Re(out[[i]]$alpha[[j]][signReference$position[j]]) < 0
        
        # if there is a NA at the given position, we try an other position
        if (!is.na(signNeg)) {
          if (signNeg != signReference$isNegative[j]) {
            out[[i]]$alpha[[j]] <- out[[i]]$alpha[[j]] * -1
            out[[i]]$pc[, j] <- out[[i]]$pc[, j] * -1
          }
        } else {
          newPos <- which.max(abs(out[[i]]$alpha[[j]]))
          signNeg <- Re(out[[i]]$alpha[[j]][newPos]) < 0
          
          if (signNeg != signReference$allNegative[[j]][newPos]) {
            out[[i]]$alpha[[j]] <- out[[i]]$alpha[[j]] * -1
            out[[i]]$pc[, j] <- out[[i]]$pc[, j] * -1
          }
        }
      }
    }
  }
  
  
  return(out)
}


# Compute the variance of alpha
#
# @param bootEncoding output of computeBootStrapEncoding function
# @param nState number of states
# @param nBasis number of basis functions
#
# @return a list (length number of harmonics) of list (length number of states) of variance matrix
#
# @author Quentin Grimonprez
computeVarianceAlpha <- function(bootEncoding, nState, nBasis) {
  nHarm <- nState * nBasis
  
  varAlpha <- list()
  for (harm in seq_len(nHarm)) {
    varAlpha[[harm]] <- list()
    for (iState in seq_len(nState)) {
      tryCatch(
        {
          varAlpha[[harm]][[iState]] <- var(do.call(
            rbind,
            lapply(bootEncoding, function(x) {
              if (!is.null(x$alpha[[harm]][, iState])) {
                return(Re(x$alpha[[harm]][, iState]))
              }
            })
          ),
          use = "pairwise.complete.obs"
          )
        },
        error = function(e) e
      )
    }
  }
  
  return(varAlpha)
}


# compute the variance of encoding
#
# a_x(t) = sum_i alpha_ix * phi_i(t)
# Var(a_x(t)) = sum_i var(alpha_ix) * phi_i^2(t) + sum_{i<j}  2 * phi_i(t) * phi_j(t) * cov(alpha_ix, alpha_jx)
#
# @author Quentin Grimonprez
computeVarianceEncoding <- function(varAlpha, basisobj, harm = 1, nx = 200) {
  nBasis <- basisobj$nbasis
  phi <- fd(diag(nBasis), basisobj)
  nState <- length(varAlpha[[harm]])
  
  timeVal <- seq(basisobj$rangeval[1], basisobj$rangeval[2], length = nx)
  
  Phi <- matrix(nrow = length(timeVal), ncol = nBasis)
  for (i in seq_len(nBasis)) {
    Phi[, i] <- eval.fd(timeVal, phi[i])
  }
  
  funcVar <- list()
  for (iState in seq_len(nState)) {
    funcVar[[iState]] <- rep(NA, nx)
    for (j in seq_along(timeVal)) {
      varAlpha[[harm]][[iState]][is.na(varAlpha[[harm]][[iState]])] <- 0
      funcVar[[iState]][j] <- Phi[j, , drop = FALSE] %*% varAlpha[[harm]][[iState]] %*% t(Phi[j, , drop = FALSE])
    }
  }
  
  return(funcVar)
}



