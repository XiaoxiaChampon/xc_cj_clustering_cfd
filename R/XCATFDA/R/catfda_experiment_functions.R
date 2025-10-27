
source("R/XCATFDA/R/catfda_main.R")

RunExperiment <- function(scenario, num_replicas, est_choice, some_identifier="noid")
{
  temp_folder <- file.path("outputs", "clustersims", paste(scenario, "_", num_replicas, "_", est_choice, "_", some_identifier, sep=""))
  # Empty the directory if it exists
  if(dir.exists(temp_folder)){
    unlink(temp_folder, recursive = TRUE)
  }
  dir.create(temp_folder)
  print(temp_folder)
  
  # scenario="C"
  # num_replicas=3
  # est_choice="probit"
  n100t300C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  # n100t750C <- ClusterSimulation(100,750,scenario,num_replicas,est_choice,TRUE,temp_folder)
  # n100t2000C <- ClusterSimulation(100,2000,scenario,num_replicas,est_choice,TRUE,temp_folder)
  # 
  # 
  # n500t300C <- ClusterSimulation(500,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  # n500t750C <- ClusterSimulation(500,750,scenario,num_replicas,est_choice,TRUE,temp_folder)
  # n500t2000C <- ClusterSimulation(500,2000,scenario,num_replicas,est_choice,TRUE,temp_folder)
  # 
  # 
  # n1000t300C <- ClusterSimulation(1000,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  # n1000t750C <- ClusterSimulation(1000,750,scenario,num_replicas,est_choice,TRUE,temp_folder)
  # n1000t2000C <- ClusterSimulation(1000,2000,scenario,num_replicas,est_choice,TRUE,temp_folder)
  
  n100t750C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n100t2000C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  
  
  n500t300C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n500t750C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n500t2000C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  
  
  n1000t300C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n1000t750C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  n1000t2000C <- ClusterSimulation(100,300,scenario,num_replicas,est_choice,TRUE,temp_folder)
  
  true_tableC <- rbind(n100t300C$cluster_table_true,n100t750C$cluster_table_true,n100t2000C$cluster_table_true,
                       n500t300C$cluster_table_true,n500t750C$cluster_table_true,n500t2000C$cluster_table_true,
                       n1000t300C$cluster_table_true,n1000t750C$cluster_table_true,n1000t2000C$cluster_table_true)
  rownames(true_tableC) <- c("n100t300","n100t750","n100t2000",
                             "n500t300","n500t750","n500t2000",
                             "n1000t300","n1000t750","n1000t2000")
  
  est_tableC <- rbind(n100t300C$cluster_table_est,n100t750C$cluster_table_est,n100t2000C$cluster_table_est,
                      n500t300C$cluster_table_est,n500t750C$cluster_table_est,n500t2000C$cluster_table_est,
                      n1000t300C$cluster_table_est,n1000t750C$cluster_table_est,n1000t2000C$cluster_table_est)
  rownames(est_tableC) <- c("n100t300","n100t750","n100t2000",
                            "n500t300","n500t750","n500t2000",
                            "n1000t300","n1000t750","n1000t2000")
  
  est_tableC_se <- rbind(n100t300C$cluster_table_est_se,n100t750C$cluster_table_est_se,n100t2000C$cluster_table_est_se,
                         n500t300C$cluster_table_est_se,n500t750C$cluster_table_est_se,n500t2000C$cluster_table_est_se,
                         n1000t300C$cluster_table_est_se,n1000t750C$cluster_table_est_se,n1000t2000C$cluster_table_est_se)
  rownames(est_tableC_se) <- c("n100t300","n100t750","n100t2000",
                               "n500t300","n500t750","n500t2000",
                               "n1000t300","n1000t750","n1000t2000")
  
  mse_tableC100 <- rbind(
    c(n100t300C$mse[1,1],n100t300C$mse[2,1],n100t300C$hellinger[1,1],n100t300C$hellinger[2,1],n100t300C$hellinger[3,1]),
    c(n100t750C$mse[1,1],n100t750C$mse[2,1],n100t750C$hellinger[1,1],n100t750C$hellinger[2,1],n100t750C$hellinger[3,1]),
    c(n100t2000C$mse[1,1],n100t2000C$mse[2,1],n100t2000C$hellinger[1,1],n100t2000C$hellinger[2,1],n100t2000C$hellinger[3,1]),
    
    c(n100t300C$mse[1,2],n100t300C$mse[2,2],n100t300C$hellinger[1,2],n100t300C$hellinger[2,2],n100t300C$hellinger[3,2]),
    c(n100t750C$mse[1,2],n100t750C$mse[2,2],n100t750C$hellinger[1,2],n100t750C$hellinger[2,2],n100t750C$hellinger[3,2]),
    c(n100t2000C$mse[1,2],n100t2000C$mse[2,2],n100t2000C$hellinger[1,2],n100t2000C$hellinger[2,2],n100t2000C$hellinger[3,2]),
    
    c(n100t300C$mse[1,3],n100t300C$mse[2,3],n100t300C$hellinger[1,3],n100t300C$hellinger[2,3],n100t300C$hellinger[3,3]),
    c(n100t750C$mse[1,3],n100t750C$mse[2,3],n100t750C$hellinger[1,3],n100t750C$hellinger[2,3],n100t750C$hellinger[3,3]),
    c(n100t2000C$mse[1,3],n100t2000C$mse[2,3],n100t2000C$hellinger[1,3],n100t2000C$hellinger[2,3],n100t2000C$hellinger[3,3])
  )
  
  rownames(mse_tableC100) <- c("s1n100t300","s1n100t750","s1n100t2000",
                               "s2n100t300","s2n100t750","s2n100t2000",
                               "s3n100t300","s3n100t750","s3n100t2000")
  colnames(mse_tableC100) <- c("z1","z2","p1","p2","p3")
  
  mse_tableC500 <- rbind(
    c(n500t300C$mse[1,1],n500t300C$mse[2,1],n500t300C$hellinger[1,1],n500t300C$hellinger[2,1],n500t300C$hellinger[3,1]),
    c(n500t750C$mse[1,1],n500t750C$mse[2,1],n500t750C$hellinger[1,1],n500t750C$hellinger[2,1],n500t750C$hellinger[3,1]),
    c(n500t2000C$mse[1,1],n500t2000C$mse[2,1],n500t2000C$hellinger[1,1],n500t2000C$hellinger[2,1],n500t2000C$hellinger[3,1]),
    
    c(n500t300C$mse[1,2],n500t300C$mse[2,2],n500t300C$hellinger[1,2],n500t300C$hellinger[2,2],n500t300C$hellinger[3,2]),
    c(n500t750C$mse[1,2],n500t750C$mse[2,2],n500t750C$hellinger[1,2],n500t750C$hellinger[2,2],n500t750C$hellinger[3,2]),
    c(n500t2000C$mse[1,2],n500t2000C$mse[2,2],n500t2000C$hellinger[1,2],n500t2000C$hellinger[2,2],n500t2000C$hellinger[3,2]),
    
    c(n500t300C$mse[1,3],n500t300C$mse[2,3],n500t300C$hellinger[1,3],n500t300C$hellinger[2,3],n500t300C$hellinger[3,3]),
    c(n500t750C$mse[1,3],n500t750C$mse[2,3],n500t750C$hellinger[1,3],n500t750C$hellinger[2,3],n500t750C$hellinger[3,3]),
    c(n500t2000C$mse[1,3],n500t2000C$mse[2,3],n500t2000C$hellinger[1,3],n500t2000C$hellinger[2,3],n500t2000C$hellinger[3,3])
  )
  
  
  rownames(mse_tableC500) <- c("s1n500t300","s1n500t750","s1n500t2000",
                               "s2n500t300","s2n500t750","s2n500t2000",
                               "s3n500t300","s3n500t750","s3n500t2000")
  colnames(mse_tableC500) <- c("z1","z2","p1","p2","p3")
  
  
  mse_tableC1000 <- rbind(
    c(n1000t300C$mse[1,1],n1000t300C$mse[2,1],n1000t300C$hellinger[1,1],n1000t300C$hellinger[2,1],n1000t300C$hellinger[3,1]),
    c(n1000t750C$mse[1,1],n1000t750C$mse[2,1],n1000t750C$hellinger[1,1],n1000t750C$hellinger[2,1],n1000t750C$hellinger[3,1]),
    c(n1000t2000C$mse[1,1],n1000t2000C$mse[2,1],n1000t2000C$hellinger[1,1],n1000t2000C$hellinger[2,1],n1000t2000C$hellinger[3,1]),
    
    c(n1000t300C$mse[1,2],n1000t300C$mse[2,2],n1000t300C$hellinger[1,2],n1000t300C$hellinger[2,2],n1000t300C$hellinger[3,2]),
    c(n1000t750C$mse[1,2],n1000t750C$mse[2,2],n1000t750C$hellinger[1,2],n1000t750C$hellinger[2,2],n1000t750C$hellinger[3,2]),
    c(n1000t2000C$mse[1,2],n1000t2000C$mse[2,2],n1000t2000C$hellinger[1,2],n1000t2000C$hellinger[2,2],n1000t2000C$hellinger[3,2]),
    
    c(n1000t300C$mse[1,3],n1000t300C$mse[2,3],n1000t300C$hellinger[1,3],n1000t300C$hellinger[2,3],n1000t300C$hellinger[3,3]),
    c(n1000t750C$mse[1,3],n1000t750C$mse[2,3],n1000t750C$hellinger[1,3],n1000t750C$hellinger[2,3],n1000t750C$hellinger[3,3]),
    c(n1000t2000C$mse[1,3],n1000t2000C$mse[2,3],n1000t2000C$hellinger[1,3],n1000t2000C$hellinger[2,3],n1000t2000C$hellinger[3,3])
  )
  
  
  rownames(mse_tableC1000) <- c("s1n1000t300","s1n1000t750","s1n1000t2000",
                                "s2n1000t300","s2n1000t750","s2n1000t2000",
                                "s3n1000t300","s3n1000t750","s3n1000t2000")
  colnames(mse_tableC1000) <- c("z1","z2","p1","p2","p3")
  
  
  
  #save(true_tableC,est_tableC,est_tableC_se,mse_tableC100,mse_tableC500,mse_tableC1000,
  #file=paste(scenario,num_replicas,est_choice,"clustering.RData",sep="_"))
  
  true_est_w_data <- (list("true_tableC"=true_tableC,
                           "est_tableC"=est_tableC,
                           "est_tableC_se"=est_tableC_se,
                           "mse_tableC100"=mse_tableC100,
                           "mse_tableC500"=mse_tableC500,
                           "mse_tableC1000"=mse_tableC1000,
                           
                           "Z_true_curves"=list(n100t300C$Z_true_curves,
                                                n100t750C$Z_true_curves,
                                                n100t2000C$Z_true_curves,
                                                n500t300C$Z_true_curves,
                                                n500t750C$Z_true_curves,
                                                n500t2000C$Z_true_curves,
                                                n1000t300C$Z_true_curves,
                                                n1000t750C$Z_true_curves,
                                                n1000t2000C$Z_true_curves),
                           "Z_est_curves"=list(n100t300C$Z_est_curves,
                                               n100t750C$Z_est_curves,
                                               n100t2000C$Z_est_curves,
                                               n500t300C$Z_est_curves,
                                               n500t750C$Z_est_curves,
                                               n500t2000C$Z_est_curves,
                                               n1000t300C$Z_est_curves,
                                               n1000t750C$Z_est_curves,
                                               n1000t2000C$Z_est_curves),
                           "p_true_curves"=list(n100t300C$p_true_curves,
                                                n100t750C$p_true_curves,
                                                n100t2000C$p_true_curves,
                                                n500t300C$p_true_curves,
                                                n500t750C$p_true_curves,
                                                n500t2000C$p_true_curves,
                                                n1000t300C$p_true_curves,
                                                n1000t750C$p_true_curves,
                                                n1000t2000C$p_true_curves),
                           "p_est_curves"=list(n100t300C$p_est_curves,
                                               n100t750C$p_est_curves,
                                               n100t2000C$p_est_curves,
                                               n500t300C$p_est_curves,
                                               n500t750C$p_est_curves,
                                               n500t2000C$p_est_curves,
                                               n1000t300C$p_est_curves,
                                               n1000t750C$p_est_curves,
                                               n1000t2000C$p_est_curves),
                           "W_cfd"=list(n100t300C$W_cfd,
                                        n100t750C$W_cfd,
                                        n100t2000C$W_cfd,
                                        n500t300C$W_cfd,
                                        n500t750C$W_cfd,
                                        n500t2000C$W_cfd,
                                        n1000t300C$W_cfd,
                                        n1000t750C$W_cfd,
                                        n1000t2000C$W_cfd)
  ))
  
  save(true_est_w_data,file=file.path("outputs", paste(scenario,num_replicas,est_choice,some_identifier,"true_est_w_data_clustering.RData",sep="_")))
  
  return(true_est_w_data)
}