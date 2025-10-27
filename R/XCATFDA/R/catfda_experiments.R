

run_parallel <- TRUE
time_elapsed <- list()
if(run_parallel)
{
  print("RUNNING PARALLEL")
  
  # For: makeCluster
  library(doParallel)
  
  # For: %dorng% or registerDoRNG for reproducable parallel random number generation
  library(doRNG)
  
  if(exists("initialized_parallel") && initialized_parallel == TRUE)
  {
    parallel::stopCluster(cl = my.cluster)
  }
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  cat("Parellel Registered: ", foreach::getDoParRegistered(), "\n")
  initialized_parallel <- TRUE
  
  # registerDoRNG(123) # ///<<<< THIS CREATES THE ERROR FOR FADPClust !!!
}

#' Create directories 
if (!dir.exists("outputs")){
  dir.create("outputs")
}
if (!dir.exists("outputs/clustersims")){
  dir.create("outputs/clustersims")
}







# EXECUTION:

# }) # profvis end


#test
# set.seed(123)
#  A_2_probit <- RunExperiment("A",2,"probit")
#  
# 
#  
#  set.seed(123)
#  A_2_mul <- RunExperiment("A",2,"multinormial")
# 

# set.seed(123)
# A_2_probit <- RunExperiment("A",2,"probit","test")
# 
# set.seed(123)
# A_2_binomial <- RunExperiment("A",2,"binomial","test")

set.seed(123)
A_2_multinomial_new <- RunExperiment("A",2,"multinomial","test")

#save(C_2_probit,file="C_2_probit.RData")
# set.seed(123)
# A_100_probit <- RunExperiment("A",100,"probit")
# 
# 
# set.seed(123)
# A_20_multinomial <- RunExperiment("A",20,"multinormial")
# 
# set.seed(123)
# C_20_probit <- RunExperiment("C",20,"probit")
# 
# set.seed(123)
# C_20_multinomial <- RunExperiment("C",20,"multinormial")
# 
# set.seed(123)
# B_20_probit <- RunExperiment("B",20,"probit")
# 
# set.seed(123)
# B_20_multinormial <- RunExperiment("B",20,"multinormial")


if(run_parallel)
{
  parallel::stopCluster(cl = my.cluster)
  initialized_parallel <- FALSE
}
