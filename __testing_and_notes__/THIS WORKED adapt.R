rm(list=ls())


## Load Package
# later on, when package is installed, use this
# library(CogIRTRUE)
# but for now, use this
library(devtools)
load_all(path = "/Users/mlthom/Dropbox/research/thomas/CogIRT/")

if(1==1){
  #update this later
  my_args <- "/Users/mlthom/Dropbox/research/thomas/rdoc/psychopy/nbackCAT/data/_nbackCAT_2021_Jan_22_1632.csv"
  #load rda file
  load("/Users/mlthom/Dropbox/research/thomas/rdoc/psychopy/nbackCAT/stimulus_generation/nbackcat.rda")
  #load response data
  dat <- read.csv(file = my_args[1])
  rda <- nbackcat
  #rda <- sdirtSS
  #rda$y <- array(NA, c(1, ncol(rda$y)))
  #rda$y[-1 *c(1:20)] <- NA
  rda$y[1:20] <- dat$test_resp.corr
  #rda$lambda <- nbackcat$lambda
  #rda$nu_mu <- nbackcat$nu_mu
  #rda$nu_sigma2 <- nbackcat$nu_sigma2
  #rda$omega_sigma2 <- nbackcat$omega_sigma2
  #rda$zeta_mu <- nbackcat$zeta_mu
  #rda$zeta_sigma2 <- nbackcat$zeta_sigma2
  #rda <- rda[-1*which(names(rda) %in% c("ystar", "nu", "omega", "zeta", "condition", "item_type", "list"))]
  #rda$list <- nbackcat$list
  #rda <- rda[-c(2:3,6:7,14:15)]
  #rda$list <- c(sapply(X = 1:(length(rda$y) / 20), FUN = rep, 20))


  #update rda file with response data
  #rda$y[,which(rda$list== 1)] <- dat$test_resp.corr
  #rda$y[,which(rda$list %in% gsub(pattern = "stim_files/lists/", "", dat$list))] <- dat$test_resp.corr
  #rda$y[which(!rda$list %in% c(1))] <- NA
  #rda$y <- matrix(as.double(rda$y), c(1, 1000))
  #rda$list <- as.integer(factor(rda$list))
  #determine update
  res <- cog_cat(rda = rda, obj_fun = dich_response_model, int_par = 1)
  cat(res$next_list)
  #cat(1)
  # Fetch command line arguments'
  #myArgs <- commandArgs(trailingOnly = TRUE)

  # Convert to numerics
  #nums = as.numeric(myArgs)

  # cat will write the result to the stdout stream
  #cat(max(nums))

}

if(1==0){
  rda <- sdirtSS
  #
  # #***
  rda <- rda[-c(2:3,6:7,14:15)]
  rda$list <- c(sapply(X = 1:(length(rda$y) / 20), FUN = rep, 20))
  #View(rda)
  #
  rda$y[which(!rda$list %in% c(1))] <- NA
  res <- cog_cat(rda = rda, obj_fun = dich_response_model, int_par = 1)
  cat(res$next_list)
}
