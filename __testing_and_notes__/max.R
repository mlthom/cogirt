
## test cog cat functionality

# later on, when package is installed, use this
#library(CogIRTRUE)

# but for now, use this
library(devtools)
load_all(path = "/Users/mlthom/Dropbox/research/thomas/CogIRT/")

rda <- sdirtSS
rda$list <- c(sapply(X = 1:(length(rda$y) / 5), FUN = rep, nrep = 5))
rda$y[which(!rda$list %in% c(1))] <- NA
res <- cog_cat(rda = rda, obj_fun = dich_response_model, int_par = 1)
cat(res$next_list)

# Fetch command line arguments'
#myArgs <- commandArgs(trailingOnly = TRUE)

# Convert to numerics
#nums = as.numeric(myArgs)

# cat will write the result to the stdout stream
#cat(max(nums))
