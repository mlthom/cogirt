# This goes with the command line argument
# Example 1
# Rscript test.R 1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1 1 20 "/Users/mlthom/Dropbox/research/thomas/CogIRT/data/sdirtSS.rda"
# args equvalent is
# myArgs=c("1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1", "1", "20, "/Users/mlthom/Dropbox/research/thomas/CogIRT/data/sdirtSS.rda")

# Example 2
# Rscript test.R 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 5,1 20 "/Users/mlthom/Dropbox/research/thomas/CogIRT/data/sdirtSS.rda"
# args equvalent is
# myArgs=c("0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0", "5,1", 20, "/Users/mlthom/Dropbox/research/thomas/CogIRT/data/sdirtSS.rda")

# Example 3
# Rscript test.R 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 5,1 20 "/Users/mlthom/Dropbox/research/thomas/CogIRT/data/sdirtSS.rda"
# args equvalent is
# myArgs=c("1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1", "5,1", 20, "/Users/mlthom/Dropbox/research/thomas/CogIRT/data/sdirtSS.rda")

# Example 3
# Rscript test.R 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 5,1 20 "/Users/mlthom/Dropbox/research/thomas/CogIRT/data/sdirtSS.rda"
# args equvalent is
# myArgs=c("1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0", "5,1", 20, "/Users/mlthom/Dropbox/research/thomas/CogIRT/data/sdirtSS.rda")


# Fetch command line arguments'
myArgs <- commandArgs(trailingOnly = TRUE)

#once package is installed, can use first option
if(1 == 0){
  library(CogIRTRUE)
  #ultimatley, may not want to load entire package. we'll see
} else {
  library(devtools)
  load_all(path = "/Users/mlthom/Dropbox/research/thomas/CogIRT/")
}

#here, I'm loading a file (sdirtSS) that is part of package being loaded above.
#in the future though. this will exist elsewhere as part of the psychopy program

load(file = myArgs[4])

#note. this line is not needed if the python file is always called rda
rda <- sdirtSS

#in the future, I the list index should be part of the rda file. still not sure how
#to best do that
rda$list <- c(sapply(X = 1:(length(rda$y) / 20), FUN = rep, 20))

#this would already be true, or some would be NA??
#this par is why it would be nice to use the rpy2 package

rda$y[c(sapply(
  match(as.numeric(unlist(strsplit(myArgs[2], split=","))),rda$list),
  function(x){
    seq(x, x + as.numeric(myArgs[3]) - 1, 1)
  }
))] <- as.numeric(unlist(strsplit(myArgs[1], split=",")))

rda$y[which(!rda$list %in% as.numeric(unlist(strsplit(myArgs[2], split=","))))] <- NA


# print(myArgs[1])
# print("joe")
# print(myArgs[2])
# print("bib")
# print(myArgs[3])
#
#print(rda$y)

rda$omega_sigma2 <- 2 * rda$omega_sigma2

res <- cog_cat(rda = rda, obj_fun = dich_response_model, int_par = 1)
print(res$next_list)
