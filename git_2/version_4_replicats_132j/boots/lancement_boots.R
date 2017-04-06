val_lanc=c(8,20,35)
for(lanc in val_lanc){
  setwd(paste("kini_",lanc,sep=""))
  source("boots.R")
  setwd("../")
}
