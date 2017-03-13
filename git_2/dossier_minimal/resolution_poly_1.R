rm(list=ls())
library(deSolve)
#set.seed(123)
pdf("poly_1.pdf", width=7, height=7)
par(mfrow=c(2,2))
print(Sys.time())

source("fonction_parametres.R")
source("resolution_systemes_nequ_v2.R")
###### INITIALISATION ######
############################
setwd("sys_equations")
source("../creation_systeme_OK.R")
setwd("../")
source("Reso_1302_v2.R")
par(mfrow=c(2,2))
total=4
for(i in 1:total){

  beta0=3
  beta_max=5
  P_beta=1
  
  a0=3.577424
  a_max=4.392827
  P_a=i/(total-1)
  
  alpha0=1
  alpha_max=1
  P_alpha=0
  
  gamma0=1
  gamma_max=1
  P_gamma=0
  
  b0=1
  b_max=1
  P_b=0
  
  q=1

  final=resolution_multiple(nbdek=40,nbmutant=250,Xini=0.1,Yini=0.1,timebfmut=10^10,dec_arrondi=2,Ymutant=10,Xmutant=25,pos_k_ini=3)

  image(final,col=gray(rev(seq(from=0,to=1,by=0.01))),main=paste("P_a=",i/(total-1)))

}
dev.off()
