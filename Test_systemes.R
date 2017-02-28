rm(list=ls())
library(deSolve)
setwd("~/Documents/Modeles_hotes") 
k=c(0.1,0.11)

source("./parametres/valeurs_parametres_puissance.R") #Contient les variables definissant les TO.
source("./parametres/expression_parametres_puissance.R") #Contient les fonctions exprimant les TO.

#beta_k#
betak=beta(beta0,beta_max,P_beta,k)

# a_k : Definition #
ak=a(a0,a_max,P_a,k)

# alpha_k #
alphak=alpha(alpha0,alpha_max,P_alpha,k)

# gamma_k #
gammak=gamma(gamma0,gamma_max,P_gamma,k)

# b_k #
bk=b(b0,b_max,P_b,k)

# q_k #
qk=cu(ak,bk,K)

# LIST !!!!!
param_k=list(betak=betak,ak=ak,alphak=alphak,gammak=gammak,bk=bk,qk=qk)


X0=c(3.5,0.35)
Y0=c(2.5,0.25)
V0 <- NULL
for (i in 1:length(X0)) V0 <- c(V0, X0[i], Y0[i])



source(paste("sys_equations/data_", length(k), sep=""), local=T)
final=ode(y=V0, t=seq(0,500,0.1), func=hote_h, parms=param_k)
final2=ode(y=V0, t=c(0,2,5,10,seq(20,500,0.1)), func=hote_h, parms=param_k)

par(mfrow=c(3,2))
plot(final[,1],final[,2],main="Sains (k=0.1)",xlab="Temps", ylab="Densité",t="l")
lines(final2[,1],final2[,2],t="p")
plot(final[,1],final[,3],main="Infectes (k=0.1)",xlab="Temps", ylab="Densité",t="l")
lines(final2[,1],final2[,3],t="p")
plot(final[,1],final[,4],main="Sains (k=0.2)",xlab="Temps", ylab="Densité",t="l")
lines(final2[,1],final2[,4],t="p")
plot(final[,1],final[,5],main="Infectes (k=0.2)",xlab="Temps", ylab="Densité",t="l")
lines(final2[,1],final2[,5],t="p")
plot(final[,1],final[,6],main="Sains (k=0.25)",xlab="Temps", ylab="Densité",t="l")
lines(final2[,1],final2[,6],t="p")
plot(final[,1],final[,7],main="Infectes (k=0.25)",xlab="Temps", ylab="Densité",t="l")
lines(final2[,1],final2[,7],t="p")