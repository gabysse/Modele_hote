### fonction permettant d'implementer un nos systemes d'equations ###

# necessite d'un vecteur k contenant nos differentes valeurs de k

setwd("~/Documents/Modeles_hotes/sys_equations")
rm(list=ls())

# creation du jeu d'equation
equa_diff <- function(i,n){
  return(paste(
    paste("dX_",i,"<- ak*X[",i,"] - qk*H*X[",i,"] - bk*X[",i,"] - sum(betak[",i,"]*X[",i,"]*sum(Y[1:",n,"])) + gammak*Y[",i,"],", sep=""),
    paste("dY_",i,"<- sum(betak[",i,"]*X[",i,"]*Y[",i,"]) - (alphak+bk+gammak)*Y[",i,"],", sep=""))
  )
}

n<-15
for(c in 1:n){
  write("systeqdif=c(", file = paste("data_", c, sep=""), append=T)
  for(j in 1:c){
   write(equa_diff(j,c), file = paste("data_", c, sep=""), append=T)
  }
  write(paste("dH=sum(ak*X[1:",c,"] - qk*H*X[1:",c,"] - bk*X[1:",c,"] - (alphak+bk)*Y[1:",c,"]))",sep=""), file = paste("data_", c, sep=""), append=T)
}



