#creation d'un fichier qui va écrire un fichier de fonction finale à exécuter hote_h en fonction du nombre de paramètres

setwd("~/Documents/Modeles_hotes/sys_equations")
rm(list=ls())

# creation du jeu d'equation
equa_diff <- function(i,n){
  return(paste(
    paste("dX_",i,"<- ak*X[",i,"] - qk*H*X[",i,"] - bk*X[",i,"] - sum(betak*X[",i,"]*sum(Y[1:",n,"])) + gammak*Y[",i,"],", sep=""),
    paste("dY_",i,"<- sum(betak*X[",i,"]*Y[1:",n,"]) - (alphak+bk+gammak)*Y[",i,"],", sep=""))
  )
}

#Creation du fichier
n<-15
for(c in 1:n){
  write(paste("hote_h <- function(t,pop , param_k){
  with(as.list(param_k), {return(list(c("), file = paste("data_", c, sep=""), append=T)
  for(j in 1:c){
    write(equa_diff(j,c), file = paste("data_", c, sep=""), append=T)
  }
  write("dH=", file = paste("data_", c, sep=""), append=T)
  for(comp in 1:c){
  write(paste("ak*X[",comp,"] - qk*H*X[",comp,"] - bk*X[",comp,"] - (alphak+bk)*Y[",c,"]+",sep=""), file = paste("data_", c, sep=""), append=T)
  }
  write("0)))
})
  }", file = paste("data_", c, sep=""), append=T)
}


#hote_h <- function(t,pop , param_k){
 # with(as.list(param_k), {
  #  dX
  #  dY
  #  dH <- dX_1 + dY_1 + dX_2 + dY_2
  #  return(list(c(dX_1, dY_1, dX_2, dY_2, dH)))
  #})
#}

hote_h <- function(t,pop , param_k){
  with(as.list(param_k), {return(list(c(
    dX_1<- ak*X[1] - qk*H*X[1] - bk*X[1] - sum(betak*X[1]*sum(Y[1:1])) + gammak*Y[1], dY_1<- sum(betak*X[1]*Y[1:1]) - (alphak+bk+gammak)*Y[1],
    dH=
      ak*X[1] - qk*H*X[1] - bk*X[1] - (alphak+bk)*Y[1]+
      0)))
  })
}