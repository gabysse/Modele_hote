#creation d'un fichier qui va écrire un fichier de fonction finale à exécuter hote_h en fonction du nombre de paramètres

setwd("~/Documents/Modeles_hotes/sys_equations")
rm(list=ls())

# Creation du jeu d'equation
# Attention, il faut remplacer les ak et autres par leurs expressions pour i 

equa_diff <- function(i,imax){
  towrite <- paste("dX_",i,"<- ak[",i,"]*pop[",i*2-1,"] - qk[",i,"]*sum(pop)*pop[",i*2-1,"] - bk[",i,"]*pop[",i*2-1,"] - betak[",i,"]*pop[",i*2-1,"]*sum(pop[2*1:",imax,"])+ gammak[",i,"]*pop[",i*2,"], dY_",i,"<- betak[",i,"]*pop[",i*2-1,"]*sum(pop[2*1:",imax,"])- (alphak[",i,"]+bk[",i,"]+gammak[",i,"])*pop[",i*2,"],", sep="")
  return(towrite)
}

equa_diff_fin <- function(i,imax){
  towrite <- paste("dX_",i,"<- ak[",i,"]*pop[",i*2-1,"] - qk[",i,"]*sum(pop)*pop[",i*2-1,"] - bk[",i,"]*pop[",i*2-1,"] - betak[",i,"]*pop[",i*2-1,"]*sum(pop[2*1:",imax,"])+ gammak[",i,"]*pop[",i*2,"], dY_",i,"<- betak[",i,"]*pop[",i*2-1,"]*sum(pop[2*1:",imax,"])- (alphak[",i,"]+bk[",i,"]+gammak[",i,"])*pop[",i*2,"]", sep="")
  return(towrite)
}

#Creation du fichier
n<-15
for(c in 1:n){
  # test d'existence du fichier ; s'il existe, on le supprime. 
  if(file.exists(paste("data_", c, sep=""))){file.remove(paste("data_", c, sep=""))}
  write(paste("hote_h <- function(t,pop , param_k){
  
  with(as.list(param_k), {return(list(c("), file = paste("data_", c, sep=""), append=T)
  if(c==1){
    write(equa_diff_fin(c,c),file="data_1",append=T)
  }else{
  for(j in 1:(c-1)){
    write(equa_diff(j,c), file = paste("data_", c, sep=""), append=T)
  }
    write(equa_diff_fin(c,c), file = paste("data_", c, sep=""), append=T)
  }
  write(")))})}", file = paste("data_", c, sep=""), append=T)
}