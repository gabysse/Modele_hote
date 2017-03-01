#creation d'un fichier qui va écrire un fichier de fonction finale à exécuter hote_h en fonction du nombre de paramètres

# Creation du jeu d'equation
# Attention, il faut remplacer les ak et autres par leurs expressions pour i 

### INITIALISATION ###
n<-(1/variationk)+1
######################

equa_diff <- function(i,imax){
  towrite <- paste("dX_",i,"<- param_k$ak[",i,"]*pop[",i*2-1,"] - param_k$qk[",i,"]*sum(pop)*pop[",i*2-1,"] - param_k$bk[",i,"]*pop[",i*2-1,"] - param_k$betak[",i,"]*pop[",i*2-1,"]*sum(pop[2*1:",imax,"])+ param_k$gammak[",i,"]*pop[",i*2,"], dY_",i,"<- param_k$betak[",i,"]*pop[",i*2-1,"]*sum(pop[2*1:",imax,"])- (param_k$alphak[",i,"]+param_k$bk[",i,"]+param_k$gammak[",i,"])*pop[",i*2,"],", sep="")
  return(towrite)
}

equa_diff_fin <- function(i,imax){
  towrite <- paste("dX_",i,"<- param_k$ak[",i,"]*pop[",i*2-1,"] - param_k$qk[",i,"]*sum(pop)*pop[",i*2-1,"] - param_k$bk[",i,"]*pop[",i*2-1,"] - param_k$betak[",i,"]*pop[",i*2-1,"]*sum(pop[2*1:",imax,"])+ param_k$gammak[",i,"]*pop[",i*2,"], dY_",i,"<- param_k$betak[",i,"]*pop[",i*2-1,"]*sum(pop[2*1:",imax,"])- (param_k$alphak[",i,"]+param_k$bk[",i,"]+param_k$gammak[",i,"])*pop[",i*2,"]", sep="")
  return(towrite)
}

#Creation du fichier

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
