rm(list=ls())
library(deSolve)

###### INITIALISATION ######
k=0.95 #Valeur d'un k initial
kinitial=k #Va servir uniquement en cas d'erreur
X=1   #Densite de susceptibles initiale
Y=0.2 #Densite d'infectes initiale
variationk=0.01
nbmutant=25 #Nb totale de mutations effectuees (+1)
pourcentpopmutantX=0.1
pourcentpopmutantY=0.1
Xfinal=matrix(0 , nrow = 1/variationk+1, ncol = nbmutant)
Yfinal=matrix(0 , nrow = 1/variationk+1, ncol = nbmutant)
############################
setwd("sys_equations")
source("../creation_systeme_OK.R")
setwd("../")

for(mut in 1:nbmutant){
  #On definit les parametres des souches en presence#
  source("fonction_parametres.R")
  param_k=func_param(k)
  
  #On resoud le systeme#
  source("resolution_systemes_nequ.R")
  resolution=resolution_systemes_nequ(param_k,X,Y)
  
  #La, on a des 0 si proche de 0, donc on peut savoir si on a disparition et stocker tout ca dans les matrices Xfinal et Yfinal
  for(i in 1:(length(resolution)/2)){
    Xfinal[k[i]*(1/variationk),mut]=resolution[-1+2*i]
    Yfinal[k[i]*(1/variationk),mut]=resolution[2*i]
  }
  
  #On supprime les k disparus 
  k_inter=NULL
  X_inter=NULL
  Y_inter=NULL
  for(j in 1:(length(resolution)/2)){ #On recupere les couples resolution[1,2], puis resolution[3,4], etc..., qui representent les valeurs de X et Y pour une souche en particulier.
    if(resolution[2*j]>0 && resolution[2*j-1]>0){
      k_inter=c(k_inter,k[j])
      X_inter=c(X_inter,resolution[2*j-1])
      Y_inter=c(Y_inter,resolution[2*j])
    }
  }
  if(is.null(X_inter)){#On a ici extinction du parasite et donc equilibre non-trivial 
    write(paste("Disparition du parasite pour k=",kinitial,"avec les parametres beta=",param_k$betak,"alpha=",param_k$alphak,"gamma=",param_k$gammak,"a=",param_k$ak,"b=",param_k$bk,"et q=",param_k$qk,sep=" "),file=paste("Rapports_erreur/Error_for_k=",k,".txt",sep="")) #On écrit un message d'erreur...
    break #...et on sort de la boucle.
  }
  #Puis on remplace les vecteurs
  k=k_inter ; X=X_inter ; Y=Y_inter
  
  #On va ajouter le mutant
  vecmutant=(X+Y)/sum(X+Y)#vecteur de proba de muter en fonction de la densite relative de la population
  signe=sample(c(-1,1),1)
  if(any(k==variationk)){signe=1}
  if(any(k==(1-variationk))){signe=-1}
  kmutant=sample(k,1,prob=vecmutant)+signe*(variationk)
  Xmutant=sum(X)*pourcentpopmutantX
  Ymutant=sum(Y)*pourcentpopmutantY
  #print(kmutant)
  #Et on fusionne ca
  k=c(k,kmutant) ; X=c(X,Xmutant) ; Y=c(Y,Ymutant)
  #print(mut)
}
par(mfrow=c(1,3))
image(Xfinal+Yfinal)
image(Xfinal)
image(Yfinal)