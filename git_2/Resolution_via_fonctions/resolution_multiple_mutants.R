rm(list=ls())
setwd("/home/gollivier/Documents/Modeles_hotes/git/Modele_hote/Resolution_via_fonctions/sys_equations")
source("../creation_systeme_OK.R")
setwd("/home/gollivier/Documents/Modeles_hotes/git/Modele_hote/Resolution_via_fonctions")

##### INITIALISATION #####
k=0.1 #Valeur d'un k initial
X=1   #Densité de susceptibles initiale
Y=0.2 #Densité d'infectés initiale
variationk=0.01
nbmutant=20 #Nb totale de mutations effectuées (-1)
pourcentpopmutantX=0.1
pourcentpopmutantY=0.1
Xfinal=matrix(0 , nrow = 1/variationk, ncol = nbmutant)
Yfinal=matrix(0 , nrow = 1/variationk, ncol = nbmutant)
############################

for(mut in 1:nbmutant){
          #On définit les paramètres des souches en présence#
  source("fonction_parametres.R")
  param_k=func_param(k)
  
          #On résoud le système#
  source("resolution_systemes_nequ.R")
  resolution=resolution_systemes_nequ(param_k,X,Y)
  
  #Là, on a des 0 si proche de 0, donc on peut savoir si on a disparition et stocker tout ça dans les matrices Xfinal et Yfinal
  for(i in 1:(length(resolution)/2)){
    Xfinal[k[i]*(1/variationk),mut]=resolution[-1+2*i]
    Yfinal[k[i]*(1/variationk),mut]=resolution[2*i]
  }
  
  #On supprime les k disparus 
  k_inter=NULL
  X_inter=NULL
  Y_inter=NULL
  for(j in 1:(length(resolution)/2)){
    if(resolution[2*j]+resolution[2*j-1]>0){
      k_inter=c(k_inter,k[j])
      X_inter=c(X_inter,resolution[2*j-1])
      Y_inter=c(Y_inter,resolution[2*j])
    }
  }
  #Puis on remplace les vecteurs
  k=k_inter ; X=X_inter ; Y=Y_inter
  
  #On va ajouter le mutant
  vecmutant=(X+Y)/sum(X+Y)#vecteur de proba de muter en fonction de la densité relative de la population
  signe=sample(c(-1,1),1)
  kmutant=sample(k,1,prob=vecmutant)+signe*(variationk)
  Xmutant=sum(X)*pourcentpopmutantX
  Ymutant=sum(Y)*pourcentpopmutantY
  #print(kmutant)
  #Et on fusionne ça
  k=c(k,kmutant) ; X=c(X,Xmutant) ; Y=c(Y,Ymutant)
  #print(mut)
}
