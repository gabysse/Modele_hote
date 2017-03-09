rm(list=ls())
library(deSolve)

###### INITIALISATION ######
k=c(0.4) #Valeur d'un k initial
kinitial=k #Va servir uniquement en cas d'erreur
X=c(1)   #Densite de susceptibles initiale
Y=c(1) #Densite d'infectes initiale
variationk=0.05
nbmutant=5 #Nb totale de mutations effectuees (+1)
pourcentpopmutantX=0.9
pourcentpopmutantY=0.9
Xfinal=matrix(0 , nrow = 1/variationk+1, ncol = nbmutant)
Yfinal=matrix(0 , nrow = 1/variationk+1, ncol = nbmutant)
kfinal=rep(0,1/variationk)

############################
# setwd("sys_equations")
# source("../creation_systeme_OK.R")
# setwd("../")

for(mut in 1:nbmutant){
  nbsouche=length(k)
  if(nbsouche>1){ #On regarde si on a bien la possibilité de faire les tests suivants, il faut pour celà avoir au moins 2 souches
    if(any(as.integer(round(k[-nbsouche]*1/variationk))==as.integer(round(k[nbsouche]*1/variationk)))){#On regarde si le dernier k ajouté (le k mutant ) est égal à un k déja existant
      #Si c'est le cas on va alors enlever cette valeur
      k=k[-nbsouche] 
      X=X[-nbsouche]
      Y=Y[-nbsouche]
      for(i in 1:nbsouche){#Et on associe les valeurs de X et Y la résolution précédente aux valeurs présentes 
        Xfinal[k[i]*(1/variationk),mut]=Xfinal[k[i]*(1/variationk),mut]+resolution[-1+2*i]
        Yfinal[k[i]*(1/variationk),mut]=Yfinal[k[i]*(1/variationk),mut]+resolution[2*i]
      }
    }else{#Si le dernier k n'est pas déjà présent, on résoud normalement.
              #On definit les parametres des souches en presence#
      source("fonction_parametres.R")
      param_k=func_param(k)
      
              #On resoud le systeme#
      source("resolution_systemes_nequ_v2.R")
      resolution=resolution_systemes_nequ(param_k,X,Y,nbsouche)
      
      #La, on a des 0 si proche de 0, donc on peut savoir si on a disparition et stocker tout ca dans les matrices Xfinal et Yfinal
      for(i in 1:nbsouche){
        Xfinal[k[i]*(1/variationk),mut]=Xfinal[k[i]*(1/variationk),mut]+resolution[-1+2*i]
        Yfinal[k[i]*(1/variationk),mut]=Yfinal[k[i]*(1/variationk),mut]+resolution[2*i]
      }
      
      #On supprime les k disparus 
      k_inter=NULL
      X_inter=NULL
      Y_inter=NULL
      for(j in 1:nbsouche){ #On recupere les couples resolution[1,2], puis resolution[3,4], etc..., qui representent les valeurs de X et Y pour une souche en particulier.
        if(resolution[2*j]>0 && resolution[2*j-1]>0){
          k_inter=c(k_inter,k[j])
          X_inter=c(X_inter,Xfinal[k[j]*(1/variationk),mut])
          Y_inter=c(Y_inter,Yfinal[k[j]*(1/variationk),mut])
        }
      }
      if(is.null(Y_inter)){#On a ici extinction du parasite et donc equilibre non-trivial 
        write(paste("Disparition du parasite pour k=",kinitial,"avec les parametres beta=",param_k$betak,"alpha=",param_k$alphak,"gamma=",param_k$gammak,"a=",param_k$ak,"b=",param_k$bk,"et q=",param_k$qk,sep=" "),file=paste("Rapports_erreur/Error_for_k=",k,".txt",sep="")) #On écrit un message d'erreur...
        break #...et on sort de la boucle.
      }
      #Puis on remplace les vecteurs
      k=k_inter ; X=X_inter ; Y=Y_inter

    }
    
    #On a traité le cas où on fait une résolution à 2 souches ou plus. Il faut également traiter le cas où l'on a un seule souche (la première résolution)
    
  }else{#On definit les parametres des souches en presence#
    source("fonction_parametres.R")
    param_k=func_param(k)
    
    #On resoud le systeme#
    source("resolution_systemes_nequ.R")
    resolution=resolution_systemes_nequ(param_k,X,Y)
    
    #La, on a des 0 si proche de 0, donc on peut savoir si on a disparition et stocker tout ca dans les matrices Xfinal et Yfinal
    for(i in 1:(length(resolution)/2)){
      Xfinal[k[i]*(1/variationk),mut]=Xfinal[k[i]*(1/variationk),mut]+resolution[-1+2*i]
      Yfinal[k[i]*(1/variationk),mut]=Yfinal[k[i]*(1/variationk),mut]+resolution[2*i]
    }
    
    #On supprime les k disparus 
    k_inter=NULL
    X_inter=NULL
    Y_inter=NULL
    for(j in 1:(length(resolution)/2)){ #On recupere les couples resolution[1,2], puis resolution[3,4], etc..., qui representent les valeurs de X et Y pour une souche en particulier.
      if(resolution[2*j]>0 && resolution[2*j-1]>0){
        k_inter=c(k_inter,k[j])
        X_inter=c(X_inter,Xfinal[k[j]*(1/variationk),mut])
        Y_inter=c(Y_inter,Yfinal[k[j]*(1/variationk),mut])
      }
    }
    if(is.null(Y_inter)){#On a ici extinction du parasite et donc equilibre non-trivial 
      write(paste("Disparition du parasite pour k=",kinitial,"avec les parametres beta=",param_k$betak,"alpha=",param_k$alphak,"gamma=",param_k$gammak,"a=",param_k$ak,"b=",param_k$bk,"et q=",param_k$qk,sep=" "),file=paste("Rapports_erreur/Error_for_k=",k,".txt",sep="")) #On écrit un message d'erreur...
      break #...et on sort de la boucle.
    }
    #Puis on remplace les vecteurs
    k=k_inter ; X=X_inter ; Y=Y_inter
    kfinal=rbind(kfinal,k)
    
  }
  kfinal=rbind(kfinal,k) #On met les k résidents dans la matrice de vérification
  #Il faut maintenant passer à l'étape suivant la résolution du système.
  #On va ajouter le mutant
  
  vecmutant=(X+Y)/sum(X+Y)#vecteur de proba de muter en fonction de la densite relative de la population
  signe=sample(c(-1,1),1)
  ktire=sample(k,1,prob=vecmutant)
  if(ktire<2*variationk){signe=1}
  if(ktire>1-variationk){signe=-1} #Pour rester entre 0 et 1
  kmutant=ktire+signe*(variationk)
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
image(matrix(as.numeric(as.logical(Xfinal+Yfinal)),nrow=1/variationk+1,ncol=nbmutant))
image(matrix(as.numeric(as.logical(Xfinal)),nrow=1/variationk+1,ncol=nbmutant))
image(matrix(as.numeric(as.logical(Yfinal)),nrow=1/variationk+1,ncol=nbmutant))
