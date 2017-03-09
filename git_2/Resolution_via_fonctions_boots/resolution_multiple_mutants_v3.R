rm(list=ls())
library(deSolve)
#set.seed(123)
print(Sys.time())
###### INITIALISATION ######
k=c(0.1) #Valeur d'un k initial
kinitial=k #Va servir uniquement en cas d'erreur
X=c(0.5)   #Densite de susceptibles initiale
Y=c(0.25) #Densite d'infectes initiale

nbdek=30
nbmutant=400 #Nb totale de mutations effectuees (+1)
pourcentpopmutantX=0.4
pourcentpopmutantY=0.2
Xfinal=matrix(0 , nrow = nbdek+1, ncol = nbmutant)
Yfinal=matrix(0 , nrow = nbdek+1, ncol = nbmutant)
kfinal=rep(0,nbdek)

############################
# setwd("sys_equations")
# source("../creation_systeme_OK.R")
# setwd("../")

for(mut in 1:nbmutant){
  nbsouche=length(k)
  # print(paste("mut",mut))
  # print(paste("nbsouche=",nbsouche))
  # print(paste("vec de k entree :",k))
  # print(paste("X a lentree :",X))
  # print(paste("Y a lentree :",Y))
  if(nbsouche>1){ #On regarde si on a bien la possibilité de faire les tests suivants, il faut pour celà avoir au moins 2 souches
    if(any(as.integer(round(k[-nbsouche]*nbdek))==as.integer(round(k[nbsouche]*nbdek)))){#On regarde si le dernier k ajouté (le k mutant ) est égal à un k déja existant
      #Si c'est le cas on va alors enlever cette valeur
      k=k[-nbsouche] 
      X=X[-nbsouche]
      Y=Y[-nbsouche]
      for(i1 in 1:nbsouche){#Et on associe les valeurs de X et Y la résolution précédente aux valeurs présentes 
        Xfinal[k[i1]*(nbdek)+1,mut]=Xfinal[k[i1]*(nbdek)+1,mut]+X[i1]
        Yfinal[k[i1]*(nbdek)+1,mut]=Yfinal[k[i1]*(nbdek)+1,mut]+Y[i1]
      }
    }else{#Si le dernier k n'est pas déjà présent, on résoud normalement.
              #On definit les parametres des souches en presence#
      source("fonction_parametres.R")
      param_k=func_param(k)
      
              #On resoud le systeme#
      source("resolution_systemes_nequ_v2.R")
      resolution=resolution_systemes_nequ(param_k,X,Y,nbsouche)
      
      #La, on a des 0 si proche de 0, donc on peut savoir si on a disparition et stocker tout ca dans les matrices Xfinal et Yfinal
      
      #On supprime les k disparus 
      k_inter=NULL
      X_inter=NULL
      Y_inter=NULL
      for(j in 1:nbsouche){ #On recupere les couples resolution[1,2], puis resolution[3,4], etc..., qui representent les valeurs de X et Y pour une souche en particulier.
        if(resolution[2*j]>0 && resolution[2*j-1]>0){
          k_inter=c(k_inter,k[j])
          X_inter=c(X_inter,resolution[-1+2*j])
          Y_inter=c(Y_inter,resolution[2*j])
        }
      }
      if(is.null(Y_inter)){#On a ici extinction du parasite et donc equilibre non-trivial 
        write(paste("Disparition du parasite pour k=",kinitial,"avec les parametres beta=",param_k$betak,"alpha=",param_k$alphak,"gamma=",param_k$gammak,"a=",param_k$ak,"b=",param_k$bk,"et q=",param_k$qk,sep=" "),file=paste("Rapports_erreur/Error_for_k=",k,".txt",sep="")) #On écrit un message d'erreur...
        break #...et on sort de la boucle.
      }
      #Puis on remplace les vecteurs
      k=k_inter ; X=X_inter ; Y=Y_inter
      for(u in 1:length(k_inter)){
        Xfinal[round(k[u]*(nbdek))+1,mut]=Xfinal[round(k[u]*(nbdek))+1,mut]+X[u]
        Yfinal[round(k[u]*(nbdek))+1,mut]=Yfinal[round(k[u]*(nbdek))+1,mut]+Y[u]
        # print(paste("numeromut tour",u,":",round(k[u]*(nbdek))))
        # print(paste("X écrit tour",u,":",X[u]))
        # print(paste("Y écrit tour",u,":",Y[u]))
      }
    }
    
    #On a traité le cas où on fait une résolution à 2 souches ou plus. Il faut également traiter le cas où l'on a un seule souche (la première résolution)
    
  }else{#On definit les parametres des souches en presence#
    source("fonction_parametres.R")
    param_k=func_param(k)
    
    #On resoud le systeme#
    source("resolution_systemes_nequ.R")
    resolution=resolution_systemes_nequ(param_k,X,Y)
    
    #La, on a des 0 si proche de 0, donc on peut savoir si on a disparition et stocker tout ca dans les matrices Xfinal et Yfinal
    #On supprime les k disparus 
    k_inter=NULL
    X_inter=NULL
    Y_inter=NULL
    for(j1 in 1:(length(resolution)/2)){ #On recupere les couples resolution[1,2], puis resolution[3,4], etc..., qui representent les valeurs de X et Y pour une souche en particulier.
      if(resolution[2*j1]>0 && resolution[2*j1-1]>0){
        k_inter=c(k_inter,k[j1])
        X_inter=c(X_inter,resolution[2*j1-1])
        Y_inter=c(Y_inter,resolution[2*j1])
      }
    }
    if(is.null(Y_inter)){#On a ici extinction du parasite et donc equilibre non-trivial 
      write(paste("Disparition du parasite pour k=",kinitial,"avec les parametres beta=",param_k$betak,"alpha=",param_k$alphak,"gamma=",param_k$gammak,"a=",param_k$ak,"b=",param_k$bk,"et q=",param_k$qk,sep=" "),file=paste("Rapports_erreur/Error_for_k=",k,".txt",sep="")) #On écrit un message d'erreur...
      break #...et on sort de la boucle.
    }
    #Puis on remplace les vecteurs
    k=k_inter ; X=X_inter ; Y=Y_inter
    for(u in 1:length(k_inter)){
      Xfinal[round(k[u]*(nbdek))+1,mut]=Xfinal[round(k[u]*(nbdek))+1,mut]+X[u]
      Yfinal[round(k[u]*(nbdek))+1,mut]=Yfinal[round(k[u]*(nbdek))+1,mut]+Y[u]
    }
    #kfinal=rbind(kfinal,k)
    
  }
  #kfinal=rbind(kfinal,k) #On met les k résidents dans la matrice de vérification
  #Il faut maintenant passer à l'étape suivant la résolution du système.
  #On va ajouter le mutant
  
  # print(paste("X a lsa sortie :",X))
  # print(paste("Y a la sortie :",Y))
  
  
  vecmutant=(X+Y)/sum(X+Y)#vecteur de proba de muter en fonction de la densite relative de la population
  signe=sample(c(-1,1),1)
  ktire=sample(k,1,prob=vecmutant)
  if(ktire<2*1/nbdek){signe=1}
  if(ktire>1-1/nbdek){signe=-1} #Pour rester entre 0 et 1
  kmutant=ktire+signe*(1/nbdek)
  Xmutant=sum(X)*pourcentpopmutantX
  Ymutant=sum(Y)*pourcentpopmutantY
  # print(paste("kmutant=",kmutant))
  # print(paste("vec de k sortie sans mut :",k))
  #Et on fusionne ca
  k=c(k,kmutant) ; X=c(X,Xmutant) ; Y=c(Y,Ymutant)

}
print(Sys.time())
par(mfrow=c(1,3))
image(Xfinal+Yfinal,col=gray(rev(seq(from=0,to=1,by=0.01))))
image(Xfinal,col=gray(rev(seq(from=0,to=1,by=0.01))))
image(Yfinal,col=gray(rev(seq(from=0,to=1,by=0.01))))
# image(matrix(as.numeric(as.logical(Xfinal+Yfinal)),nrow=nbdek+1,ncol=nbmutant))
# image(matrix(as.numeric(as.logical(Xfinal)),nrow=nbdek+1,ncol=nbmutant))
# image(matrix(as.numeric(as.logical(Yfinal)),nrow=nbdek+1,ncol=nbmutant))
