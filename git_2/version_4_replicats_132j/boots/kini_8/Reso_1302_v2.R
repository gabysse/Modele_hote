

resolution_multiple=function(nbdek,nbmutant,Xini,Yini,timebfmut,dec_arrondi,Ymutant,Xmutant,pos_k_ini){
	Y=Yini
	X=Xini
	nbdek=nbdek
	variationk=1/nbdek
	nbmutant=nbmutant
	Xfinal=matrix(0 , nrow = nbdek+1, ncol = nbmutant)
	Yfinal=matrix(0 , nrow = nbdek+1, ncol = nbmutant)
	resolmax=timebfmut
	arrondi=dec_arrondi
	Xmutant=Xmutant*10^(-arrondi)
	Ymutant=Ymutant*10^(-arrondi)
	posrelativek=pos_k_ini
	k=(posrelativek-1)*variationk
	kinitial=k
	  
	  
	for(mut in 1:nbmutant){
	  nbsouche=length(k)
	  # print(paste("mut",mut))
	  # print(paste("nbsouche=",nbsouche))
	  # print(paste("vec de k entree :",k))
	  # print(paste("X a lentree :",X))
	  # print(paste("Y a lentree :",Y))
	  if(nbsouche>1){ #On regarde si on a bien la possibilité de faire les tests suivants, il faut pour celà avoir au moins 2 souches
	    if(sum(posrelativek[nbsouche]==posrelativek[-nbsouche])>0){#On regarde si le dernier k ajouté (le k mutant ) est égal à un k déja existant
	      #Si c'est le cas on va alors enlever cette valeur
	      k=k[-nbsouche] 
	      posrelativek=(k/variationk)+1
	      X=X[-nbsouche]
	      Y=Y[-nbsouche]
	      nbsouche=nbsouche-1
	      Xfinal[posrelativek,mut]=Xfinal[posrelativek,mut]+X
	      Yfinal[posrelativek,mut]=Yfinal[posrelativek,mut]+Y
	      
	      
	    }else{#Si le dernier k n'est pas déjà présent, on résoud normalement.
	      #On definit les parametres des souches en presence#
	      #source("fonction_parametres.R")
	      param_k=func_param(k)
	      
	      #On resoud le systeme#
	      #source("resolution_systemes_nequ_v2.R")
	      resolution=resolution_systemes_nequ(param_k,X,Y,nbsouche,resolmax,arrondi)
	      
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
	      posrelativek=round((k/variationk)+1)
	      Xfinal[posrelativek,mut]=Xfinal[posrelativek,mut]+X
	      Yfinal[posrelativek,mut]=Yfinal[posrelativek,mut]+Y
	    }
	    
	    
	    #On a traité le cas où on fait une résolution à 2 souches ou plus. Il faut également traiter le cas où l'on a un seule souche (la première résolution)
	    
	  }else{#On definit les parametres des souches en presence#
	    #source("fonction_parametres.R")
	    param_k=func_param(k)
	    
	    #On resoud le systeme#
	    #source("resolution_systemes_nequ_v2.R")
	    resolution=resolution_systemes_nequ(param_k,X,Y,nbsouche,resolmax,arrondi)
	    
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
	    posrelativek=round((k/variationk)+1)
	    Xfinal[posrelativek,mut]=Xfinal[posrelativek,mut]+X
	    Yfinal[posrelativek,mut]=Yfinal[posrelativek,mut]+Y
	    
	    #kfinal=rbind(kfinal,k)
	    
	  }
	  #kfinal=rbind(kfinal,k) #On met les k résidents dans la matrice de vérification
	  #Il faut maintenant passer à l'étape suivant la résolution du système.
	  #On va ajouter le mutant
	  
	  # print(paste("X a lsa sortie :",X))
	  # print(paste("Y a la sortie :",Y))
	  
	  #print(paste("Simulation terminée à",mut/nbmutant*100,"%"))
	  vecmutant=(X+Y)/sum(X+Y)#vecteur de proba de muter en fonction de la densite relative de la population
	  signe=sample(c(-1,1),1)
	  ktire=sample(k,1,prob=vecmutant)
	  if(ktire<2*variationk){signe=1}
	  if(ktire>1-variationk){signe=-1} #Pour rester entre 0 et 1
	  kmutant=ktire+signe*(variationk)
	  # print(paste("kmutant=",kmutant))
	  # print(paste("vec de k sortie sans mut :",k))
	  #Et on fusionne ca
	  k=c(k,kmutant) ; X=c(X,Xmutant) ; Y=c(Y,Ymutant)
	  posrelativek=round((k/variationk)+1)
	}

	resultats=list(Xfinal,Yfinal)

	return(resultats)
}
