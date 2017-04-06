resolution_systemes_nequ=function(param_k,X,Y,nbsouche,resolmax,arrondi){

	##### INITIALISATION #####

#times <- c(0,seq(resolmax-1,resolmax,0.1))
times <- c(0,resolmax)
#zero <- 1e-4	
#compteurtourmax <-10
#times_cycle_sup <- c(resolmax,seq(resolmax*2-1,resolmax*2,0.1))



  #On crée le vecteur de valeurs initiales pour la résolution
 V0=NULL
 for (i in 1:length(X)) {V0 <- c(V0, X[i], Y[i])}
  
 
  source(paste("sys_equations/data_", length(param_k$ak), sep=""), local=T)
  final=ode(y=V0, t=times, func=hote_h, parms=param_k)
          
  #Test de l'équilibre : si on n'y est pas, il faut recommencer la résolution de là où on en était 
 #  compteurtour=0 #Si on dépasse un certain nombre de répétitions, ça viens peut-être d'un mvt cyclique, auquel cas il faut s'arrêter.
 # # verif=0 #On initialise le compteur de vérification d'équilibre...
 #  for(i in 1:length(param_k$ak)){
 #    verif=verif+sum(abs(final[nrow(final),i+1]-final[-c(1, nrow(final)),i+1]))#...et on le calcule.
 #  }
 # 
 #  #Si on n'a pas atteint l'équilibre et qu'on n'a pas atteint notre "seuil de cyclicité", i faut continuer la résolution.
 #  while(verif>zero && compteurtour<compteurtourmax){
 #    final=ode(y=final[nrow(final),-1], t=times_cycle_sup, func=hote_h, parms=param_k)
 #    resolmax=resolmax*2
 #    compteurtour=compteurtour+1
 #    times_cycle_sup <- c(resolmax,seq(resolmax*2-1,resolmax*2,0.1))
 #    verif=0
 # 
 #    for(i in 1:length(param_k$ak)){#On recalcule la vérification pour savoir si on refait un tour.
 #      verif=verif+sum(abs(final[nrow(final),i]-final[-c(1, nrow(final)),i]))
 #    }
 # 
 #  }
  
  #Et au final, on renvoie le résultat de notre système, donné par la dernière ligne de "final".

  return(round(final[nrow(final),-1],digits=arrondi))
}
