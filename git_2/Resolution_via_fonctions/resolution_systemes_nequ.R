resolution_systemes_nequ=function(param_k,X,Y){
  #On crée le vecteur de valeurs initiales pour la résolution
 V0=NULL
 for (i in 1:length(X)) {V0 <- c(V0, X[i], Y[i])}
  
  resolmax=200
  source(paste("sys_equations/data_", length(param_k$ak), sep=""), local=T)
  final=ode(y=V0, t=c(0,seq(resolmax-1,resolmax,0.1)), func=hote_h, parms=param_k)
          
  #Test de l'équilibre : si on n'y est pas, il faut recommencer la résolution de là où on en était 
  compteurtourmax=0 #Si on dépasse un certain nombre de répétitions, ça viens peut-être d'un mvt cyclique, auquel cas il faut s'arrêter.
  verif=0 #On initialise le compteur de vérification d'équilibre...
  for(i in 1:length(param_k$ak)){
    verif=verif+sum(abs(final[11,i+1]-final[2:10,i+1]))#...et on le calcule.
  }

  #Si on n'a pas atteint l'équilibre et qu'on n'a pas atteint notre "seuil de cyclicité", i faut continuer la résolution.
  while(verif>10^-4 && compteurtourmax<10){
    final=ode(y=final[nrow(final),-1], t=c(resolmax,seq(resolmax*2-1,resolmax*2,0.1)), func=hote_h, parms=param_k)
    resolmax=resolmax*2
    compteurtourmax=compteurtourmax+1
    verif=0

    for(i in 1:length(param_k$ak)){#On recalcule la vérification pour savoir si on refait un tour.
      verif=verif+sum(abs(final[11,i]-final[2:10,i]))
    }

  }
  
  #Et au final, on renvoie le résultat de notre système, donné par la dernière ligne de "final".

  return(round(final[nrow(final),-1],digits=2))#On arrondi à la deuxième décimale pour avoir des 0 si une équation tends vers 0.
}
