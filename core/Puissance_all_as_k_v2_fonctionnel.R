#Resolution "optimisee" du systeme a n equations en prenant en compte une expression de chaque parametre en fonction de k selon notre expression en puissance.
rm(list = ls())
library(deSolve)
setwd("~/Documents/Modeles_hotes") #On set le working directory a la racine de la serie de fichier

source("./parametres/valeurs_parametres_puissance.R") #Contient les variables definissant les TO.
source("./parametres/expression_parametres_puissance.R") #Contient les fonctions exprimant les TO.
k=0.1
variationk=0.01
nbsouchemax=15 
duree=25000 #Durée de la résolution numérique via ode()
nombre_mutations=1000 #Nombre de mutations effectuées


### Initialisation des vecteurs de stockage de résultats
Xfinal=matrix(0,nombre_mutations,nbsouchemax)
Yfinal=matrix(0,nombre_mutations,nbsouchemax)
kfinal=matrix(0,nombre_mutations,nbsouchemax)
### Initialisation des paramètres ###

Xi=1
Yi=0.1

                      ### Resolution ###
compt_couple=c(-1,0) #Sert pour la sélection des futurs résidents. Nécessairement défini comme içi.
vecteurini=c(Xi,Yi)
#Le vecteur vecteurini doit comprendre les valeurs de X et Y des résidents et du mutant de la prochaine itération. Première fois : uniquement celles du premier résident.
for(i in 1:nombre_mutations){
  number_strain=(length(vecteurini)/2)
  
  #beta_k#
  betak=beta(beta0,beta_max,P_beta,k)
  
  # a_k : Definition #
  ak=a(a0,a_max,P_a,k)
  
  # alpha_k #
  alphak=alpha(alpha0,alpha_max,P_alpha,k)
  
  # gamma_k #
  gammak=gamma(gamma0,gamma_max,P_gamma,k)
  
  # b_k #
  bk=b(b0,b_max,P_b,k)
  
  # q_k #
  qk=cu(ak,bk,K)
  
  param_k=list(betak=betak,ak=ak,alphak=alphak,gammak=gammak,bk=bk,qk=qk)
  source(paste("sys_equations/data_", number_strain, sep=""), local=T)
  out <- ode(y=vecteurini,times=c(0, seq(duree-0.1, duree, 0.01)) ,func=hote_h, parms=param_k)
  out=out[,-1]
  if(sum((out[11,1]-out[2:10,1])^2)>10^-4){
    sink("error_report.txt",append=TRUE)
    print(paste("La population residente n'etait pas a l'equilibre a l'iteration",i,". La population avait pour parametre",k[i],".   "))
    sink()
  } #Va remplir un fichier d'erreur si on n'atteint pas l'équilibre
  
  #Le vecteur out aura 2n+1 sorties. Il va sortir les résultats qui nous intéressent par couple au vu de la structure du fichier de création du système d'équations. Du coup, on va vouloir sélectionner les couples de valeurs (Xk,Yk) dont au moins un des deux membres est non-nul. 
  vecteurini=NULL
  X=rep(0,nbsouchemax)
  Y=rep(0,nbsouchemax)
  kinter=rep(0,nbsouchemax)
  for(j in 1:((ncol(out)/2))){
    couple=out[11,compt_couple+2*j]
    if(couple[1]>10^-2&&couple[2]>10^-2){
      X[j]=couple[1]
      Y[j]=couple[2]
      kinter[j]=k[j]
      vecteurini=c(vecteurini,couple[1],couple[2])
    }
  }
  k=kinter #Il est nécessaire de différencier k et kinter puisqu'on nomme l'un en fonction de l'autre juste avant
  
  kfinal[i,]=k #on stocke les valeurs de k encore présente
  Yfinal[i,]=Y
  Xfinal[i,]=X
  #On va maintenant intégrer le mutant : il faut déterminer son k et ses paramètres initiaux.
  vecmutant=X+Y/sum(X+Y)#vecteur de proba de muter en fonction de la densité relative de la population
  signe=sample(c(-1,1),1)
  kmutant=sample(k,1,prob=vecmutant)+signe*(variationk)
  Xmutant=sum(X)/10
  Ymutant=sum(Y)/10
  for(i in 1:length(k)) if(k[i]==0) {k[i]=kmutant;break}
  vecteurini=c(vecteurini,Xmutant,Ymutant)
}
