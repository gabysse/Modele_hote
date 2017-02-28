#Resolution "optimisee" du systeme a n equations en prenant en compte une expression de chaque parametre en fonction de k selon notre expression en puissance.
rm(list = ls())
library(deSolve)
setwd("~/Documents/Modeles_hotes") #On set le working directory a la racine de la serie de fichier

source("./parametres/valeurs_parametres_puissance.R") #Contient les variables definissant les TO.
source("./parametres/expression_parametres_puissance.R") #Contient les fonctions exprimant les TO.
k=0.5
variationk=0.01
nbsouchemax=15

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
qk=q(ak,bk,K)

### Initialisation des vecteurs de stockage de résultats
Xfinal=NULL
Yfinal=NULL
kfinal=NULL
### Parametres initiaux ###

Xi <- 1
Yi <- 0.1
pop=Xi+Yi
H <- sum(pop)
                      ### Resolution ###
compt_couple=c(-1,0) #Sert pour la sélection des futurs résidents
duree=10
nombre_mutations=100
vecteurini=c(Xi,Yi,H)
#Le vecteur vecteurini doit comprendre les valeurs de X et Y des résidents et du mutant de la prochaine itération
for(i in 1:nombre_mutations){
  number_strain=(lenght(vecteurini)/2)-1
  out <- ode(y=c(vecteurini),times=c(0, seq(duree-0.1, duree, 0.01)) ,func="hote_h", parms=param_h)
  if(sum((out[11]-out[2:10])^2)>10^-4){
    sink("error_report.txt",append=TRUE)
    print(paste("La population residente n'etait pas a l'equilibre a l'iteration",i,". La population avait pour parametre",k[i],".   "))
    sink()
  } #Va remplir un fichier d'erreur si on n'atteint pas l'équilibre 
  #Le vecteur out aura 2n+1 sorties. Il va sortir les résultats qui nous intéressent par couple au vu de la structure du fichier de création du système d'équations. Du coup, on va vouloir sélectionner les couples de valeurs (Xk,Yk) dont au moins un des deux membres est non-nul. 
  vecteurini=NULL
  X=rep(0,nbsouchemax)
  Y=rep(0,nbsouchemax)
  kinter=rep(0,nbsouchemax)
  for(j in 1:((lenght(out)/2)-1)){
    couple=out[compt_couple+2*j]
    if(couple[1]>10^-5&&couple[2]>10^-5){
      X[j]=couple[1]
      Y[j]=couple[2]
      kinter[j]=k[j]
      vecteurini=c(vecteurini,couple[1],couple[2])
    }
  }
  k=kinter #Il est nécessaire de différencier k et kinter puisqu'on nomme l'un en fonction de l'autre juste avant
  
  kfinal=rbind(kfinal,k) #on stocke les valeurs de k encore présente
  Yfinal=rbind(Yfinal,Y)
  Xfinal=rbind(Xfinal,X)
  #On va maintenant intégrer le mutant : il faut déterminer son k et ses paramètres initiaux.
  vecmutant=X+Y/sum(X+Y)#vecteur de proba de muter en fonction de la densité relative de la population
  signe=sample(c(-1,1),1)
  kmutant=sample(k,1,prob=vecmutant)+signe*(variationk)
  Xmutant=sum(X)/1000
  Ymutant=sum(Y)/1000
  X=c(X,Xmutant)
  Y=c(Y,Ymutant)
  H=sum(X,Y)
  
  vecteurini=c(vecteurini,Xmutant,Ymutant,H)
}
