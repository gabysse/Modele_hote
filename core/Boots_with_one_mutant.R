#Résolution "optimisée" du système à 2 équations selon la méthode de Runge-Kutta en prenant en compte une expression de chaque paramètre en fonction de beta selon l'article de Boots 2012 et introduction d'un mutant
rm(list = ls())

setwd("~/Documents/Modeles_hotes") #On set le working directory à la racine de la ésrie de fichier

### Paramètres initiaux pour le futur résident ###
k=1 

Xi <- 1
Yi <- 0.1
pop=Xi+Yi
H <- sum(pop)
betah <- 5-k/5
q <- 1

km=k+0.5#+rnorm(1,0,0.1) # On a une variation de k autour du k résident
Xim <- 0.1 
Yim <- 0
betahm <- 5-km/5 #Beta ressemble toujours au beta précédent mais le km sera donc changé.

# a_h : Définition des paramètres nécessaires #
source("./parametres/a/a_boots.R") #Va définir la fonction d'expression de a en fonction de beta exprimé dans le fichier a_boots.R
a1=4
a2=-0.84
a3=-4.98
a4=-1
ah=a(a1,a2,a3,a4,betah)
ahm=a(a1,a2,a3,a4,betahm)
#La forme des trades-off ne change pas, on garde donc les paramètres définisant le trade-off pour la valeur des paramètres du mutant.

# alpha_h : Définition des paramètres #
source("./parametres/alpha/alpha_boots.R") 
alpha1=0.1
alpha2=-0.05
alpha3=-10
alpha4=-1
alphah=alpha(a1,a2,a3,a4,betah)
alphahm=alpha(a1,a2,a3,a4,betahm)

# gamma_h : Définition #
source("./parametres/gamma/gamma_boots.R")
gamma1=4
gamma2=-0.5
gamma3=2
gamma4=-1
gammah=gamma(gamma1,gamma2,gamma3,gamma4,betah)
gammahm=gamma(gamma1,gamma2,gamma3,gamma4,betahm)

# b_h : Définition #
source("./parametres/b/b_boots.R")
b1=4
b2=-0.5
b3=2
b4=-1
bh=b(b1,b2,b3,b4,betah)
bhm=b(b1,b2,b3,b4,betahm)

### Résolution via Runge-Kutta ###

#On va tenter de résoudre ça simplement avec un Runge-Kutta simple qui permettrait d'être optimisé sur la résolution en limitant le temps de calcul au minimum en étant sûr d'atteindre l'équilibre.
h=0.05
i=1
j=1
X=Xi
Y=Yi
Xfin=rep(0,10)
Yfin=rep(0,10)
testX=0
testY=0
while(i<20|(testX>0.0000000001|testY>0.0000000001)){
  Xint=X+h/2*(ah*X-q*H*X-bh*X-betah*X*Y+gammah*Y)
  Yint=Y+h/2*(betah*X*Y-Y*(alphah+bh+gammah))
  X=X+h*(ah*Xint-q*H*Xint-bh*Xint-betah*Xint*Yint+gammah*Yint)
  Y=Y+h*(betah*Xint*Yint-Yint*(alphah+bh+gammah))
  Xfin[j]=X
  Yfin[j]=Y
  #On modifie j entre 1 et 10. L'ordre importe peu, ce qui importe c'est d'avoir les 10 dernières valeurs de Y et X pour calculer les tests d'équilibre. 
  if(j==10){j=1
  }else{j=j+1
  }
  #On calcule la condition 
  if(i<15){
    testX=0;testY=0
  }else{
    testX=sum((X-Xfin[1:10])^2);testY=sum((Y-Yfin[1:10])^2)
    #On n'aura pas de soucis de calcul puisque nos Xfin et Yfin se rempliront avant que la condition i>=15 ne soit atteinte
  }
  
  H=sum(X,Y)#H est un paramètre qui évolue avec la population, il est nécessaire de le recalculer à chaque pas de temps.
  
  i=i+1 #On incrémente le compteur
}

                 ### Arrivée du mutant ###
pop=X+Y+Xim+Yim
H <- sum(pop)
#Le paramètre q ne dépends pas de h pour nous dans un premier temps, on le laisse donc tel quel.

        # Résolution selon la méthode de Runge-Kutta
h=0.05
i=1
Xm=Xim
Ym=Yim
#Attention à ne pas réinitialiser X et Y : on se situe à l'équilibre ! L'interaction va se faire au niveau de H, dont l'expression change.
Xfin=X
Yfin=Y
Xmfin=Xim
Ymfin=Yim
testX=0
testY=0
testXm=0
testYm=0
#On envoie la boucle
while(i<20|(testX>0.0000000001|testY>0.0000000001|testXm>0.0000000001|testYm>0.0000000001)){
      #Définition des valeurs au demi-pas
  Xint=X+h/2*(ah*X-q*H*X-bh*X-betah*X*(Ym+Y)+gammah*Y) #Attention, l'expression de l'infectivité change puisqu'on a désormais 2 types d'infection possible selon le type d'hôte.
  Yint=Y+h/2*(betah*X*(Y+Ym)-Y*(alphah+bh+gammah))
  Xmint=Xm+h/2*(ahm*Xm-q*H*Xm-bhm*Xm-betahm*Xm*(Y+Ym)+gammahm*Ym)
  Ymint=Ym+h/2*(betahm*Xm*(Ym+Y)-Ym*(alphahm+bhm+gammahm))
  
      #Application pour le calcul au pas
  X=X+h*(ah*Xint-q*H*Xint-bh*Xint-betah*Xint*(Ymint+Yint)+gammah*Yint)
  Y=Y+h*(betah*Xint*(Ymint+Yint)-Yint*(alphah+bh+gammah))
  Xm=Xm+h*(ahm*Xmint-q*H*Xmint-bhm*Xmint-betahm*Xmint*(Yint+Ymint)+gammahm*Ymint)
  Ym=Ym+h*(betahm*Xmint*(Ymint*Yint)-Ymint*(alphahm+bhm+gammahm))
      
  
  Xfin=c(Xfin,X)
  Yfin=c(Yfin,Y)
  Xmfin=c(Xmfin,Xm)
  Ymfin=c(Ymfin,Ym) #On pourra changer tout ça quand on sera convaincu que ça marche pour éviter de créer de trop gros vecteurs (cf calcul du R-K pour le résident seul), et ne garder que la dernière valeur.
  
  #On calcule la condition 
  if(i<15){
    testX=0;testY=0;testXm=0;testYm=0
  }else{
    testX=sum((X-Xfin[(i-10):i])^2);testY=sum((Y-Yfin[(i-10):i])^2)
    testXm=sum((Xm-Xmfin[(i-10):i])^2);testYm=sum((Ym-Ymfin[(i-10):i])^2)
  }
  
  H=sum(X,Y,Xm,Ym)#H est un paramètre qui évolue avec la population, il est nécessaire de le recalculer à chaque pas de temps.
  
  i=i+1 #On incrémente le compteur
}

        ##Représentation graphique##

par(mfrow=c(2,2))
plot(1:length(Xfin),Xfin,t="l",main="Résident sain")
plot(1:length(Yfin),Yfin,t="l",main="Résident infecté")
plot(1:length(Xmfin),Xmfin,t="l",main="Mutant sain")
plot(1:length(Ymfin),Ymfin,t="l",main="Mutant infecté")
