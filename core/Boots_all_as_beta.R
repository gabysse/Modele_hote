#Résolution "optimisée" du système à 2 équations selon la méthode de Runge-Kutta en prenant en compte une expression de chaque paramètre en fonction de beta selon l'article de Boots 2012
rm(list = ls())

setwd("~/Documents/Modeles_hotes") #On set le working directory à la racine de la 

                      ### Paramètres initiaux ###
k=1 #k peut varier sans problème si nécessaire. 
Xi <- 1
Yi <- 0.1
pop=Xi+Yi
H <- sum(pop)

betah <- 5-k/5
q <- 0.0001

# a_h : Définition des paramètres nécessaires #
source("./parametres/a/a_boots.R") #Va définir la fonction d'expression de a en fonction de beta exprimé dans le fichier a_boots.R
a1=4
a2=-0.84
a3=-4.98
a4=-1
ah=a(a1,a2,a3,a4,betah)

# alpha_h : Définition des paramètres #
source("./parametres/alpha/alpha_boots.R") 
alpha1=0.1
alpha2=-0.05
alpha3=-10
alpha4=-1
alphah=alpha(a1,a2,a3,a4,betah)

# gamma_h : Définition #
source("./parametres/gamma/gamma_boots.R")
gamma1=4
gamma2=-0.5
gamma3=2
gamma4=-1
gammah=gamma(gamma1,gamma2,gamma3,gamma4,betah)

# b_h : Définition #
source("./parametres/b/b_boots.R")
b1=4
b2=-0.5
b3=2
b4=-1
bh=b(b1,b2,b3,b4,betah)

                  ### Résolution via Runge-Kutta ###

#On va tenter de résoudre ça simplement avec un Runge-Kutta simple qui permettrait d'être optimisé sur la résolution en limitant le temps de calcul au minimum en étant sûr d'atteindre l'équilibre.
h=0.05
i=1
X=Xi
Y=Yi
Xfin=Xi
Yfin=Yi
testX=0
testY=0
while(i<20|(testX>0.00001|testY>0.00001)){
  Xint=X+h/2*(ah*X-q*H*X-bh*X-betah*X*Y+gammah*Y)
  Yint=Y+h/2*(betah*X*Y-Y*(alphah+bh+gammah))
  X=X+h*(ah*Xint-q*H*Xint-bh*Xint-betah*Xint*Yint+gammah*Yint)
  Y=Y+h*(betah*Xint*Yint-Yint*(alphah+bh+gammah))
  Xfin=c(Xfin,X)
  Yfin=c(Yfin,Y)
  
  #On calcule la condition 
  if(i<15){
    testX=0;testY=0
  }else{
    testX=sum((X-Xfin[(i-10):i])^2);testY=sum((Y-Yfin[(i-10):i])^2)
  }
  
  H=sum(X,Y)#H est un paramètre qui évolue avec la population, il est nécessaire de le recalculer à chaque pas de temps.
  
  i=i+1 #On incrémente le compteur
}

par(mfrow=c(1,2))
plot(1:length(Xfin),Xfin,t="l")
plot(1:length(Yfin),Yfin,t="l")

