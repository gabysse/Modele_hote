#Résolution "optimisée" du système à 2 équations selon la méthode de Runge-Kutta en prenant en compte une expression de chaque paramètre en fonction de k selon notre expression en puissance.
rm(list = ls())

setwd("~/Documents/Modeles_hotes") #On set le working directory à la racine de la série de fichier

### Paramètres initiaux ###
k=0.000001
Xi <- 1
Yi <- 0.1
pop=Xi+Yi
H <- sum(pop)
q <- 1

#beta_h#
source("./parametres/beta/beta_puissance.R")
beta0=0
beta_max=2
P_beta=1
betah=beta(beta0,beta_max,P_beta,k)

# a_h : Définition des paramètres nécessaires #
source("./parametres/a/a_puissance.R") #Va définir la fonction d'expression de a en fonction de beta exprimé dans le fichier a_boots.R
a0=0
a_max=10
P_a=2
ah=a(a0,a_max,P_a,k)
#La forme des trades-off ne change pas, on garde donc les paramètres définisant le trade-off pour la valeur des paramètres du mutant.

# alpha_h : Définition des paramètres #
source("./parametres/alpha/alpha_puissance.R") 
alpha0=0
alpha_max=5
P_alpha=2
alphah=alpha(alpha0,alpha_max,P_alpha,k)

# gamma_h : Définition #
source("./parametres/gamma/gamma_puissance.R")
gamma0=0
gamma_max=1
P_gamma=0.5
gammah=gamma(gamma0,gamma_max,P_gamma,k)

# b_h : Définition #
source("./parametres/b/b_puissance.R")
b0=0.2
b_max=1
P_b=0.5
bh=b(b0,b_max,P_b,k)

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
while(i<20|(testX>10^-8|testY>10^-8)){
  Xint=X+h/2*(ah*X-q*H*X-bh*X-betah*X*Y+gammah*Y)
  Yint=Y+h/2*(betah*X*Y-Y*(alphah+bh+gammah))
  X=X+h*(ah*Xint-q*H*Xint-bh*Xint-betah*Xint*Yint+gammah*Yint)
  Y=Y+h*(betah*Xint*Yint-Yint*(alphah+bh+gammah))
  Xfin=c(Xfin,X)
  Yfin=c(Yfin,Y)
  
  #On calcule les conditions de ressemblance des dernières valeurs
  if(i<15){
    testX=0;testY=0 #Ainsi, on a pas de calcul pendant les 15 premières itérations le temps de remplir les vecteurs Xfin et Yfin. De toute façon, on a la condition de la boucle qui ne permet pas d'avoir moins de 20 tours de boucle. 
  }else{
    testX=sum((X-Xfin[(i-10):i])^2);testY=sum((Y-Yfin[(i-10):i])^2)
  }
  
  H=sum(X,Y)#H est un paramètre qui évolue avec la population, il est nécessaire de le recalculer à chaque pas de temps.
  
  i=i+1 #On incrémente le compteur
}

par(mfrow=c(1,2))
plot(1:length(Xfin),Xfin,t="l",main="Susceptible hosts", xlab="Time",ylab="Density")
plot(1:length(Yfin),Yfin,t="l",main="Infected hosts", xlab="Time",ylab="Density")
