#Résolution "optimisée" du système à 2 équations selon la méthode de Runge-Kunta en prenant en compte une expression de a et alpha en fonction de beta selon l'article de Boots 2012

setwd("~/Documents/Modeles_hotes") #On set le working directory à la racine de la 

### Paramètres initiaux
h <- 1 
Xi <- 1
Yi <- 0.1
pop = c(Xi, Yi)

source("./parametres/a/a_boots.R") #Va définir la fonction d'expression de a en fonction de beta exprimé dans le fichier a_boots.R
source("./parametres/alpha/alpha_boots.R") 
betah <- 5 - h/5
b <- 1
q <- 1
gamma <- 1
H <- sum(pop)


ah=a(betah)
alphah=alpha(betah)
L=alphah + betah +gamma
#Résolution via Runge-Kutta

duree <- 5 # temps que dure la simulation
ta <- seq(1, 30, by=1/50) # d?coupage des pas de temps

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
  Xint=X+h/2*(ah*X-q*H*X-b*X-betah*X*Y+gamma*Y)
  Yint=Y+h/2*(betah*X*Y-Y*(alphah+b+gamma))
  X=X+h*(ah*Xint-q*H*Xint-b*Xint-betah*Xint*Yint+gamma*Yint)
  Y=Y+h*(betah*Xint*Yint-Yint*(alphah+b+gamma))
  Xfin=c(Xfin,X)
  Yfin=c(Yfin,Y)
  if(i<15){
    testX=0;testY=0
  }else{
    testX=sum((X-Xfin[(i-10):i])^2);testY=sum((Y-Yfin[(i-10):i])^2)
  }
  i=i+1
}
par(mfrow=c(2,2))
plot(1:length(Xfin),Xfin,t="l")
plot(1:length(Yfin),Yfin,t="l")
plot(Xfin,Yfin,t="l")
