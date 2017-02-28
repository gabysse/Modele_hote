rm(list=ls())
k=seq(0.001,5,0.001)
betah <- 5-k/5
par(mfrow=c(2,3))
# a_h : Définition des paramètres nécessaires #
source("./parametres/a/a_boots.R") #Va définir la fonction d'expression de a en fonction de beta exprimé dans le fichier a_boots.R
a1=4
a2=-0.84
a3=-4.98
a4=-1
ah=a(a1,a2,a3,a4,betah)
plot(k,ah,t="l",main="a") 

# alpha_h : Définition des paramètres #
source("./parametres/alpha/alpha_boots.R") 
alpha1=0.1
alpha2=-0.05
alpha3=-10
alpha4=-1
alphah=alpha(a1,a2,a3,a4,betah)
plot(k,alphah,t="l",main="alpha")

# gamma_h : Définition #
source("./parametres/gamma/gamma_boots.R")
gamma1=4
gamma2=-0.5
gamma3=2
gamma4=-1
gammah=gamma(gamma1,gamma2,gamma3,gamma4,betah)
plot(k,gammah,t="l",main="gamma")

# b_h : Définition #
source("./parametres/b/b_boots.R")
b1=4
b2=-0.5
b3=2
b4=-1
bh=b(b1,b2,b3,b4,betah)
plot(k,bh,t="l",main="b")

plot(k,betah,t="l",main="beta")