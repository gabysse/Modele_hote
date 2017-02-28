### fonction permettant d'implémenter un nos systèmes d'équations ###

# nécessité d'un vecteur k contenant nos différentes valeurs de k

setwd("~/Documents/Modeles_hotes")


# Paramètres pour réaliser des vérifications
k <- c(0.5, 1)
X <- c(1, 1.3)
Y <- c(0.1, 1)
pop <- c(X, Y)

a <- 4.8
b <- 1
q <- 1
alpha <- 1
gamma <- 1
beta <- function(k){return(5-k/5)}
L <- alpha+beta(k)+gamma
H <- sum(X+Y)
param_k <- c(a, b,q, alpha, beta, gamma,L, H)


# création d'un jeu d'equation
equa_diff <- function(i){
  return(paste(
    paste("dX_",i,"<- a*X[",i,"] - q*H*X[",i,"] - b*X[",i,"] - sum(beta(k[",i,"])*X[",i,"]*Y[",i,"]) + gamma*Y[",i,"];", sep=""),
    
    paste("dY_",i,"<- sum(beta(k[",i,"])*X[",i,"]*Y[",i,"]) - L[",i,"]*Y[",i,"];", sep=""))
  )
}


# création jeux deonnées sous forme dX_...
equa_diff <- function(i){
  return(paste(
    paste("dX_",i,"<- a*X_",i," - q*H*X_",i," - b*X_",i," - sum(beta(k[",i,"])*X_",i,"*Y_",i,") + gamma*Y_",i,";", sep=""),
    
    paste("dY_",i,"<- sum(beta(k[",i,"])*X_",i,"*Y_",i,") - L[",i,"]*Y_",i,";", sep=""))
  )
}


equa_diff(2)
a*X[1] - q*H*X[1] - b*X[1] - sum(beta(k[1])*X[1]*Y[1]) + gamma*Y[1]
sum(beta(k[1])*X[1]*Y[1]) - L[1]*Y[1]
a*X[2] - q*H*X[2] - b*X[2] - sum(beta(k[2])*X[2]*Y[2]) + gamma*Y[2]
sum(beta(k[2])*X[2]*Y[2]) - L[2]*Y[2]

is.character(equa_diff(1))

n<-15
for(i in 1:n){
  for(j in 1:i){
   #write(equa_diff(j), file = paste("data_", i, sep=""), append=T)
  }
}


# Pour en créer plusieurs (dépendant du vecteur k)
for(i in 1:length(k)){
  equa_diff(i)
  print(equa_diff(i))# créé système Xk et Yk pour chaque valeur de k
}



# obtenir une somme avec une boucle for à l'intérieur
somme <- 0
dH <-0
dH_<-  for(i in 1:length(k)){
  dH <- dH + sum(X[i], Y[i])
}
dH
sum(somme)


## Fonction finale à faire tournée :
hote_h <- function(t,pop , param_k){
  
  # initialisation des paramètres
  #for(i in 1:length(k)){
   # X[i] <- X[i]
    #Y[i] <- Y[i]
  #}
  X_1 <- X[1]
  X_2 <- X[2]
  
  Y_1 <- Y[1]
  Y_2 <- Y[2]
  H <- sum(pop)
  
  # début des équa diff
  with(as.list(param_k), {
    
    source(paste("data_", length(k), sep=""), local=T)
    dH <- dX_1 + dY_1 + dX_2 + dY_2
    
    return(list(c(dX_1, dY_1, dX_2, dY_2, dH))) # pose le même soucis d'appel que pour dH
  })
}
out <- ode(y=c(X[1], Y[1], X[2], Y[2], H),times=seq(0,20,0.1) ,func=hote_h, parms=param_k)
head(out)
tail(out)
diagnostics(out)

plot(out[,1], out[,2])
plot(out[,1], out[,3])
plot(out[,2], out[,3])

is.character(paste(1))
paste("data_", paste(length(k)), sep="")
")))
