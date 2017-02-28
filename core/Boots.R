#### corescript modèle Boots 2012 ####

# /!\ # Non fonctionnel
# /!\ # En chantier


## spécification des paramètres initiaux :

h <- 1 
Xi <- 1
Yi <- 0.1



## Appel des paramètres du modèle :

yini <- c(h, Xi, Yi)

# /!\ # appeler le fichier "compilé" via programme (non encore possible)

L <- fontion(h){
  return(alpha(h) + beta(h) + gamma(h))
}
# Il s'agit du paramètre traduisant de ce qui quite le groupe des Y, mis en fonction pour prendre en compte qu'il varie avec h
# Inutile de l'appeler d'un autre fichier car ne change normalement pas

q <- 1
# Paramètre de densité dépendance
# Valeur fixe ?



## Résolution numérique :

# /!\ # C'est le format de l'ancien fichier donc tout est à modifier même si ça permet d'avoir le squelette


duree <- 30 # temps que dure la simulation
t <- seq(0, duree, by=0.1) # découpage des pas de temps qu'on veut sur le graph



hote_h <- function(yini, param_h){
  Xh <- yini[2]
  Yh <- yini[3]
  H <- sum(yini[2:3])
  with(as.list(param_h), {
    
    dXh <- a(h)*Xh - q*H*Xh - b*Xh - sum(betah*Xh*Yh) + gamma*Yh
    dYh <- sum(betah*Xh*Yh) - L*Yh
    dH <- dXh + dYh
    
    return(list(c(dXh, dYh, dH)))
  })
}

out <- ode(y=yini,times=t ,func=hote_h, parms=param)


### Repr?sentation graphique :

par(mfrow=c(1, 3))


plot(out[, 1], out[, 2], type="l")
plot(out[, 1], out[, 3], type="l")
plot(out[, 2], out[, 3], type="l")
