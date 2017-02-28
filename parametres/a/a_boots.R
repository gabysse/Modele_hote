#### Paramètre de fécondité, "a", trade-off en fonction de beta "à la Boots12" ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant

a <- function(a1,a2,a3,a4,beta){
  return(a1+a2*(exp(a3*(beta-a1))+a4))
}

