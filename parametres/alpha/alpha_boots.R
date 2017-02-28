#### Paramètre de mortalité induite, "alpha", trade-off en fonction de beta "à la Boots12" ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant

alpha <- function(alpha1,alpha2,alpha3,alpha4,betah){
  return(alpha1+alpha2*(exp(alpha3*(betah-alpha1))+alpha4))
}

