#### Paramètre de recovery, "gamma", trade-off en fonction de beta "à la Boots12" ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant


gamma <- function(gamma1,gamma2,gamma3,gamma4,betah){
  return(gamma1+gamma2*(exp(gamma3*(betah-gamma1))+gamma4))
}

