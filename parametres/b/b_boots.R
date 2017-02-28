#### Paramètre de mortalité naturelle, "b", trade-off en fonction de beta "à la Boots12" ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant

b <- function(b1,b2,b3,b4,betah){
  return(b1+b2*(exp(b3*(betah-b1))+b4))
}

