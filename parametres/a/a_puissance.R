#### Paramètre de fécondité, "a", trade-off puissance ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant


## Fonction "a" dépendant de k au format trade-off puissance

a <- function(a0,a_max,P_a,k){
  return( (a0 - a_max)*k^P_a + a_max )
}

