#### Paramètre de mortalité naturelle, "b", trade-off puissance ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant

## Fonction "b" dépendant de k au format trade-off puissance

b <- function(b0, b_max, P_b, k){
  return( (b_max - b0)*k^P_b + b0 )
}

