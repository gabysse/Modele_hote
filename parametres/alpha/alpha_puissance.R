#### Paramètre de mortalité induite, "alpha", trade-off puissance ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant


## Fonction "alpha" dépendant de k au format trade-off puissance

alpha <- function(alpha0, alpha_max, P_alpha,k){
  return( (alpha0 - alpha_max)*k^P_alpha + alpha_max )
}

