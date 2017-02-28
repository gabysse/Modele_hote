#### Paramètre de transmission, "beta", trade-off puissance ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant

## Fonction "beta" dépendant de k au format trade-off puissance

beta <- function(beta0, beta_max, P_beta, k){
  return( (beta0 - beta_max)*k^P_beta + beta_max )
}

