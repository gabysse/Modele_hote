#### Paramètre de guérison "gamma", trade-off puissance ####

# /!\ # Paramètres initiaux non renseignés
# /!\ # Pas de valeurs correspondant à Boots 2012 pour l'instant

## Fonction "b" dépendant de k au format trade-off puissance

gamma <- function(gamma0, gamma_max, P_gamma, k){
  return( (gamma_max - gamma0)*k^P_gamma + gamma0 )
}

