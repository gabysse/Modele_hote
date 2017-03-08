#### Fonction permettant de créer la liste param_k

func_param <- function(k,valeurs_param_puissance){
  # Source des types de fonctions et variables deffinissant TO 
  source("parametres/expression_parametres_puissance.R") #Contient les fonctions exprimant les TO.
  # A terme, remplacé par param_to dans la liste des arguments de la fonction
  
  betak <- beta(valeurs_param_puissance$beta0,valeurs_param_puissance$beta_max,valeurs_param_puissance$P_beta,k)
  ak <- a(valeurs_param_puissance$a0,valeurs_param_puissance$a_max,valeurs_param_puissance$P_a,k)
  alphak <- alpha(valeurs_param_puissance$alpha0,valeurs_param_puissance$alpha_max,valeurs_param_puissance$P_alpha,k)
  gammak <- gamma(valeurs_param_puissance$gamma0,valeurs_param_puissance$gamma_max,valeurs_param_puissance$P_gamma,k)
  bk <-b(valeurs_param_puissance$b0,valeurs_param_puissance$b_max,valeurs_param_puissance$P_b,k)
  qk <- cu(valeurs_param_puissance$ak,valeurs_param_puissance$bk,valeurs_param_puissance$K)
  
  return(list(betak=betak,ak=ak,alphak=alphak,gammak=gammak,bk=bk,qk=qk))
}

### Appel de param_k :
# source("fonction_parametres.R")
# param_k <- func_param(c(1.1,1.5))
# param_k

