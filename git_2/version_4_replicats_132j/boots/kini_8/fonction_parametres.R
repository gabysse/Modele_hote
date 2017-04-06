#### Fonction permettant de créer la liste param_k

func_param <- function(k){
  # Source des types de fonctions et variables deffinissant TO 
  #source("parametres/valeurs_parametres_puissance.R") #Contient les variables definissant les TO.
  source("parametres/expression_parametres_puissance.R") #Contient les fonctions exprimant les TO.
  # A terme, remplacé par param_to dans la liste des arguments de la fonction
  
  betak <- beta(beta0,beta_max,P_beta,k)
  ak <- a(a0,a_max,P_a,k)
  alphak <- alpha(alpha0,alpha_max,P_alpha,k)
  gammak <- gamma(gamma0,gamma_max,P_gamma,k)
  bk <-b(b0,b_max,P_b,k)
  qk <- rep(q,length(k))
  
  return(list(betak=betak,ak=ak,alphak=alphak,gammak=gammak,bk=bk,qk=qk))
}

### Appel de param_k :
# source("fonction_parametres.R")
# param_k <- func_param(c(1.1,1.5))
# param_k
