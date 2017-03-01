#### Fonction permettant de créer la liste param_k
## Source des types de fonctions et variables deffinissant TO 
# source("parametres/valeurs_parametres_puissance.R") #Contient les variables definissant les TO.
# source("parametres/expression_parametres_puissance.R") #Contient les fonctions exprimant les TO.
## A terme, remplacé par param_to dans la liste des arguments de la fonction


func_param <- function(k, param_to){
 
  
  betak <- beta(param_to$beta0,param_to$beta_max,param_to$P_beta,k)
  ak <- a(param_to$a0,param_to$a_max,param_to$P_a,k)
  alphak <- alpha(param_to$alpha0,param_to$alpha_max,param_to$P_alpha,k)
  gammak <- gamma(param_to$gamma0,param_to$gamma_max,param_to$P_gamma,k)
  bk <-b(param_to$b0,param_to$b_max,param_to$P_b,k)
  qk <- cu(param_to$ak,param_to$bk,param_to$K)
  
  return(list(betak=betak,ak=ak,alphak=alphak,gammak=gammak,bk=bk,qk=qk))
}

### Appel de param_k :
# source("fonction_parametres.R")
# param_k <- func_param(c(1.1,1.5))
# param_k