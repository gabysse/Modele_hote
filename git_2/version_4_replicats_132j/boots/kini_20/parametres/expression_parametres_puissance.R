#Paramètres exprimées en fonction des trade-off "puissance"

a <- function(a0,a_max,P_a,k){
  	return((a0 - a_max)*k^P_a + a_max )
}

alpha <- function(alpha0, alpha_max, P_alpha,k){
  	return((alpha0 - alpha_max)*k^P_alpha + alpha_max )
}

b <- function(b0, b_max, P_b, k){
  	return((b_max - b0)*k^P_b + b0 )
}

beta <- function(beta0, beta_max, P_beta, k){
	return((beta0 - beta_max)*k^P_beta + beta_max )
}

gamma <- function(gamma0, gamma_max, P_gamma, k){
  	return((gamma_max - gamma0)*k^P_gamma + gamma0 )
}
