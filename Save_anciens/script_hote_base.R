## pas besoin de faire toute une matrice, on peut direct faire calculer les valeurs sur le tas



# library n?cessaires :
library(deSolve)
# Debnote : il faut update R pour avoir deSolve ##


####?Premi?re Etape : population initiale ####

### Param?tres initiaux ? sp?cifier
h <- 1 
Xi <- 1
Yi <- 0.1
pop = c(Xi, Yi)

### Automatisation pour les autres param?tres

betah <- 5 - h/5

ah <- 4 - 0.84*(exp(-4.89*(betah - 4)) -1)

b <- 1

q <- 1

alpha <- 1

gamma <- 1

H <- sum(pop)

L <- alpha + betah +gamma


### premi?re r?solution num?rique :

duree <- 5 # temps que dure la simulation
ta <- seq(1, 30, by=1/50) # d?coupage des pas de temps


param_h <- c(ah, b, H, q, alpha, betah, gamma, pop, L, t )

hote_h <- function(t, pop, param_h){
  Xh <- pop[1]
  Yh <- pop[2]
  H <- sum(pop)
  with(as.list(param_h), {
    
    dXh <- ah*Xh - q*H*Xh - b*Xh - sum(betah*Xh*Yh) + gamma*Yh
    dYh <- sum(betah*Xh*Yh) - L*Yh
    dH <- dXh + dYh
    
    return(list(c(dXh, dYh, dH)))
  })
}

out <- ode(y=c(Xi, Yi, H),times=c(seq(0,duree,0.01)) ,func=hote_h, parms=param_h)

out1 <- ode(y=c(Xi, Yi, H),times=seq(0,duree, 0.5) ,func=hote_h, parms=param_h)

out2 <- ode(y=c(Xi, Yi, H),times=c(seq(0,duree, 0.2)) ,func=hote_h, parms=param_h)

out3 <- ode(y=c(Xi, Yi, H),times=c(seq(0,duree, 1)) ,func=hote_h, parms=param_h)

### Repr?sentation graphique :

par(mfrow=c(2, 2))


plot(out[, 1], out[, 2], type="l")
plot(out[, 1], out[, 3], type="l")
plot(out[, 2], out[, 3], type="l")

diagnostics(out)


plot(out[, 2], out[, 3], type="l", main="step = 0.01")
plot(out2[, 2], out2[, 3], type="l", main="step = 0.2")
plot(out1[, 2], out1[, 3], type="l", main="step = 0.5")
plot(out3[, 2], out3[, 3], type="l", main="step = 1")



diagnostics(out)
diagnostics(out1)
diagnostics(out2)
diagnostics(out3)






