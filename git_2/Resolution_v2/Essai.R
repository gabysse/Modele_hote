rm(list=ls())
library(deSolve)
setwd("sys_equations")
source("../creation_systeme_OK.R")
setwd("../")

#####  Initialisation  #####
kini=0.1
Xini=1
Yini=0.1
variationk=0.01
nbmutant=20
pourcentpopmutantX=0.1
pourcentpopmutantY=0.1
para_multiples=list(k=kini,X=X,Y=Y,variationk=variationk,nbmutant=nbmutant,pourcentpopmutantY=pourcentpopmutantY,pourcentpopmutantX=pourcentpopmutantX)

#On va commencer par tester uniquement des variations sur les paramètres P_X (forme du trade off)
valeurp=c(0,1,0.5,2)

for(compteura in valeurp){
  for(compteurb in valeurp){
    for(compteurbeta in valeurp){
      for(compteurgamma in valeurp){
        for(compteuralpha in valeurp){

        }
      }
    }
  }
}