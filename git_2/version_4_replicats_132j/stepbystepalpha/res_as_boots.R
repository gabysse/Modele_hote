
library(deSolve)
#set.seed(123)

source("fonction_parametres.R")
source("resolution_systemes_nequ_v2.R")
source("para_simu.R")
############################
###### INITIALISATION ######
############################

nbdek=40
variationk=1/nbdek
valeurs_para=c(seq(0.1,0.9,by=0.2),seq(0.95,1.05,by=0.01),seq(1.1,1.9,by=0.2),seq(2,5,by=0.5))
decimax=100#valeur de la précision de la décimale max testée dans la variation de Pa et Px (ex : Pa->0.91 : centième : decimax=100)
total=length(valeurs_para)
#setwd("sys_equations")
#source("../creation_systeme_OK.R")
#setwd("../")
source("Reso_1302_v2.R")
comptx=seq(0,10,by=1)
longcomptx=length(comptx)




dir.create(paste("results/Pa_",Pa*decimax,"/",sep=""))
dir.create(paste("results/Pa_",Pa*decimax,"/X/",sep=""))
dir.create(paste("results/Pa_",Pa*decimax,"/Y/",sep=""))
for(i in 1:total){
	valseed="random..."#.Random.seed #On enregistre la seed en cas d'envie de reproduction parfaite des résultats.
	#pdf(paste("results/Pa_",Pa,"/num_alphaz_",x0,"_and_num_alpham_",xmax,"/compt_P_alpha",i,".pdf",sep=""), width=7, height=7) #On redirige la sorte vers le fichier pdf
	beta0=3
	beta_max=5
	P_beta=1
	  
	a0=3.577424
	a_max=4.392827
	P_a=Pa
	  
	gamma0=1
	gamma_max=1
	P_gamma=0
	  
	alpha0=comptx[x0]
	alpha_max=comptx[xmax]
	P_alpha=valeurs_para[i]
	  
	b0=1
	b_max=1
	P_b=0
	  
	q=1

	final=resolution_multiple(nbdek=nbdek,nbmutant=1000,Xini=0.1,Yini=0.1,timebfmut=10^10,dec_arrondi=2,Ymutant=10,Xmutant=25,pos_k_ini=3)
	Xf=final[[1]] ; attr(Xf,"seed")=valseed #On donne la seed en attribut aux résultats
	Yf=final[[2]] ; attr(Yf,"seed")=valseed
	#image(Xf+Yf,col=gray(rev(seq(from=0,to=1,by=0.05))),main=paste("P_alpha=",valeurs_para[i],"_alphaz_",longcomptx[x0],"_alpham_",longcomptx[xmax]))
	  
	#dev.off()#On repasse en sortie normale
	dput(Xf,file = paste("results/Pa_",Pa*decimax,"/X/compt_P_alpha",i,".data",sep=""))
	dput(Yf,file = paste("results/Pa_",Pa*decimax,"/Y/compt_P_alpha",i,".data",sep=""))
}

#source("analyse_res_diag_bif.R")
