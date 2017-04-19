#On va chercher a rassembler les diagrammes de bifurcation par valeur de Pa. Le soucis est qu'il va falloir comparer dans un premier temps les replicats pour verifier que l'on obtiens bien tjrs la meme chose.
decimax=100#valeur de la precision de la decimale max testee dans la variation de Pa et Px (ex : Pa->0.91 : centieme : decimax=100)
totPa=c(seq(0.1,0.9,by=0.2),seq(0.95,1.05,by=0.01),seq(1.1,1.9,by=0.2),seq(2,5,by=0.5))*decimax
lengthtotPa=length(totPa)
proptest=0.05 #Proportion de valeurs testees pour l'equilibre evolutif
ncoldata=41 #nbdek+1 (cf res_as_boots)
nbreplic=10 #Nombre de replicats
posKini=c(8,20,35) #Resistance initiale
trace_TO=seq(0,1,by=0.01)
comp_a_gnon=seq(3,5,by=0.01)

vecdecoul=rainbow(lengthtotPa)


par(cex.lab=1.5)
dir.create("resultatsjpeg")

for(j in c("alpha","gamma")){
  if(j=="alpha"){
    fon1=function(P_alpha,sequence){return((3-5)*sequence^P_alpha+5)}
    fon2=function(P_a,sequence){return((3-5)*sequence^P_a+5)}
    fon3=function(P_alpha,sequence){return((3-5)*((sequence-5)/(3-5))^(1/P_alpha)+5)}
    fon4=function(P_alpha,sequence,P_a){return((3-5)*((sequence-5)/(3-5))^(P_a/P_alpha)+5)}
    fon5=function(P_a,sequence){return((3-5)*((sequence-5)/(3-5))^(P_a)+5)}
  }
  if(j=="gamma"){
    fon1=function(P_gamma,sequence){return((5-3)*sequence^P_gamma+3)}
    fon2=function(P_a,sequence){return((3-5)*sequence^P_a+5)}
    fon3=function(P_gamma,sequence){return((3-5)*((sequence-3)/(5-3))^(1/P_gamma)+5)}
    fon4=function(P_gamma,sequence,P_a){return((3-5)*((sequence-3)/(5-3))^(P_a/P_gamma)+5)}
    fon5=function(P_a,sequence){return((3-5)*((sequence-5)/(3-5))^(P_a)+5)}
  }
  
  dir.create(paste("resultatsjpeg/",j,sep=""))
  for(k in 1:lengthtotPa){#On fait varier le Pa
    jpeg(paste("resultatsjpeg/",j,"/res_Pa_",totPa[k],".jpeg",sep=""),width=1500,height=750)
    par(mfrow=c(2,4))
    for(i in posKini){  
      plot(seq(0.1,5.1,by=1),seq(0,1,by=0.2),t="n",xlog="T",log="x",xlab=paste("Valeur de P",j,sep=""),ylab="Valeur de resistance k finale",main=paste("Resistance initiale = ",i,sep=""),cex.lab=1.5,lwd=1)
      Yboots=dget(paste("boots/kini_",i,"/results/Y/compt_Pa_",k,".data",sep=""))
      Xboots=dget(paste("boots/kini_",i,"/results/X/compt_Pa_",k,".data",sep=""))#On recupere le resultat selon Boots, calcule a part.
      XYboots=Xboots+Yboots
      for(ligne in 1:ncoldata){
        if(XYboots[ligne,ncol(XYboots)]!=0){ # On prends chaque valeur de notre derniere ligne de donnees, si c'est pas un 0 on abline(a=val).
          abline(a=ligne/ncoldata,b=0,lwd=1.2)
        }
      }
      vecOK=rep(0,lengthtotPa)
      for(m in 1:lengthtotPa){# variation de Pj
        XYfinal=matrix(0,ncol=ncoldata,nrow=10) #Matrice qui comportera les resultats des =/= replicats 
        for(l in 1:10){
          #Pour chaque valeur de Pa (forme de TO entre beta et a), on va chercher les replicats et verifier que tous sont egaux dans les dernieres lignes (equilibre evolutif independant des etapes aleatoires).
          X=dget(paste("poskini_",i,"/ajout_",j,"/ajout_",j,"_Pa_",totPa[k],"/stepbystep",j,l,"/results/Pa_",totPa[k],"/X/compt_P_",j,m,".data",sep=""))
          Y=dget(paste("poskini_",i,"/ajout_",j,"/ajout_",j,"_Pa_",totPa[k],"/stepbystep",j,l,"/results/Pa_",totPa[k],"/Y/compt_P_",j,m,".data",sep=""))
          XY=X+Y
          XYsave=XY[,ncol(XY)] #On recupere la derniere ligne des donnees (resultats apres la derniere mutation).
          XYfinal[l,]=XY[,ncol(XY)]

          #On teste l'equilibre
          testeq=XY[,(ncol(XY)-round(ncol(XY)*proptest)):(ncol(XY)-1)]
          if(sum(XY[,ncol(XY)]!=testeq)>0){# On n'atteint pas l'equilibre evolutif
            for(remplacement in 1:ncoldata){
              if(XYsave[remplacement]!=0){ # On prends chaque valeur de notre derniere ligne de donnees, si c'est pas un 0 on la plot.
                points(totPa[m]/decimax,remplacement/ncoldata,col=vecdecoul[m],pch=20,cex.lab=1.5,lwd=1)# On aura un triangle rouge vers le haut
                arrows(x0=totPa[m]/decimax,x1=totPa[m]/decimax,y0=-1,y1=0.03,col="red",length=0.025,lwd=5)#... accompagne d'une fleche rouge.
              }
            }
          }else{#On est a l'equilibre evolutif
            compteurlocal=NULL 
            for(remplacement in 1:ncoldata){
              if(XYsave[remplacement]!=0){ # On prends chaque valeur de notre derniere ligne de donnees, si c'est pas un 0 on la plot.
                compteurlocal=c(compteurlocal,remplacement)
              }
            }
            #if(length(compteurlocal==1)){#Ne marche pas pour une raison indeterminee actuellement
              points(rep(totPa[m]/decimax,times=length(compteurlocal)),compteurlocal/ncoldata,col=vecdecoul[m],pch=20,cex.lab=1.5,lwd=1)
            #}
            #if(length(compteurlocal>1)){
            #  points(rep(totPa[m]/decimax,times=length(compteurlocal)),compteurlocal/ncoldata,col="black",pch="|")
            #}
          }
          #On a recupere la derniere valeur de la l-ieme iteration en sachant si c'etait ou non a l'equilibre.
        }
        #On a maintenant des matrices resumant nos 10 repetitions, dont les lignes doivent etre egales.
        testvalid=0
        for(testeur in 1:(nbreplic-1)){ #Testons donc cela ! 
          test=XYfinal[testeur,]
          if(sum(XYfinal[nbreplic,]!=test)==0){#Si toutes les colonnes de XYfinal sont identiques a elle-meme, ça signifie qu'on a pas de problemes de replication. On teste l'equivalence entre chaque replicats par rapport a un autre.
            testvalid=testvalid+1
          }
        }
        if(testvalid!=nbreplic-1){ # Dans ce cas, on n'est pas a l'equilibre 
          arrows(x0=totPa[m]/decimax,x1=totPa[m]/decimax,y0=1.1,y1=1,col="darkgreen",length=0.025,lwd=5)
        }
      }
    }

    plot(c(2.95,5.05),c(2.95,5.05),t="n",xlab="beta",ylab="a",cex.lab=1.5) #On regarde le trade-off a en fonction de beta.
    
    TO=NULL
      for(poulp in comp_a_gnon){
        TO=c(TO,fon5(totPa[k]/decimax,poulp))
      }
    points(comp_a_gnon,TO,t="l",col="black",lwd=1)
    

    plot(c(0,1),c(2.95,5.05),t="n",xlab="k",ylab=paste("P_",j,sep=""),cex.lab=1.5) #On va tracer les relations P_x en fonction de k
    for(trace in 1:lengthtotPa){
      TO=NULL
      for(poulpy in trace_TO){
        TO=c(TO,fon1(totPa[trace]/decimax,poulpy))
      }
      points(trace_TO,TO,t="l",col=vecdecoul[trace],lwd=1)
    }
 
    plot(c(0,1),c(2.95,5.05),t="n",xlab="k",ylab="P_a",cex.lab=1.5) #On trace la relation P_a en fonction de k
    TO=NULL
    for(poulpon in trace_TO){
      TO=c(TO,fon2(totPa[k]/decimax,poulpon))
    }
    points(trace_TO,TO,t="l",col="black",lwd=1)


    plot(c(2.95,5.05),c(2.95,5.05),t="n",xlab=j,ylab="beta",cex.lab=1.5) #On trace le Trade-Off entre x et a.
    for(trace2 in 1:lengthtotPa){
    TO=NULL
      for(poulp in comp_a_gnon){
        TO=c(TO,fon3(totPa[trace2]/decimax,poulp))
      }
    points(comp_a_gnon,TO,t="l",col=vecdecoul[trace2],lwd=1)
    }

    plot(c(2.95,5.05),c(2.95,5.05),t="n",xlab=j,ylab="a",cex.lab=1.5) #On trace le Trade-Off entre x et a.
    for(trace3 in 1:lengthtotPa){
    TO=NULL
      for(poulp in comp_a_gnon){
        TO=c(TO,fon4(totPa[trace3]/decimax,poulp,totPa[k]/decimax))
      }
    points(comp_a_gnon,TO,t="l",col=vecdecoul[trace3],lwd=1)
    }

  dev.off()
  }
}
