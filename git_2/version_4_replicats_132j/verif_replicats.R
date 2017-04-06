#On va chercher a rassembler les diagrammes de bifurcation par valeur de Pa. Le soucis est qu'il va falloir comparer dans un premier temps les replicats pour verifier que l'on obtiens bien tjrs la meme chose.
decimax=10#valeur de la précision de la décimale max testée dans la variation de Pa et Px (ex : Pa->0.91 : centième : decimax=100)
totPa=c(seq(0.1,1.5,by=0.1),seq(2,5,by=0.5))*decimax
lengthtotPa=length(totPa)
proptest=0.05 #Proportion de valeurs testees pour l'equilibre evolutif
ncoldata=41 #nbdek+1 (cf res_as_boots)
nbreplic=10 #Nombre de replicats
posKini=c(8,20,35) #Resistance initiale

dir.create("resultatsjpeg")

for(j in c("alpha","gamma")){
  dir.create(paste("resultatsjpeg/",j,sep=""))
  for(k in 1:lengthtotPa){#On fait varier le Pa
    jpeg(paste("resultatsjpeg/",j,"/res_Pa_",totPa[k],".jpeg",sep=""),width=1500,height=500)
    par(mfrow=c(1,3))
    for(i in posKini){  
      plot(c(0,5),0:1.05,t="n",xlab=paste("Valeur de P",j,sep=""),ylab="Valeur de résistance k finale")
      Yboots=dget(paste("boots/kini_",i,"/results/Y/compt_Pa_",k,".data",sep=""))
      Xboots=dget(paste("boots/kini_",i,"/results/X/compt_Pa_",k,".data",sep=""))#On récupère le résultat selon Boots, calculé à part.
      XYboots=Xboots+Yboots
      for(ligne in 1:ncoldata){
        if(XYboots[ligne,ncol(XYboots)]!=0){ # On prends chaque valeur de notre derniere ligne de donnees, si c'est pas un 0 on abline(a=val).
          abline(a=ligne/ncoldata,b=0)
        }
      }
      vecOK=rep(0,lengthtotPa)
      for(m in 1:lengthtotPa){# variation de Pj
        XYfinal=matrix(0,ncol=ncoldata,nrow=10) #Matrice qui comportera les résultats des =/= réplicats 
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
                points(totPa[m]/decimax,remplacement/ncoldata,col="red",pch=2)# On aura un triangle rouge vers le haut
                arrows(x0=totPa[m]/decimax,x1=totPa[m],y0=-1,y1=0,col="red",length=0.02)#... accompagné d'une flèche rouge.
              }
            }
          }else{#On est a l'equilibre evolutif 
            for(remplacement in 1:ncoldata){
              if(XYsave[remplacement]!=0){ # On prends chaque valeur de notre derniere ligne de donnees, si c'est pas un 0 on la plot.
                points(totPa[m]/decimax,remplacement/ncoldata,col="black",pch=20)
              }
            }
          }
          #On a recuperer la derniere valeur de la l-ieme iteration en sachant si c'etait ou non a l'equilibre.
        }
        #On a maintenant des matrices résumant nos 10 répétitions, dont les lignes doivent être egales.
        testvalid=0
        for(testeur in 1:(nbreplic-1)){ #Testons donc celà ! 
          test=XYfinal[testeur,]
          if(sum(XYfinal[nbreplic,]!=test)==0){#Si toutes les colonnes de XYfinal sont identiques a elle-même, ça signifie qu'on a pas de problemes de replication. On teste l'équivalence entre chaque réplicats par rapport à un autre.
            testvalid=testvalid+1
          }
        }
        if(testvalid!=nbreplic-1){ # Dans ce cas, on n'est pas à l'équilibre 
          arrows(x0=totPa[m]/decimax,x1=totPa[m]/decimax,y0=1.1,y1=1,col="green",length=0.02)
        }
      }
    }
  dev.off()
  }
}
