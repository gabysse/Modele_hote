comptcrea=c(seq(0.1,0.9,by=0.2),seq(0.95,1.05,by=0.01),seq(1.1,1.9,by=0.2),seq(2,5,by=0.5))
variable=c("alpha","gamma")
max=5 # Valeurs de variablemax et variable0
zero=3
posk=c(8,20,35)
replicats=10

for(posik in posk){
  dir.create(paste("poskini_",posik,sep=""))
  for(var in variable){
    dir.create(paste("poskini_",posik,"/ajout_",var,sep=""))
    for(creat in comptcrea){
	dir.create(paste("poskini_",posik,"/ajout_",var,"/ajout_",var,"_Pa_",creat*100,"/",sep=""))
	for(replic in 1:replicats){
	      file.copy(paste("stepbystep",var,sep=""),paste("poskini_",posik,"/ajout_",var,"/ajout_",var,"_Pa_",creat*100,sep=""),recursive=T)
	      file.rename(from=paste("poskini_",posik,"/ajout_",var,"/ajout_",var,"_Pa_",creat*100,"/stepbystep",var,sep=""),to=paste("poskini_",posik,"/ajout_",var,"/ajout_",var,"_Pa_",creat*100,"/stepbystep",var,replic,sep=""))
	      cat(paste("Pa=",creat,";xmax=5;x0=3;poskini=",posik,sep=""),file=paste("poskini_",posik,"/ajout_",var,"/ajout_",var,"_Pa_",creat*100,"/stepbystep",var,replic,"/para_simu.R",sep=""))
	      }		
     file.copy(paste("res_poly_",var,".R",sep=""),paste("poskini_",posik,"/ajout_",var,"/ajout_",var,"_Pa_",creat*100,sep=""))
     cat(paste("#! /bin/sh","\n","#$ -S /bin/sh","\n","#$ -j y","\n","#$ -N ",var,creat*100,"_",posik,"\n","#$ -cwd","\n","R CMD BATCH res_poly_",var,".R",sep=""),file=paste("poskini_",posik,"/ajout_",var,"/ajout_",var,"_Pa_",creat*100,"/run_job.sh",sep=""))			
    }
  }
}
