fichier <- "run_jobs_gamma.sh"
write("nvv=(8 20 35)",file=fichier)
write("nvvle=(10 30 50 70 90 95 96 97 98 99 100 101 102 103 104 105 110 130 150 170 190 200 250 300 350 400 450 500 )",file=fichier,append=TRUE)

write("for i in ${nvv[@]}", file=fichier, append=TRUE)
write("do", file=fichier, append=TRUE)
write("for j in ${nvvle[@]}", file=fichier, append=TRUE)
write("do", file=fichier, append=TRUE)

write("cd poskini_${i}/ajout_gamma/ajout_gamma_Pa_${j}/", file=fichier, append=TRUE)
write("qsub run_job.sh", file=fichier, append=TRUE)
write("cd ../../../", file=fichier, append=TRUE)
write("done", file=fichier, append=TRUE)
write("done", file=fichier, append=TRUE)

system(paste("chmod +x ", fichier, sep=""))
