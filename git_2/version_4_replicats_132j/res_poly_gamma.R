nbreplicats=10
var="gamma"
for(rep in 1:nbreplicats){
	setwd(paste("stepbystep",var,rep,sep=""))
	source("res_as_boots.R")
	setwd("../")
}
