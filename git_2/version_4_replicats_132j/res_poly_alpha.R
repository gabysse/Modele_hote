nbreplicats=10
var="alpha"
for(rep in 1:nbreplicats){
	setwd(paste("stepbystep",var,rep,sep=""))
	source("res_as_boots.R")
	setwd("../")
}
