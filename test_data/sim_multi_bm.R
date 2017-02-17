require(phytools)

tree <- read.tree("hap.tre")

rate<-0.05
Q<-matrix(c(-rate,rate,rate,-rate),2,2)
mtree<-sim.history(tree,Q)
plotSimmap(mtree)
rates<-c(0.05,1.5); names(rates)<-c(1,2)
x<-sim.rates(mtree,sig2=rates,plot=F)


