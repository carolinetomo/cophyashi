require(phytools)

tree <- read.tree("0.mcc.tre")

Q<-matrix(c(-0.3,0.3,0.3,-0.3),2,2)
mtree<-sim.history(tree,Q)
plotSimmap(mtree)
rates<-c(0.25,0.75); names(rates)<-c(1,2)
x<-sim.rates(mtree,sig2=rates,plot=F)


