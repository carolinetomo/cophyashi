require(phytools)

tree = read.tree("./hap.tre")

rate = 0.05
Q = matrix(c(-rate,rate,rate,-rate),2,2)
mtree = sim.history(tree,Q)
rates = c(0.05,3.5); names(rates)<-c(1,2)
x = sim.rates(mtree,sig2=rates,plot=F)
x = data.frame(x)

for (i in 2:30){
    t = sim.rates(mtree,sig2=rates,plot=F)
    x[,i] = t
}

write.table(x, sep="\t",quot=F,file="hap.30.multi.csv")
