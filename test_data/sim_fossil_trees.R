require("paleotree")

#set.seed(444);
record<-simFossilRecord(p=0.1, q=0.1, nruns=1,prop.bifurc=1.0,prop.cryptic=0.0,
                        totalTime=c(5,50), nTotalTaxa=c(55,55), nExtant=c(10,10))
taxa<-fossilRecord2fossilTaxa(record)
rangesCont <- sampleRanges(taxa,r=1.0)
cladogram <- taxa2cladogram(taxa,plot=F)
ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=F,timeres=T)

