require("paleotree")

#set.seed(444);
record<-simFossilRecord(p=0.1, q=0.1, nruns=1,prop.bifurc=1.0,prop.cryptic=0.0,
                        totalTime=10, nTotalTaxa=20, nExtant=c(5,10))
taxa<-fossilRecord2fossilTaxa(record)
rangesCont <- sampleRanges(taxa,r=1.0)
cladogram <- taxa2cladogram(taxa,plot=F)
ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=F,timeres=T)
rangeSamp <- sampleRanges(taxa,r = 0.5,ranges.only=F)

