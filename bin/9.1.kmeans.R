library(ConsensusClusterPlus)
load("../9.litingting/predict.seq.Rdata")
maxk=15
mydata=as.matrix(nlsmtr)
title="NLS"
print("-------------------NLS----------------")
pdf("../9.litingting/kmeans.nls.pdf")
nlsresults=ConsensusClusterPlus(mydata,maxK=maxk,reps = 50,pItem = 0.8,pFeature = 1,
                             title=title,clusterAlg="km",distance="euclidean",seed=1262118388.71279)
nesicl=calcICL(nlsresults,title = title)
dev.off()
mydata=as.matrix(nesmtr)
title="NES"
print("-------------------NES----------------")
pdf("../9.litingting/kmeans.nes.pdf")
nesresults=ConsensusClusterPlus(mydata,maxK=maxk,reps = 50,pItem = 0.8,pFeature = 1,
                             title=title,clusterAlg="km",distance="euclidean",seed=1262118388.71279)
nesicl=calcICL(nesresults,title = title)
dev.off()
mydata=as.matrix(unknownmtr)
title="NA"
print("-------------------unknown----------------")
pdf("../9.litingting/kmeans.unknown.pdf")
unknownresults=ConsensusClusterPlus(mydata,maxK=maxk,reps = 50,pItem = 0.8,pFeature = 1,
                             title=title,clusterAlg="km",distance="euclidean",seed=1262118388.71279)
icl=calcICL(unknownresults,title = title)
dev.off()

save(nlsresults, nesresults, unknownresults, file = "../9.litingting/kmeans.Rdata")


