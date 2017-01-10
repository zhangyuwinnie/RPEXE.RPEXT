library(RPEXE.RPEXT)
simple<-read.csv("Data_raw/simple.txt", sep="\t",header=FALSE,na.string=c("","null","NaN","X"))
save(simple, file = "data/simple.rdata")