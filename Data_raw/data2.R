library(RPEXE.RPEXT)
data2<-read.csv("Data_raw/data2.csv", sep=";",header=FALSE,na.string=c("","null","NaN","X"))
save(data2, file = "data/data2.rdata")