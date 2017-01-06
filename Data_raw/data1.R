library(RPEXE.RPEXT)
data1<-read.csv("Data_raw/data1.csv", sep=";",header=FALSE,na.string=c("","null","NaN","X"))
save(data1, file = "data/data1.rdata")