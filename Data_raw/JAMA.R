library(RPEXE.RPEXT)
df<-read.csv("Data_raw/JAMABreast.CSV", sep=";",header=TRUE,na.string=c("","null","NaN","X"))
save(df, file = "data/df.rdata")