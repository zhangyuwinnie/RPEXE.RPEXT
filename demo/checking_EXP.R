require(RPEXE.RPEXT)
library(RPEXE.RPEXT)

# we load in the data1 dataset
data(data1)
times = data1[,1]
censor = data1[,2]
sum(censor)

# Fit the rpexe with monotonic order restriction;
pexeoutA = RPEXEv1_2(EventTime = 'EventTime',times,
                     Censor = 'Censor',censor,Trend = 'Trend',
                     trend = 3,Criticalp = 'Criticalp',criticalps=0.03)

