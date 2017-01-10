# Example of the RPEXE in R
# None Small Cell Lung cancer example

require(RPEXE.RPEXT)
library(RPEXE.RPEXT)

# we load in the data2 dataset
data(simple)
times = simple[,2]

cens  = simple[,1]

test2.decrease = RPEXEv1_2(times,cens,trend = 1,criticalps=0.05)
test2.umbrella = RPEXEv1_2(times,cens,trend = 4,criticalps=0.05)

# Results:  
test2.decrease$times
# $times
# [1]  1.8931507 23.4602740 47.8547945 12.7095890 24.7095890  0.0849315  2.3863014 18.5945205 28.0301370

test2.decrease$pvalues
# $pvalues
# [1] 9.779920e-01 9.607891e-01 7.538870e-01 7.045244e-01 6.876457e-01 6.012978e-01 3.059765e-01 8.523531e-02
# [9] 9.100265e-07

test2.decrease$trend
# $trend
# [1] "Decreasing falilure rate"

test2.decrease$changet
# $changet
# NULL



test2.umbrella$times
# $times
# [1] 23.4602740 47.8547945  5.2794521 24.7095890 12.7095890  5.3452055  0.6438356 18.5945205  5.3780822
# [10]  5.4109589 28.0301370

test2.umbrella$pvalues
# $pvalues
# [1] 9.607891e-01 7.538870e-01 7.207672e-01 6.876457e-01 6.499540e-01 4.152210e-01 1.854924e-01 1.441013e-01
# [9] 4.383775e-03 3.085773e-01 9.100265e-07

test2.umbrella$trend
# $trend
# [1] "Increasing-decreasing failure rate"

test2.umbrella$changet
# $changet
# [1] 5.378082
