# One-sided nominal error levels can deviate markedly 
# from the exact levels when the logrank test or the 
# generalizations of the Wilcoxon test are applied to 
# small samples. The errors are generally less serious 
# for the logrank test with continuity correction; however, 
#     in unbalanced trials when higher event rates 
#     are inferred in the smaller group, this test, 
#     too, can be highly nonconservative. On the other 
#     hand, the tests are overly conservative when lower
#     effect rates are inferred in the smaller sample. 
#     The mere utilization of a sufficiently increased 
#     continuity correction does not adequately improve 
#     the small-sample properties of the tests. These 
#     results confirm the statements of Prentice and Marek (1979) 
#     who, applying various rank tests to a highly 
#     unbalanced trial, expressed doubts concerning the 
#     validity of the significance levels for higher 
#     event rates in the smaller sample. Methods that 
#     are not based on the standard normal approximation 
#     are therefore desirable.

require(RPEXE.RPEXT)
library(RPEXE.RPEXT)

# we load in the JAMABreast dataset
data(df)
#data("JAMABreast", package = 'RPEXE.RPEXT')

# seperate variable
# Three columns in the number: indicator of the validation, censor indicator, time to drfs; 
# Two columns in the text: er status ("P" or "N"), chemo prediction: ("Rx Sensitive" or "Rx Insensitive");
vali_indi = df[,1]
bstcens = df[,2]
bsttime = df[,3]
ertxt = df[,4]
chemotxt = df[,5]
nodaltxt = df[,6]
resptxt = df[,7]
groupnum = df[,8]
t_stage = df[,9]
predresptxt = df[,10]

# Get the indices for er, nodal, group, tumor size,
# pathological response, predicted pathological response(dlda 30 predictor), 
er_p_indi = as.integer((ertxt == 'P'))
er_n_indi = as.integer((ertxt == 'N'))


nodal_p_indi = as.integer((nodaltxt == 'N1' |nodaltxt == 'N2'|nodaltxt == 'N3'))
nodal_n_indi = as.integer((nodaltxt == 'N0'))
group_p_indi= as.integer((groupnum == 3))
group_n_indi= as.integer((groupnum == 1) | (groupnum == 2))
tsize_n_indi= as.integer((t_stage == 'T1') | (t_stage == 'T2'))
tsize_p_indi= as.integer((t_stage=='T3') | (t_stage=='T4'))
res_p_indi  = as.integer((resptxt == 'pCR'))
res_n_indi  = as.integer((resptxt == 'RD'))
pres_p_indi = as.integer((predresptxt == 'pCR'))
pres_n_indi = as.integer((predresptxt == 'RD'))
chemo_p_indi = as.integer((chemotxt == 'Rx Sensitive'))
chemo_n_indi = as.integer((chemotxt == 'Rx Insensitive'))

# calculate the sample size
table(vali_indi)
table(er_p_indi[vali_indi == 0])
table(er_p_indi[vali_indi == 1])
# table(nodal_p_indi[vali_indi==0])
# table(nodal_p_indi[vali_indi==1])
# table(group_p_indi[vali_indi==0])
# table(group_p_indi[vali_indi==1])
# table(tsize_p_indi[vali_indi==0])
# table(tsize_p_indi[vali_indi==1])

res_tab = table(resptxt, predresptxt)
table(resptxt)
table(predresptxt)
# confusionMatrix(predresptxt, resptxt,positive = NULL)

rest <- vector()
for (i in 1: 508)
{
  rest[i] = paste(resptxt[i],predresptxt[i],sep = "")
}
table(rest)

# create all index
er_P_train = intersect(which(vali_indi == 0), which(er_p_indi == 1))
er_N_train = intersect(which(vali_indi == 0), which(er_n_indi == 1))
er_P_valid = intersect(which(vali_indi == 1), which(er_p_indi == 1))
er_N_valid = intersect(which(vali_indi == 1), which(er_n_indi == 1))
# nodal
er_P_train_nodal_p = intersect(er_P_train, which(nodal_p_indi == 1))
er_P_train_nodal_n = intersect(er_P_train, which(nodal_n_indi == 1))
er_N_train_nodal_p = intersect(er_N_train, which(nodal_p_indi == 1))
er_N_train_nodal_n = intersect(er_N_train, which(nodal_n_indi == 1))
er_P_valid_nodal_p = intersect(er_P_valid, which(nodal_p_indi == 1))
er_P_valid_nodal_n = intersect(er_P_valid, which(nodal_n_indi == 1))
er_N_valid_nodal_p = intersect(er_N_valid, which(nodal_p_indi == 1))
er_N_valid_nodal_n = intersect(er_N_valid, which(nodal_n_indi == 1))
length(er_P_train_nodal_p) + length(er_N_train_nodal_p)
length(er_P_train_nodal_n) + length(er_N_train_nodal_n)
length(er_P_valid_nodal_p) + length(er_N_valid_nodal_p)
length(er_P_valid_nodal_n) + length(er_N_valid_nodal_n)

# group
er_P_train_group_p = intersect(er_P_train, which(group_p_indi == 1))
er_P_train_group_n = intersect(er_P_train, which(group_n_indi == 1))
er_N_train_group_p = intersect(er_N_train, which(group_p_indi == 1))
er_N_train_group_n = intersect(er_N_train, which(group_n_indi == 1))
er_P_valid_group_p = intersect(er_P_valid, which(group_p_indi == 1))
er_P_valid_group_n = intersect(er_P_valid, which(group_n_indi == 1))
er_N_valid_group_p = intersect(er_N_valid, which(group_p_indi == 1))
er_N_valid_group_n = intersect(er_N_valid, which(group_n_indi == 1))
length(er_P_train_group_p) + length(er_N_train_group_p)
length(er_P_train_group_n) + length(er_N_train_group_n)
length(er_P_valid_group_p) + length(er_N_valid_group_p)
length(er_P_valid_group_n) + length(er_N_valid_group_n)

# tumor size
er_P_train_tsize_p = intersect(er_P_train, which(tsize_p_indi == 1))
er_P_train_tsize_n = intersect(er_P_train, which(tsize_n_indi == 1))
er_N_train_tsize_p = intersect(er_N_train, which(tsize_p_indi == 1))
er_N_train_tsize_n = intersect(er_N_train, which(tsize_n_indi == 1))
er_P_valid_tsize_p = intersect(er_P_valid, which(tsize_p_indi == 1))
er_P_valid_tsize_n = intersect(er_P_valid, which(tsize_n_indi == 1))
er_N_valid_tsize_p = intersect(er_N_valid, which(tsize_p_indi == 1))
er_N_valid_tsize_n = intersect(er_N_valid, which(tsize_n_indi == 1))
length(er_P_train_tsize_p) + length(er_N_train_tsize_p)
length(er_P_train_tsize_n) + length(er_N_train_tsize_n)
length(er_P_valid_tsize_p) + length(er_N_valid_tsize_p)
length(er_P_valid_tsize_n) + length(er_N_valid_tsize_n)

# true pathological response
er_P_train_res_p = intersect(er_P_train, which(res_p_indi == 1))
er_P_train_res_n = intersect(er_P_train, which(res_n_indi == 1))
er_N_train_res_p = intersect(er_N_train, which(res_p_indi == 1))
er_N_train_res_n = intersect(er_N_train, which(res_n_indi == 1))
er_P_valid_res_p = intersect(er_P_valid, which(res_p_indi == 1))
er_P_valid_res_n = intersect(er_P_valid, which(res_n_indi == 1))
er_N_valid_res_p = intersect(er_N_valid, which(res_p_indi == 1))
er_N_valid_res_n = intersect(er_N_valid, which(res_n_indi == 1))
length(er_P_train_res_p) + length(er_N_train_res_p)
length(er_P_train_res_n) + length(er_N_train_res_n)
length(er_P_valid_res_p) + length(er_N_valid_res_p)
length(er_P_valid_res_n) + length(er_N_valid_res_n)

# predicted pathological response from dlda 30 predictor
er_P_train_pres_p = intersect(er_P_train, which(pres_p_indi == 1))
er_P_train_pres_n = intersect(er_P_train, which(pres_n_indi == 1))
er_N_train_pres_p = intersect(er_N_train, which(pres_p_indi == 1))
er_N_train_pres_n = intersect(er_N_train, which(pres_n_indi == 1))
er_P_valid_pres_p = intersect(er_P_valid, which(pres_p_indi == 1))
er_P_valid_pres_n = intersect(er_P_valid, which(pres_n_indi == 1))
er_N_valid_pres_p = intersect(er_N_valid, which(pres_p_indi == 1))
er_N_valid_pres_n = intersect(er_N_valid, which(pres_n_indi == 1))
length(er_P_train_pres_p) + length(er_N_train_pres_p)
length(er_P_train_pres_n) + length(er_N_train_pres_n)
length(er_P_valid_pres_p) + length(er_N_valid_pres_p)
length(er_P_valid_pres_n) + length(er_N_valid_pres_n)

# Predicted good prognosis vs Poor prognosis groups
er_P_train_prechem_p = intersect(er_P_train, which(chemo_p_indi==1))
er_P_train_prechem_n = intersect(er_P_train, which(chemo_n_indi==1))

er_N_train_prechem_p = intersect(er_N_train, which(chemo_p_indi==1))
er_N_train_prechem_n = intersect(er_N_train, which(chemo_n_indi==1))

er_P_valid_prechem_p = intersect(er_P_valid, which(chemo_p_indi==1))
er_P_valid_prechem_n = intersect(er_P_valid, which(chemo_n_indi==1))

er_N_valid_prechem_p = intersect(er_N_valid, which(chemo_p_indi==1))
er_N_valid_prechem_n = intersect(er_N_valid, which(chemo_n_indi==1))


## Analysis in the paper Ex 4.3;
# Note: we used 0.02 as the critical p-value cutoff. One can 
# use different values to control the significance level.
#
# DLDA30
# ER- training, res+
# [length(er_N_train_pres_p) sum(bstcens[er_N_train_pres_p])]
#   112    36
train.er_n_pres_p = RPEXEv1_2(bsttime[er_N_train_pres_p],
                              bstcens[er_N_train_pres_p],
                              trend = 3,criticalps=0.02)

train.er_n_pres_p$times 
train.er_n_pres_p$pvalues
#     1.3771    0.1106
#     2.7324    0.0007

# ER- training, res-
train.er_n_pres_n = RPEXEv1_2(bsttime[er_N_train_pres_n],
                             bstcens[er_N_train_pres_n],
                             trend = 3,criticalps=0.02)

train.er_n_pres_n$times 
train.er_n_pres_n$pvalues

# ER- validation, res+
valid.er_n_pres_p = RPEXEv1_2(bsttime[er_N_valid_pres_p],
                              bstcens[er_N_valid_pres_p],
                              trend = 3,criticalps=0.02)

valid.er_n_pres_p$times
valid.er_n_pres_p$pvalues

# ER- validation, res-
valid.er_n_pres_n = RPEXEv1_2(bsttime[er_N_valid_pres_n],
                              bstcens[er_N_valid_pres_n],
                              trend = 3,criticalps=0.02)

valid.er_n_pres_n$times
valid.er_n_pres_n$pvalues

# ER all; res+ 
all.er_n_pres_p = RPEXEv1_2(bsttime[c(er_N_train_pres_p,er_N_valid_pres_p)],
                            bstcens[c(er_N_train_pres_p,er_N_valid_pres_p)],
                            trend = 3,criticalps=0.02)
all.er_n_pres_p$times
all.er_n_pres_p$pvalues

# ER all; res- 
all.er_n_pres_n = RPEXEv1_2(bsttime[c(er_N_train_pres_n,er_N_valid_pres_n)],
                            bstcens[c(er_N_train_pres_n,er_N_valid_pres_n)],
                            trend = 3,criticalps=0.02)
all.er_n_pres_n$times
all.er_n_pres_n$pvalues

# ACES
# Train, ER-, chemo+
train.er_n_prechem_p = RPEXEv1_2(bsttime[c(er_N_train_prechem_p)],
                                 bstcens[c(er_N_train_prechem_p)],
                                 trend = 3,criticalps=0.02)
train.er_n_prechem_p$times 
train.er_n_prechem_p$pvalues

# Train, ER-, chemo-
train.er_n_prechem_n = RPEXEv1_2(bsttime[c(er_N_train_prechem_n)],
                                bstcens[c(er_N_train_prechem_n)],
                                trend = 3,criticalps=0.02)
train.er_n_prechem_n$times 
train.er_n_prechem_n$pvalues          

# Valid, ER-, chemo+
valid.er_n_prechem_p = RPEXEv1_2(bsttime[c(er_N_valid_prechem_p)],
                                 bstcens[c(er_N_valid_prechem_p)],
                                 trend = 3,criticalps=0.02)
valid.er_n_prechem_p$times 
valid.er_n_prechem_p$pvalues  

# Valid, ER-, chemo-
valid.er_n_prechem_n = RPEXEv1_2(bsttime[c(er_N_valid_prechem_n)],
                                 bstcens[c(er_N_valid_prechem_n)],
                                 trend = 3,criticalps=0.02)
valid.er_n_prechem_n$times 
valid.er_n_prechem_n$pvalues   

# ER all; chemo+ 
all.er_n_prechem_p = RPEXEv1_2(bsttime[c(er_N_train_prechem_p,er_N_valid_prechem_p)],
                               bstcens[c(er_N_train_prechem_p,er_N_valid_prechem_p)],
                               trend = 3,criticalps=0.02)
all.er_n_prechem_p$times 
all.er_n_prechem_p$pvalues 

# ER all; chemo- 
all.er_n_prechem_n = RPEXEv1_2(bsttime[c(er_N_train_prechem_n,er_N_valid_prechem_n)],
                               bstcens[c(er_N_train_prechem_n,er_N_valid_prechem_n)],
                               trend = 3,criticalps=0.02)
all.er_n_prechem_n$times 
all.er_n_prechem_n$pvalues 

## Run the Kaplan-Meier's plot and the exponential plot;

# Example, DLDA, training, cPR
train.er_n_pres_p = RPEXEv1_2(bsttime[er_N_train_pres_p],
                              bstcens[er_N_train_pres_p],
                              trend = 3,criticalps=0.05)
# DLDA, training, RD
train.er_n_pres_n = RPEXEv1_2(bsttime[er_N_train_pres_n],
                              bstcens[er_N_train_pres_n],
                              trend = 3,criticalps=0.05)


# plot the two overlaid Kaplan-Meier curves;
time_p = train.er_n_pres_p$plotdatakme_times
# bsttime[er_N_train_pres_p]
cens_p = train.er_n_pres_p$plotdatakme_censoring
# bstcens[er_N_train_pres_p]
cens_indi = 0
time_n = train.er_n_pres_n$plotdatakme_times
cens_n = train.er_n_pres_n$plotdatakme_censoring

# figure(1)
km_blacksolid(time_n, cens_n, cens_indi)
km_redsolid(time_p, cens_p, cens_indi)
legend("topright",legend=c("DLDA30 -","DLDA30 +"),col=c("black", "red"), lty=1:1, cex=0.8)
x1 = cbind(time_p,cens_p)
x2 = cbind(time_n,cens_n)
# figure (2)
logrank(x1,x2)

#
# LOG-RANK TEST FOR KAPLAN-MEIER SURVIVAL FUNCTIONS 
# --------------------------------------------------------------------------------------------------------------
# HAZARD RATE IS AN EXPERIMENTAL FUNCTION!!!!
# Treatment 1: Hazard rate: 1.4982
# Treatment 2: Hazard rate: 1.2750
# Hazard ratio: 1.1750
# --------------------------------------------------------------------------------------------------------------
# UL			S.E.		z		p-value (2-tailed test)		alpha
# --------------------------------------------------------------------------------------------------------------
# 5.78849			3.96296		1.33448		0.18205				0.050
# --------------------------------------------------------------------------------------------------------------
# 		The survival functions are not statistically different
# 
#figure(3)
plot(train.er_n_pres_n$plotdatapexe_t100, train.er_n_pres_n$plotdatapexe_pred100, 
     xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "black",ylim = c(0,1),lwd = 2)
lines(train.er_n_pres_p$plotdatapexe_t100, train.er_n_pres_p$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "red",ylim = c(0,1),lwd = 2)
points(train.er_n_pres_p$plotdatapexe_tchange, train.er_n_pres_p$plotdatapexe_predc, 
       cex=1,pch=21, col = "red")
legend("topright",legend=c('DLDA30 -','DLDA30 +'),col=c("black", "red"), lty=2:2, cex=0.8)

##
# plot DLDA30 training;
# figure(11)
km_blacksolid(train.er_n_pres_n$plotdatakme_times, 
              train.er_n_pres_n$plotdatakme_censoring, 0)
km_redsolid(train.er_n_pres_p$plotdatakme_times, 
            train.er_n_pres_p$plotdatakme_censoring, 0)
lines(train.er_n_pres_n$plotdatapexe_t100, train.er_n_pres_n$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "black",ylim = c(0,1),lwd = 2)
lines(train.er_n_pres_p$plotdatapexe_t100, train.er_n_pres_p$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "red",ylim = c(0,1),lwd = 2)
points(train.er_n_pres_p$plotdatapexe_tchange, train.er_n_pres_p$plotdatapexe_predc, 
       cex=1,pch=24, col = "red",lwd=4)
legend("bottomleft",legend=c('KME, DLDA30 RD, training','KME, DLDA30 pCR, training', 
                             'RPEXE, DLDA30 RD, training','RPEXE, DLDA30 pCR, training',
                             'Significant failure rate change'),col=c("black", "red","black", "red","red"), 
       lty=c(1,1,2,2,NA), pch = c(NA,NA,NA,NA,24), cex=0.8)

# figure(12)
# plot DLDA30 testing;
km_blacksolid(valid.er_n_pres_n$plotdatakme_times, 
              valid.er_n_pres_n$plotdatakme_censoring, 0)
km_redsolid(valid.er_n_pres_p$plotdatakme_times, 
            valid.er_n_pres_p$plotdatakme_censoring, 0)
lines(valid.er_n_pres_n$plotdatapexe_t100, valid.er_n_pres_n$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "black",ylim = c(0,1),lwd = 2)
lines(valid.er_n_pres_p$plotdatapexe_t100, valid.er_n_pres_p$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "red",ylim = c(0,1),lwd = 2)
points(valid.er_n_pres_p$plotdatapexe_tchange, valid.er_n_pres_p$plotdatapexe_predc, 
       cex=1,pch=24, col = "red", lwd = 4)
legend("bottomleft",legend=c('KME, DLDA30 RD, validation','KME, DLDA30 pCR, validation', 
                             'RPEXE, DLDA30 RD, validation','RPEXE, DLDA30 pCR, validation',
                             'Significant failure rate change'),col=c("black", "red","black", "red","red"), 
       lty=c(1,1,2,2,NA), pch = c(NA,NA,NA,NA,24), cex=0.8)

# figure(13);
# plot ACES training;
km_blacksolid(train.er_n_prechem_n$plotdatakme_times, 
              train.er_n_prechem_n$plotdatakme_censoring, 0)
km_redsolid(train.er_n_prechem_p$plotdatakme_times, 
            train.er_n_prechem_p$plotdatakme_censoring, 0)
lines(train.er_n_prechem_n$plotdatapexe_t100, train.er_n_prechem_n$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "black",ylim = c(0,1),lwd = 2)
lines(train.er_n_prechem_p$plotdatapexe_t100, train.er_n_prechem_p$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "red",ylim = c(0,1),lwd = 2)
points(train.er_n_prechem_p$plotdatapexe_tchange, train.er_n_prechem_p$plotdatapexe_predc, 
       cex=1,pch=24, col = "red", lwd = 4)
legend("bottomleft",legend=c('KME, ACES Insensitive, training','KME, ACES Sensitive, training', 
                             'RPEXE, ACES Insensitive, training','RPEXE, ACES Sensitive, training',
                             'Significant failure rate change'),col=c("black", "red","black", "red","red"), 
       lty=c(1,1,2,2,NA), pch = c(NA,NA,NA,NA,24), cex=0.8)

# figure(14);
# plot ACES testing;
km_blacksolid(valid.er_n_prechem_n$plotdatakme_times, 
              valid.er_n_prechem_n$plotdatakme_censoring, 0)
km_redsolid(valid.er_n_prechem_p$plotdatakme_times, 
            valid.er_n_prechem_p$plotdatakme_censoring, 0)
lines(valid.er_n_prechem_n$plotdatapexe_t100, valid.er_n_prechem_n$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "black",ylim = c(0,1),lwd = 2)
lines(valid.er_n_prechem_p$plotdatapexe_t100, valid.er_n_prechem_p$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "red",ylim = c(0,1),lwd = 2)
points(valid.er_n_prechem_p$plotdatapexe_tchange, valid.er_n_prechem_p$plotdatapexe_predc, 
       cex=1,pch=24, col = "red", lwd = 4)
legend("bottomleft",legend=c('KME, ACES Insensitive, validation','KME, ACES Sensitive, validation', 
                             'RPEXE, ACES Insensitive, validation','RPEXE, ACES Sensitive, validation',
                             'Significant failure rate change'),col=c("black", "red","black", "red","red"), 
       lty=c(1,1,2,2,NA), pch = c(NA,NA,NA,NA,24), cex=0.8)

# figure(15);
# plot DLDA30 combined;
km_blacksolid(all.er_n_pres_n$plotdatakme_times, 
              all.er_n_pres_n$plotdatakme_censoring, 0)
km_redsolid(all.er_n_pres_p$plotdatakme_times, 
            all.er_n_pres_p$plotdatakme_censoring, 0)
lines(all.er_n_pres_n$plotdatapexe_t100, all.er_n_pres_n$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "black",ylim = c(0,1),lwd = 2)
lines(all.er_n_pres_p$plotdatapexe_t100, all.er_n_pres_p$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "red",ylim = c(0,1),lwd = 2)
points(all.er_n_pres_p$plotdatapexe_tchange, all.er_n_pres_p$plotdatapexe_predc, 
       cex=1,pch=24, col = "red", lwd = 4)
legend("bottomleft",legend=c('KME, DLDA30 RD, combined','KME, DLDA30 pCR, combined', 
                             'RPEXE, DLDA30 RD, combined','RPEXE, DLDA30 pCR, combined',
                             'Significant failure rate change'),col=c("black", "red","black", "red","red"), 
       lty=c(1,1,2,2,NA), pch = c(NA,NA,NA,NA,24), cex=0.8)

# figure(16);
# plot ACES combined;
km_blacksolid(all.er_n_prechem_n$plotdatakme_times, 
              all.er_n_prechem_n$plotdatakme_censoring, 0)
km_redsolid(all.er_n_prechem_p$plotdatakme_times, 
            all.er_n_prechem_p$plotdatakme_censoring, 0)
lines(all.er_n_prechem_n$plotdatapexe_t100, all.er_n_prechem_n$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "black",ylim = c(0,1),lwd = 2)
lines(all.er_n_prechem_p$plotdatapexe_t100, all.er_n_prechem_p$plotdatapexe_pred100, 
      xlab="Time", ylab="Estimated survival functions", type = "s",lty = 2, col = "red",ylim = c(0,1),lwd = 2)
points(all.er_n_prechem_p$plotdatapexe_tchange, all.er_n_prechem_p$plotdatapexe_predc, 
       cex=1,pch=24, col = "red", lwd = 4)
legend("bottomleft",legend=c('KME, ACES Insensitive, combined','KME, ACES Sensitive, combined', 
                             'RPEXE, ACES Insensitive, combined','RPEXE, ACES Sensitive, combined',
                             'Significant failure rate change'),col=c("black", "red","black", "red","red"), 
       lty=c(1,1,2,2,NA), pch = c(NA,NA,NA,NA,24), cex=0.8)

################################
## Logrank tests illustration;
timesn = train.er_n_pres_n$plotdatakme_times
censn  = 1-train.er_n_pres_n$plotdatakme_censoring
timesp = train.er_n_pres_p$plotdatakme_times
censp  = 1-train.er_n_pres_p$plotdatakme_censoring
x1 = cbind(timesn, censn)
x2 = cbind(timesp, censp)
logrank_Yates228(x1,x2)
# prepare data for sas;
# data1 = [timesn 1-censn ones(length(timesn),1);timesp 1-censp 2*ones(length(timesp),1)]

# Run the exponential test to compare the failure rate before 
train.er_n_pres_p$times
censn = 1-censn
censp = 1-censp
indn = which(timesn>=2.7325)
indp = which(timesp>=2.7325)
censn[indn] = 0
censp[indp] = 0
timesn[indn] = 2.7325
timesp[indp] = 2.7325
ttot1 = sum(timesn) 
d1 = sum(censn)
ttot2 = sum(timesp) 
d2 = sum(censp)
reout=exact_pvalue(ttot1,ttot2,d1,d2,0)
a2 = reout[1]
p2 = reout[2]
# two sided exact test p= 0.4949;
reout=exact_pvalue(ttot1,ttot2,d1,d2,2)
a2 = reout[1]
p2 = reout[2]
# one-sided test exact p-value 0.4647;

timesn = valid.er_n_pres_n$plotdatakme_times
censn  = 1-valid.er_n_pres_n$plotdatakme_censoring
timesp = valid.er_n_pres_p$plotdatakme_times
censp  = 1-valid.er_n_pres_p$plotdatakme_censoring
x1 = cbind(timesn, censn)
x2 = cbind(timesp, censp)
logrank_Yates228(x1,x2)
# Run the exponential test to compare the failure rate before 
valid.er_n_pres_p$times_c
censn = 1-censn
censp = 1-censp
indn = which(timesn>=2.1794)
indp = which(timesp>=2.1794)
censn[indn] = 0
censp[indp] = 0
timesn[indn] = 2.1794
timesp[indp] = 2.1794
ttot1 = sum(timesn) 
d1 = sum(censn)
ttot2 = sum(timesp) 
d2 = sum(censp)
reout=exact_pvalue(ttot1,ttot2,d1,d2,0)
a2 = reout[1]
p2 = reout[2]
# two sided exact test p= 0.4949;
reout=exact_pvalue(ttot1,ttot2,d1,d2,2)
a2 = reout[1]
p2 = reout[2]
# one-sided test exact p-value 0.4647;

timesn = train.er_n_prechem_n$plotdatakme_times
censn  = 1-train.er_n_prechem_n$plotdatakme_censoring
timesp = train.er_n_prechem_p$plotdatakme_times
censp  = 1-train.er_n_prechem_p$plotdatakme_censoring
x1 = cbind(timesn, censn)
x2 = cbind(timesp, censp)
logrank_Yates228(x1,x2)
# Run the exponential test to compare the failure rate before 
train.er_n_prechem_n$times_c
censn = 1-censn
censp = 1-censp
indn = which(timesn>=2.7325)
indp = which(timesp>=2.7325)
censn[indn] = 0
censp[indp] = 0
timesn[indn] = 2.7325
timesp[indp] = 2.7325
ttot1 = sum(timesn) 
d1 = sum(censn)
ttot2 = sum(timesp) 
d2 = sum(censp)
reout=exact_pvalue(ttot1,ttot2,d1,d2,0)
a2 = reout[1]
p2 = reout[2]
# two sided exact test p= 0.4949;
reout=exact_pvalue(ttot1,ttot2,d1,d2,2)
a2 = reout[1]
p2 = reout[2]
# one-sided test exact p-value 0.4647;

timesn = valid.er_n_prechem_n$plotdatakme_times
censn  = 1-valid.er_n_prechem_n$plotdatakme_censoring
timesp = valid.er_n_prechem_p$plotdatakme_times
censp  = 1-valid.er_n_prechem_p$plotdatakme_censoring
x1 = cbind(timesn, censn)
x2 = cbind(timesp, censp)
logrank(x1,x2)
logrank_Yates228(x1,x2)

# Run the exponential test to compare the failure rate before 
# the change point time 2.17933;
censn = 1-censn
censp = 1-censp
indn = which(timesn>=2.17933)
indp = which(timesp>=2.17933)
censn[indn] = 0
censp[indp] = 0
timesn[indn] = 2.17933
timesp[indp] = 2.17933
ttot1 = sum(timesn)
d1 = sum(censn)
ttot2 = sum(timesp) 
d2 = sum(censp)
reout=exact_pvalue(ttot1,ttot2,d1,d2,0)
a2 = reout[1]
p2 = reout[2]
# two sided exact test p= 0.4949;
reout=exact_pvalue(ttot1,ttot2,d1,d2,2)
a2 = reout[1]
p2 = reout[2]
# one-sided test exact p-value 0.0167

# Note: Thelog rank test at cut off 2.17933 gives p-value 0.0348;

ttot1
ttot2
d1
d2
ttot1/d1
ttot2/d2

