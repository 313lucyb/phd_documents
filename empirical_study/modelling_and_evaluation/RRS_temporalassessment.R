library(foreign)
library(rms)
library(ROCR)
library(pROC)
library(pROC)
library(plyr)
library(colorRamps)
library(grDevices)
library(DescTools)
library(dummies)
library(RColorBrewer)
cols <-brewer.pal(11, "RdBu")

# Setting up temporal assessment datasets

test0<-read.dta("baseline_binary.dta")
test0$smoke<-factor(test0$smoke)
test0$prev_sinf<-0
test0$prev_tnf<-1
test0$age2<-0
test0$treatfail<-0
test0$steroids2<-0
test0$lowsteroids<-0
for (i in (1:length(test0$studyno))){
  if (test0$age[i] > 60) {test0$age2[i] = 1}
  if (test0$previous_dmards[i] > 5) {test0$treatfail[i] = 1}
  if ((test0$steroids[i] == 1) & 
      (test0$diabetes[i] == 0) & (test0$smoke[i] == 0)){
    test0$steroids2[i] = 1}
  if ((test0$steroids[i] == 1) & (test0$steroids2[i] == 0)){
    test0$lowsteroids[i] = 1}
}

test1<-read.dta("validation_6binary.dta")
test1$smoke<-factor(test1$smoke)
test1$prev_tnf<-1
test1$age2<-0
test1$treatfail<-0
test1$steroids2<-0
test1$lowsteroids<-0
for (i in (1:length(test1$studyno))){
  if (test1$age[i] > 60) {test1$age2[i] = 1}
  if (test1$previous_dmards[i] > 5) {test1$treatfail[i] = 1}
  if ((test1$steroids[i] == 1) & 
      (test1$diabetes[i] == 0) & (test1$smoke[i] == 0)){
    test1$steroids2[i] = 1}
  if ((test1$steroids[i] == 1) & (test1$steroids2[i] == 0)){
    test1$lowsteroids[i] = 1}
}

test2<-read.dta("validation_12binary.dta")
test2$smoke<-factor(test2$smoke)
test2$prev_tnf<-1
test2$age2<-0
test2$treatfail<-0
test2$steroids2<-0
test2$lowsteroids<-0
for (i in (1:length(test2$studyno))){
  if (test2$age[i] > 60) {test2$age2[i] = 1}
  if (test2$previous_dmards[i] > 5) {test2$treatfail[i] = 1}
  if ((test2$steroids[i] == 1) & 
      (test2$diabetes[i] == 0) & (test2$smoke[i] == 0)){
    test2$steroids2[i] = 1}
  if ((test2$steroids[i] == 1) & (test2$steroids2[i] == 0)){
    test2$lowsteroids[i] = 1}
}

test3<-read.dta("validation_18binary.dta")
test3$smoke<-factor(test3$smoke)
test3$prev_tnf<-1
test3$age2<-0
test3$treatfail<-0
test3$steroids2<-0
test3$lowsteroids<-0
for (i in (1:length(test3$studyno))){
  if (test3$age[i] > 60) {test3$age2[i] = 1}
  if (test3$previous_dmards[i] > 5) {test3$treatfail[i] = 1}
  if ((test3$steroids[i] == 1) & 
      (test3$diabetes[i] == 0) & (test3$smoke[i] == 0)){
    test3$steroids2[i] = 1}
  if ((test3$steroids[i] == 1) & (test3$steroids2[i] == 0)){
    test3$lowsteroids[i] = 1}
}

test4<-read.dta("validation_24binary.dta")
test4$smoke<-factor(test4$smoke)
test4$prev_tnf<-1
test4$age2<-0
test4$treatfail<-0
test4$steroids2<-0
test4$lowsteroids<-0
for (i in (1:length(test4$studyno))){
  if (test4$age[i] > 60) {test4$age2[i] = 1}
  if (test4$previous_dmards[i] > 5) {test4$treatfail[i] = 1}
  if ((test4$steroids[i] == 1) & 
      (test4$diabetes[i] == 0) & (test4$smoke[i] == 0)){
    test4$steroids2[i] = 1}
  if ((test4$steroids[i] == 1) & (test4$steroids2[i] == 0)){
    test4$lowsteroids[i] = 1}
}

test5<-read.dta("validation_30binary.dta")
test5$smoke<-factor(test5$smoke)
test5$prev_tnf<-1
test5$age2<-0
test5$treatfail<-0
test5$steroids2<-0
test5$lowsteroids<-0
for (i in (1:length(test5$studyno))){
  if (test5$age[i] > 60) {test5$age2[i] = 1}
  if (test5$previous_dmards[i] > 5) {test5$treatfail[i] = 1}
  if ((test5$steroids[i] == 1) & 
      (test5$diabetes[i] == 0) & (test5$smoke[i] == 0)){
    test5$steroids2[i] = 1}
  if ((test5$steroids[i] == 1) & (test5$steroids2[i] == 0)){
    test5$lowsteroids[i] = 1}
}

test6<-read.dta("validation_36binary.dta")
test6$smoke<-factor(test6$smoke)
test6$prev_tnf<-1
test6$age2<-0
test6$treatfail<-0
test6$steroids2<-0
test6$lowsteroids<-0
for (i in (1:length(test6$studyno))){
  if (test6$age[i] > 60) {test6$age2[i] = 1}
  if (test6$previous_dmards[i] > 5) {test6$treatfail[i] = 1}
  if ((test6$steroids[i] == 1) & 
      (test6$diabetes[i] == 0) & (test6$smoke[i] == 0)){
    test6$steroids2[i] = 1}
  if ((test6$steroids[i] == 1) & (test6$steroids2[i] == 0)){
    test6$lowsteroids[i] = 1 
  }
}


# temporal assessment - predictions (1) [all low dose]

sum10<- (-4.191) + 0.470*test0$age2 + 0.309*(test0$ovmean - 3.16) +
  0.484*test0$lung + 0.415*test0$renal + 0.992*test0$prev_sinf +
  0.397*test0$treatfail + 0.782*test0$steroids + 0.589*test0$prev_tnf
test0$pred10<-1-exp(-exp(sum10))

sum11<- (-4.191) + 0.470*test1$age2 + 0.309*(test1$ovmean - 3.16) +
  0.484*test1$lung + 0.415*test1$renal + 0.992*test1$prev_sinf +
  0.397*test1$treatfail + 0.782*test1$steroids + 0.589*test1$prev_tnf
test1$pred11<-1-exp(-exp(sum11))

sum12<- (-4.191) + 0.470*test2$age2 + 0.309*(test2$ovmean - 3.16) +
  0.484*test2$lung + 0.415*test2$renal + 0.992*test2$prev_sinf +
  0.397*test2$treatfail + 0.782*test2$steroids + 0.589*test2$prev_tnf
test2$pred12<-1-exp(-exp(sum12))

sum13<- (-4.191) + 0.470*test3$age2 + 0.309*(test3$ovmean - 3.16) +
  0.484*test3$lung + 0.415*test3$renal + 0.992*test3$prev_sinf +
  0.397*test3$treatfail + 0.782*test3$steroids + 0.589*test3$prev_tnf
test3$pred13<-1-exp(-exp(sum13))

sum14<- (-4.191) + 0.470*test4$age2 + 0.309*(test4$ovmean - 3.16) +
  0.484*test4$lung + 0.415*test4$renal + 0.992*test4$prev_sinf +
  0.397*test4$treatfail + 0.782*test4$steroids + 0.589*test4$prev_tnf
test4$pred14<-1-exp(-exp(sum14))

sum15<- (-4.191) + 0.470*test5$age2 + 0.309*(test5$ovmean - 3.16) +
  0.484*test5$lung + 0.415*test5$renal + 0.992*test5$prev_sinf +
  0.397*test5$treatfail + 0.782*test5$steroids + 0.589*test5$prev_tnf
test5$pred15<-1-exp(-exp(sum15))

sum16<- (-4.191) + 0.470*test6$age2 + 0.309*(test6$ovmean - 3.16) +
  0.484*test6$lung + 0.415*test6$renal + 0.992*test6$prev_sinf +
  0.397*test6$treatfail + 0.782*test6$steroids + 0.589*test6$prev_tnf
test6$pred16<-1-exp(-exp(sum16))

### SCENARIO 2


# temporal assessment - predictions (2) [all high dose]

sum20<- (-4.191) + 0.470*test0$age2 + 0.309*(test0$ovmean - 3.16) +
  0.484*test0$lung + 0.415*test0$renal + 0.992*test0$prev_sinf +
  0.397*test0$treatfail + 1.355*test0$steroids + 0.589*test0$prev_tnf
test0$pred20<-1-exp(-exp(sum20))

sum21<- (-4.191) + 0.470*test1$age2 + 0.309*(test1$ovmean - 3.16) +
  0.484*test1$lung + 0.415*test1$renal + 0.992*test1$prev_sinf +
  0.397*test1$treatfail + 1.355*test1$steroids + 0.589*test1$prev_tnf
test1$pred21<-1-exp(-exp(sum21))

sum22<- (-4.191) + 0.470*test2$age2 + 0.309*(test2$ovmean - 3.16) +
  0.484*test2$lung + 0.415*test2$renal + 0.992*test2$prev_sinf +
  0.397*test2$treatfail + 1.355*test2$steroids + 0.589*test2$prev_tnf
test2$pred22<-1-exp(-exp(sum22))

sum23<- (-4.191) + 0.470*test3$age2 + 0.309*(test3$ovmean - 3.16) +
  0.484*test3$lung + 0.415*test3$renal + 0.992*test3$prev_sinf +
  0.397*test3$treatfail + 1.355*test3$steroids + 0.589*test3$prev_tnf
test3$pred23<-1-exp(-exp(sum23))

sum24<- (-4.191) + 0.470*test4$age2 + 0.309*(test4$ovmean - 3.16) +
  0.484*test4$lung + 0.415*test4$renal + 0.992*test4$prev_sinf +
  0.397*test4$treatfail + 1.355*test4$steroids + 0.589*test4$prev_tnf
test4$pred24<-1-exp(-exp(sum24))

sum25<- (-4.191) + 0.470*test5$age2 + 0.309*(test5$ovmean - 3.16) +
  0.484*test5$lung + 0.415*test5$renal + 0.992*test5$prev_sinf +
  0.397*test5$treatfail + 1.355*test5$steroids + 0.589*test5$prev_tnf
test5$pred25<-1-exp(-exp(sum25))

sum26<- (-4.191) + 0.470*test6$age2 + 0.309*(test6$ovmean - 3.16) +
  0.484*test6$lung + 0.415*test6$renal + 0.992*test6$prev_sinf +
  0.397*test6$treatfail + 1.355*test6$steroids + 0.589*test6$prev_tnf
test6$pred26<-1-exp(-exp(sum26))


### SCENARIO 3 

sum30<- (-4.191) + 0.470*test0$age2 + 0.309*(test0$ovmean - 3.16) +
  0.484*test0$lung + 0.415*test0$renal + 0.992*test0$prev_sinf +
  0.397*test0$treatfail + 1.355*test0$steroids2 + 0.782*test0$lowsteroids + 
  0.589*test0$prev_tnf
test0$pred30<-1-exp(-exp(sum30))

sum31<- (-4.191) + 0.470*test1$age2 + 0.309*(test1$ovmean - 3.16) +
  0.484*test1$lung + 0.415*test1$renal + 0.992*test1$prev_sinf +
  0.397*test1$treatfail + 1.355*test1$steroids2 + 0.782*test1$lowsteroids + 
  0.589*test1$prev_tnf
test1$pred31<-1-exp(-exp(sum31))

sum32<- (-4.191) + 0.470*test2$age2 + 0.309*(test2$ovmean - 3.16) +
  0.484*test2$lung + 0.415*test2$renal + 0.992*test2$prev_sinf +
  0.397*test2$treatfail + 1.355*test2$steroids2 + 0.782*test2$lowsteroids + 
  0.589*test2$prev_tnf
test2$pred32<-1-exp(-exp(sum32))

sum33<- (-4.191) + 0.470*test3$age2 + 0.309*(test3$ovmean - 3.16) +
  0.484*test3$lung + 0.415*test3$renal + 0.992*test3$prev_sinf +
  0.397*test3$treatfail + 1.355*test3$steroids2 + 0.782*test3$lowsteroids + 
  0.589*test3$prev_tnf
test3$pred33<-1-exp(-exp(sum33))

sum34<- (-4.191) + 0.470*test4$age2 + 0.309*(test4$ovmean - 3.16) +
  0.484*test4$lung + 0.415*test4$renal + 0.992*test4$prev_sinf +
  0.397*test4$treatfail + 1.355*test4$steroids2 + 0.782*test4$lowsteroids + 
  0.589*test4$prev_tnf
test4$pred34<-1-exp(-exp(sum34))

sum35<- (-4.191) + 0.470*test5$age2 + 0.309*(test5$ovmean - 3.16) +
  0.484*test5$lung + 0.415*test5$renal + 0.992*test5$prev_sinf +
  0.397*test5$treatfail + 1.355*test5$steroids2 + 0.782*test5$lowsteroids + 
  0.589*test5$prev_tnf
test5$pred35<-1-exp(-exp(sum35))

sum36<- (-4.191) + 0.470*test6$age2 + 0.309*(test6$ovmean - 3.16) +
  0.484*test6$lung + 0.415*test6$renal + 0.992*test6$prev_sinf +
  0.397*test6$treatfail + 1.355*test6$steroids2 + 0.782*test6$lowsteroids + 
  0.589*test6$prev_tnf
test6$pred36<-1-exp(-exp(sum36))

#### NUMERICAL MEASURES

# Brier score

### SCENARIO 1


BS0<-0
BS0<-BrierScore(resp= test0$event, pred = test0$pred10, 
                   scaled = FALSE)

BS1<-0
BS1<-BrierScore(resp= test1$event, pred = test1$pred11, 
                   scaled = FALSE)

BS2<-0
BS2<-BrierScore(resp= test2$event, pred = test2$pred12, 
                   scaled = FALSE)

BS3<-0
BS3<-BrierScore(resp= test3$event, pred = test3$pred13, 
                   scaled = FALSE)

BS4<-0
BS4<-BrierScore(resp= test4$event, pred = test4$pred14, 
                   scaled = FALSE)

BS5<-0
BS5<-BrierScore(resp= test5$event, pred = test5$pred15, 
                   scaled = FALSE)

BS6<-0
BS6<-BrierScore(resp= test6$event, pred = test6$pred16, 
                   scaled = FALSE)

bs.one<-c(BS0, BS1, BS2, BS3, BS4, BS5, BS6)

#Calibration-in-the-large

### SCENARIO 1

m1 <- glm(test0$event~offset(sum10), family="binomial")
Clge0 <- as.numeric(m1$coef)
x10<-confint(m1)
Clge0.L<-as.numeric(x10[1])
Clge0.U<-as.numeric(x10[2])


m2 <- glm(test1$event~offset(sum11), family="binomial")
Clge1 <- as.numeric(m2$coef)
x11<-confint(m2)
Clge1.L<-as.numeric(x11[1])
Clge1.U<-as.numeric(x11[2])


m3 <- glm(test2$event~offset(sum12), family="binomial")
Clge2 <- as.numeric(m3$coef)
x12<-confint(m3)
Clge2.L<-as.numeric(x12[1])
Clge2.U<-as.numeric(x12[2])

m4 <- glm(test3$event~offset(sum13), family="binomial")
Clge3 <- as.numeric(m4$coef)
x13<-confint(m4)
Clge3.L<-as.numeric(x13[1])
Clge3.U<-as.numeric(x13[2])

m5 <- glm(test4$event~offset(sum14), family="binomial")
Clge4 <- as.numeric(m5$coef)
x14<-confint(m5)
Clge4.L<-as.numeric(x14[1])
Clge4.U<-as.numeric(x14[2])

m6 <- glm(test5$event~offset(sum15), family="binomial")
Clge5 <- as.numeric(m6$coef)
x15<-confint(m6)
Clge5.L<-as.numeric(x15[1])
Clge5.U<-as.numeric(x15[2])

m7 <- glm(test6$event~offset(sum16), family="binomial")
Clge6 <- as.numeric(m7$coef)
x16<-confint(m7)
Clge6.L<-as.numeric(x16[1])
Clge6.U<-as.numeric(x16[2])

clarge.one<-c(Clge0, Clge1, Clge2, Clge3, Clge4, Clge5, Clge6)
clargeL.one<-c(Clge0.L, Clge1.L, Clge2.L, Clge3.L, Clge4.L, Clge5.L, Clge6.L)
clargeU.one<-c(Clge0.U, Clge1.U, Clge2.U, Clge3.U, Clge4.U, Clge5.U, Clge6.U)

# C-statistic / AUC

### SCENARIO 1

C1 <- roc(test0$event~test0$pred10,ci=TRUE)
c0<-as.numeric(C1$auc)
c0.L<-as.numeric(C1$ci[1])
c0.U<-as.numeric(C1$ci[3])

C2 <- roc(test1$event~test1$pred11,ci=TRUE)
c1<-as.numeric(C2$auc)
c1.L<-as.numeric(C2$ci[1])
c1.U<-as.numeric(C2$ci[3])

C3 <- roc(test2$event~test2$pred12,ci=TRUE)
c2<-as.numeric(C3$auc)
c2.L<-as.numeric(C3$ci[1])
c2.U<-as.numeric(C3$ci[3])

C4 <- roc(test3$event~test3$pred13,ci=TRUE)
c3<-as.numeric(C4$auc)
c3.L<-as.numeric(C4$ci[1])
c3.U<-as.numeric(C4$ci[3])

C5 <- roc(test4$event~test4$pred14,ci=TRUE)
c4<-as.numeric(C5$auc)
c4.L<-as.numeric(C5$ci[1])
c4.U<-as.numeric(C5$ci[3])

C6 <- roc(test5$event~test5$pred15,ci=TRUE)
c5<-as.numeric(C6$auc)
c5.L<-as.numeric(C6$ci[1])
c5.U<-as.numeric(C6$ci[3])

C7 <- roc(test6$event~test6$pred16,ci=TRUE)
c6<-as.numeric(C7$auc)
c6.L<-as.numeric(C7$ci[1])
c6.U<-as.numeric(C7$ci[3])

c.one<-c(c0, c1, c2, c3, c4, c5, c6)
cL.one<-c(c0.L, c1.L, c2.L, c3.L, c4.L, c5.L, c6.L)
cU.one<-c(c0.U, c1.U, c2.U, c3.U, c4.U, c5.U, c6.U)


RRS_TA_1 <- matrix(NA, nrow= 7, ncol = 10)
RRS_TA_1[,1] <- rep("RRS", times = 7)
RRS_TA_1[,2] <- rep(1, times = 7)
RRS_TA_1[,3] <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
RRS_TA_1[,4] <- c.one
RRS_TA_1[,5] <- cL.one
RRS_TA_1[,6] <- cU.one
RRS_TA_1[,7] <- clarge.one
RRS_TA_1[,8] <- clargeL.one
RRS_TA_1[,9] <- clargeU.one
RRS_TA_1[,10] <- bs.one
RRS_TA_1 <- as.data.frame(RRS_TA_1)
names(RRS_TA_1)<- c("Mod", "Scenario", "LT", "C", "C.low", "C.high", "Clarge", "Clarge.low", "Clarge.high", "Brier")


### SCENARIO 2


# Brier score


BS0<-0
BS0<-BrierScore(resp= test0$event, pred = test0$pred20, 
                scaled = FALSE)

BS1<-0
BS1<-BrierScore(resp= test1$event, pred = test1$pred21, 
                scaled = FALSE)

BS2<-0
BS2<-BrierScore(resp= test2$event, pred = test2$pred22, 
                scaled = FALSE)

BS3<-0
BS3<-BrierScore(resp= test3$event, pred = test3$pred23, 
                scaled = FALSE)

BS4<-0
BS4<-BrierScore(resp= test4$event, pred = test4$pred24, 
                scaled = FALSE)

BS5<-0
BS5<-BrierScore(resp= test5$event, pred = test5$pred25, 
                scaled = FALSE)

BS6<-0
BS6<-BrierScore(resp= test6$event, pred = test6$pred26, 
                scaled = FALSE)

bs.one<-c(BS0, BS1, BS2, BS3, BS4, BS5, BS6)

#Calibration-in-the-large


m1 <- glm(test0$event~offset(sum20), family="binomial")
Clge0 <- as.numeric(m1$coef)
x10<-confint(m1)
Clge0.L<-as.numeric(x10[1])
Clge0.U<-as.numeric(x10[2])


m2 <- glm(test1$event~offset(sum21), family="binomial")
Clge1 <- as.numeric(m2$coef)
x11<-confint(m2)
Clge1.L<-as.numeric(x11[1])
Clge1.U<-as.numeric(x11[2])


m3 <- glm(test2$event~offset(sum22), family="binomial")
Clge2 <- as.numeric(m3$coef)
x12<-confint(m3)
Clge2.L<-as.numeric(x12[1])
Clge2.U<-as.numeric(x12[2])

m4 <- glm(test3$event~offset(sum23), family="binomial")
Clge3 <- as.numeric(m4$coef)
x13<-confint(m4)
Clge3.L<-as.numeric(x13[1])
Clge3.U<-as.numeric(x13[2])

m5 <- glm(test4$event~offset(sum24), family="binomial")
Clge4 <- as.numeric(m5$coef)
x14<-confint(m5)
Clge4.L<-as.numeric(x14[1])
Clge4.U<-as.numeric(x14[2])

m6 <- glm(test5$event~offset(sum25), family="binomial")
Clge5 <- as.numeric(m6$coef)
x15<-confint(m6)
Clge5.L<-as.numeric(x15[1])
Clge5.U<-as.numeric(x15[2])

m7 <- glm(test6$event~offset(sum26), family="binomial")
Clge6 <- as.numeric(m7$coef)
x16<-confint(m7)
Clge6.L<-as.numeric(x16[1])
Clge6.U<-as.numeric(x16[2])

clarge.one<-c(Clge0, Clge1, Clge2, Clge3, Clge4, Clge5, Clge6)
clargeL.one<-c(Clge0.L, Clge1.L, Clge2.L, Clge3.L, Clge4.L, Clge5.L, Clge6.L)
clargeU.one<-c(Clge0.U, Clge1.U, Clge2.U, Clge3.U, Clge4.U, Clge5.U, Clge6.U)

# C-statistic / AUC

C1 <- roc(test0$event~test0$pred20,ci=TRUE)
c0<-as.numeric(C1$auc)
c0.L<-as.numeric(C1$ci[1])
c0.U<-as.numeric(C1$ci[3])

C2 <- roc(test1$event~test1$pred21,ci=TRUE)
c1<-as.numeric(C2$auc)
c1.L<-as.numeric(C2$ci[1])
c1.U<-as.numeric(C2$ci[3])

C3 <- roc(test2$event~test2$pred22,ci=TRUE)
c2<-as.numeric(C3$auc)
c2.L<-as.numeric(C3$ci[1])
c2.U<-as.numeric(C3$ci[3])

C4 <- roc(test3$event~test3$pred23,ci=TRUE)
c3<-as.numeric(C4$auc)
c3.L<-as.numeric(C4$ci[1])
c3.U<-as.numeric(C4$ci[3])

C5 <- roc(test4$event~test4$pred24,ci=TRUE)
c4<-as.numeric(C5$auc)
c4.L<-as.numeric(C5$ci[1])
c4.U<-as.numeric(C5$ci[3])

C6 <- roc(test5$event~test5$pred25,ci=TRUE)
c5<-as.numeric(C6$auc)
c5.L<-as.numeric(C6$ci[1])
c5.U<-as.numeric(C6$ci[3])

C7 <- roc(test6$event~test6$pred26,ci=TRUE)
c6<-as.numeric(C7$auc)
c6.L<-as.numeric(C7$ci[1])
c6.U<-as.numeric(C7$ci[3])

c.one<-c(c0, c1, c2, c3, c4, c5, c6)
cL.one<-c(c0.L, c1.L, c2.L, c3.L, c4.L, c5.L, c6.L)
cU.one<-c(c0.U, c1.U, c2.U, c3.U, c4.U, c5.U, c6.U)


RRS_TA_2 <- matrix(NA, nrow= 7, ncol = 10)
RRS_TA_2[,1] <- rep("RRS", times = 7)
RRS_TA_2[,2] <- rep(2, times = 7)
RRS_TA_2[,3] <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
RRS_TA_2[,4] <- c.one
RRS_TA_2[,5] <- cL.one
RRS_TA_2[,6] <- cU.one
RRS_TA_2[,7] <- clarge.one
RRS_TA_2[,8] <- clargeL.one
RRS_TA_2[,9] <- clargeU.one
RRS_TA_2[,10] <- bs.one
RRS_TA_2 <- as.data.frame(RRS_TA_2)
names(RRS_TA_2)<- c("Mod", "Scenario", "LT", "C", "C.low", "C.high", "Clarge", "Clarge.low", "Clarge.high", "Brier")


### SCENARIO 3


# Brier score

BS0<-0
BS0<-BrierScore(resp= test0$event, pred = test0$pred30, 
                scaled = FALSE)

BS1<-0
BS1<-BrierScore(resp= test1$event, pred = test1$pred31, 
                scaled = FALSE)

BS2<-0
BS2<-BrierScore(resp= test2$event, pred = test2$pred32, 
                scaled = FALSE)

BS3<-0
BS3<-BrierScore(resp= test3$event, pred = test3$pred33, 
                scaled = FALSE)

BS4<-0
BS4<-BrierScore(resp= test4$event, pred = test4$pred34, 
                scaled = FALSE)

BS5<-0
BS5<-BrierScore(resp= test5$event, pred = test5$pred35, 
                scaled = FALSE)

BS6<-0
BS6<-BrierScore(resp= test6$event, pred = test6$pred36, 
                scaled = FALSE)

bs.one<-c(BS0, BS1, BS2, BS3, BS4, BS5, BS6)

#Calibration-in-the-large


m1 <- glm(test0$event~offset(sum30), family="binomial")
Clge0 <- as.numeric(m1$coef)
x10<-confint(m1)
Clge0.L<-as.numeric(x10[1])
Clge0.U<-as.numeric(x10[2])


m2 <- glm(test1$event~offset(sum31), family="binomial")
Clge1 <- as.numeric(m2$coef)
x11<-confint(m2)
Clge1.L<-as.numeric(x11[1])
Clge1.U<-as.numeric(x11[2])


m3 <- glm(test2$event~offset(sum32), family="binomial")
Clge2 <- as.numeric(m3$coef)
x12<-confint(m3)
Clge2.L<-as.numeric(x12[1])
Clge2.U<-as.numeric(x12[2])

m4 <- glm(test3$event~offset(sum33), family="binomial")
Clge3 <- as.numeric(m4$coef)
x13<-confint(m4)
Clge3.L<-as.numeric(x13[1])
Clge3.U<-as.numeric(x13[2])

m5 <- glm(test4$event~offset(sum34), family="binomial")
Clge4 <- as.numeric(m5$coef)
x14<-confint(m5)
Clge4.L<-as.numeric(x14[1])
Clge4.U<-as.numeric(x14[2])

m6 <- glm(test5$event~offset(sum35), family="binomial")
Clge5 <- as.numeric(m6$coef)
x15<-confint(m6)
Clge5.L<-as.numeric(x15[1])
Clge5.U<-as.numeric(x15[2])

m7 <- glm(test6$event~offset(sum36), family="binomial")
Clge6 <- as.numeric(m7$coef)
x16<-confint(m7)
Clge6.L<-as.numeric(x16[1])
Clge6.U<-as.numeric(x16[2])

clarge.one<-c(Clge0, Clge1, Clge2, Clge3, Clge4, Clge5, Clge6)
clargeL.one<-c(Clge0.L, Clge1.L, Clge2.L, Clge3.L, Clge4.L, Clge5.L, Clge6.L)
clargeU.one<-c(Clge0.U, Clge1.U, Clge2.U, Clge3.U, Clge4.U, Clge5.U, Clge6.U)

# C-statistic / AUC


C1 <- roc(test0$event~test0$pred30,ci=TRUE)
c0<-as.numeric(C1$auc)
c0.L<-as.numeric(C1$ci[1])
c0.U<-as.numeric(C1$ci[3])

C2 <- roc(test1$event~test1$pred31,ci=TRUE)
c1<-as.numeric(C2$auc)
c1.L<-as.numeric(C2$ci[1])
c1.U<-as.numeric(C2$ci[3])

C3 <- roc(test2$event~test2$pred32,ci=TRUE)
c2<-as.numeric(C3$auc)
c2.L<-as.numeric(C3$ci[1])
c2.U<-as.numeric(C3$ci[3])

C4 <- roc(test3$event~test3$pred33,ci=TRUE)
c3<-as.numeric(C4$auc)
c3.L<-as.numeric(C4$ci[1])
c3.U<-as.numeric(C4$ci[3])

C5 <- roc(test4$event~test4$pred34,ci=TRUE)
c4<-as.numeric(C5$auc)
c4.L<-as.numeric(C5$ci[1])
c4.U<-as.numeric(C5$ci[3])

C6 <- roc(test5$event~test5$pred35,ci=TRUE)
c5<-as.numeric(C6$auc)
c5.L<-as.numeric(C6$ci[1])
c5.U<-as.numeric(C6$ci[3])

C7 <- roc(test6$event~test6$pred36,ci=TRUE)
c6<-as.numeric(C7$auc)
c6.L<-as.numeric(C7$ci[1])
c6.U<-as.numeric(C7$ci[3])

c.one<-c(c0, c1, c2, c3, c4, c5, c6)
cL.one<-c(c0.L, c1.L, c2.L, c3.L, c4.L, c5.L, c6.L)
cU.one<-c(c0.U, c1.U, c2.U, c3.U, c4.U, c5.U, c6.U)


RRS_TA_3 <- matrix(NA, nrow= 7, ncol = 10)
RRS_TA_3[,1] <- rep("RRS", times = 7)
RRS_TA_3[,2] <- rep(3, times = 7)
RRS_TA_3[,3] <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
RRS_TA_3[,4] <- c.one
RRS_TA_3[,5] <- cL.one
RRS_TA_3[,6] <- cU.one
RRS_TA_3[,7] <- clarge.one
RRS_TA_3[,8] <- clargeL.one
RRS_TA_3[,9] <- clargeU.one
RRS_TA_3[,10] <- bs.one
RRS_TA_3 <- as.data.frame(RRS_TA_3)
names(RRS_TA_3)<- c("Mod", "Scenario", "LT", "C", "C.low", "C.high", "Clarge", "Clarge.low", "Clarge.high", "Brier")


RRS_TA <- rbind(RRS_TA_1, RRS_TA_2, RRS_TA_3)
write.table(RRS_TA, file = "rrsTA.csv", sep = ",", col.names = c("Mod", "Scenario", "LT", "C", "C.low", "C.high", "Clarge", "Clarge.low", "Clarge.high", "Brier"),
            qmethod = "double")



# Calibration plots 

par(mfrow=c(1,3))
par(pty="s")
# SCENARIO 1

# LT = 0 

groups0 <- cut(test0$pred10,breaks=quantile(test0$pred10, 
                                          prob = c(0,0.1,
                                                   0.2,0.3,0.4,
                                                   0.5,0.6,0.7,0.8,
                                                   0.9,1.0)),
              labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test0,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred10))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[1], lwd = 2)
}
h0 <- hist(test0$pred10, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                 1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test0$event))~test0$pred10,span=1))
lines_data0 <- data.frame(test0$pred10,obs_all0)
lines_data20 <- lines_data0[order(test0$pred10),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("bottomright", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test0$pred20,breaks=quantile(test0$pred20, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test0,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred20))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[2],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[2], lwd = 2)
}
h0 <- hist(test0$pred20, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test0$event))~test0$pred20,span=1))
lines_data0 <- data.frame(test0$pred20,obs_all0)
lines_data20 <- lines_data0[order(test0$pred20),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("bottomright", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test0$pred30,breaks=quantile(test0$pred30, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test0,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred30))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[3],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[3], lwd = 2)
}
h0 <- hist(test0$pred30, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test0$event))~test0$pred30,span=1))
lines_data0 <- data.frame(test0$pred30,obs_all0)
lines_data20 <- lines_data0[order(test0$pred30),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("bottomright", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")

 
# LT= 0.5


groups0 <- cut(test1$pred11,breaks=quantile(test1$pred11, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test1,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred11))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[1], lwd = 2)
}
h0 <- hist(test1$pred11, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test1$event))~test1$pred11,span=1))
lines_data0 <- data.frame(test1$pred11,obs_all0)
lines_data20 <- lines_data0[order(test1$pred11),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("bottomright", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test1$pred21,breaks=quantile(test1$pred21, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test1,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred21))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[2],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[2], lwd = 2)
}
h0 <- hist(test1$pred21, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test1$event))~test1$pred21,span=1))
lines_data0 <- data.frame(test1$pred21,obs_all0)
lines_data20 <- lines_data0[order(test1$pred21),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("bottomright", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test1$pred31,breaks=quantile(test1$pred31, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test1,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred31))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[3],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[3], lwd = 2)
}
h0 <- hist(test1$pred31, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test1$event))~test1$pred31,span=1))
lines_data0 <- data.frame(test1$pred31,obs_all0)
lines_data20 <- lines_data0[order(test1$pred31),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("bottomright", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")



## LT= 1.0


groups0 <- cut(test2$pred12,breaks=quantile(test2$pred12, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test2,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred12))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[1], lwd = 2)
}
h0 <- hist(test2$pred12, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test2$event))~test2$pred12,span=1))
lines_data0 <- data.frame(test2$pred12,obs_all0)
lines_data20 <- lines_data0[order(test2$pred12),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test2$pred22,breaks=quantile(test2$pred22, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test2,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred22))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[2],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[2], lwd = 2)
}
h0 <- hist(test2$pred22, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test2$event))~test2$pred22,span=1))
lines_data0 <- data.frame(test2$pred22,obs_all0)
lines_data20 <- lines_data0[order(test2$pred22),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test2$pred32,breaks=quantile(test2$pred32, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test2,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred32))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[3],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[3], lwd = 2)
}
h0 <- hist(test2$pred32, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test2$event))~test2$pred32,span=1))
lines_data0 <- data.frame(test2$pred32,obs_all0)
lines_data20 <- lines_data0[order(test2$pred32),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


## LT = 1.5

groups0 <- cut(test3$pred13,breaks=quantile(test3$pred13, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test3,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred13))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[1], lwd = 2)
}
h0 <- hist(test3$pred13, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test3$event))~test3$pred13,span=1))
lines_data0 <- data.frame(test3$pred13,obs_all0)
lines_data20 <- lines_data0[order(test3$pred13),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test3$pred23,breaks=quantile(test3$pred23, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test3,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred23))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[2],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[2], lwd = 2)
}
h0 <- hist(test3$pred23, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test3$event))~test3$pred23,span=1))
lines_data0 <- data.frame(test3$pred23,obs_all0)
lines_data20 <- lines_data0[order(test3$pred23),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test3$pred33,breaks=quantile(test3$pred33, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test3,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred33))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[3],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[3], lwd = 2)
}
h0 <- hist(test3$pred33, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test3$event))~test3$pred33,span=1))
lines_data0 <- data.frame(test3$pred33,obs_all0)
lines_data20 <- lines_data0[order(test3$pred33),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


## LT = 2.0

groups0 <- cut(test4$pred14,breaks=quantile(test4$pred14, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test4,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred14))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[1], lwd = 2)
}
h0 <- hist(test4$pred14, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test4$event))~test4$pred14,span=1))
lines_data0 <- data.frame(test4$pred14,obs_all0)
lines_data20 <- lines_data0[order(test4$pred14),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test4$pred24,breaks=quantile(test4$pred24, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test4,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred24))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[2],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[2], lwd = 2)
}
h0 <- hist(test4$pred24, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test4$event))~test4$pred24,span=1))
lines_data0 <- data.frame(test4$pred24,obs_all0)
lines_data20 <- lines_data0[order(test4$pred24),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test4$pred34,breaks=quantile(test4$pred34, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test4,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred34))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[3],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[3], lwd = 2)
}
h0 <- hist(test4$pred34, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test4$event))~test4$pred34,span=1))
lines_data0 <- data.frame(test4$pred34,obs_all0)
lines_data20 <- lines_data0[order(test4$pred34),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")




## LT= 2.5

groups0 <- cut(test5$pred15,breaks=quantile(test5$pred15, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test5,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred15))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[1], lwd = 2)
}
h0 <- hist(test5$pred15, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test5$event))~test5$pred15,span=1))
lines_data0 <- data.frame(test5$pred15,obs_all0)
lines_data20 <- lines_data0[order(test5$pred15),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test5$pred25,breaks=quantile(test5$pred25, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test5,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred25))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[2],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[2], lwd = 2)
}
h0 <- hist(test5$pred25, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test5$event))~test5$pred25,span=1))
lines_data0 <- data.frame(test5$pred25,obs_all0)
lines_data20 <- lines_data0[order(test5$pred25),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test5$pred35,breaks=quantile(test5$pred35, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test5,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred35))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[3],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[3], lwd = 2)
}
h0 <- hist(test5$pred35, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test5$event))~test5$pred35,span=1))
lines_data0 <- data.frame(test5$pred35,obs_all0)
lines_data20 <- lines_data0[order(test5$pred35),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


## LT = 3

groups0 <- cut(test6$pred16,breaks=quantile(test6$pred16, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test6,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred16))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[1],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[1], lwd = 2)
}
h0 <- hist(test6$pred16, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test6$event))~test6$pred16,span=1))
lines_data0 <- data.frame(test6$pred16,obs_all0)
lines_data20 <- lines_data0[order(test6$pred16),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[1],"black",cols[1],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test6$pred26,breaks=quantile(test6$pred26, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test6,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred26))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[2],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[2], lwd = 2)
}
h0 <- hist(test6$pred26, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test6$event))~test6$pred26,span=1))
lines_data0 <- data.frame(test6$pred26,obs_all0)
lines_data20 <- lines_data0[order(test6$pred26),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[2],"black",cols[2],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")


groups0 <- cut(test6$pred36,breaks=quantile(test6$pred36, 
                                            prob = c(0,0.1,
                                                     0.2,0.3,0.4,
                                                     0.5,0.6,0.7,0.8,
                                                     0.9,1.0)),
               labels=c(1:10),include.lowest=TRUE)
gpdata0 <- cbind(test6,groups0)
obs0 <- ddply(gpdata0,~groups0,summarise,mean=mean(as.numeric(event)))[,2]
exp0 <- ddply(gpdata0,~groups0,summarise,mean=mean(pred36))
attach(gpdata0)
obsn0 <- table(event,groups0)[1,] 
lci0 = pmax(0,(obs0 - (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
uci0 = pmin(1,(obs0 + (1.96*(((obs0*(1-obs0))/obsn0)^.5))))
plot(obs0~exp0[,2],xlim=c(0,0.25),ylim=c(0,0.25),
     col=cols[3],ylab="Observed",xlab="Expected", pch = 16)
lines(c(0,0.25),c(0,0.25),lty=2)
for(i in 1:10){
  lines(c(exp0[i,2],exp0[i,2]),c(lci0[i],uci0[i]),col=cols[3], lwd = 2)
}
h0 <- hist(test6$pred36, breaks=50, plot=FALSE)
for(i in 1:length(h0$mids)){
  lines(c(h0$mids[i],h0$mids[i]),c(rep(1,length(h0$mids))[i],
                                   1-((h0$counts[i]/max(h0$counts))/10)))}
obs_all0 <- predict(loess((as.numeric(test6$event))~test6$pred36,span=1))
lines_data0 <- data.frame(test6$pred36,obs_all0)
lines_data20 <- lines_data0[order(test6$pred36),] 
lines(lines_data20[,1],lines_data20[,2],col="grey")
legend("topleft", c("Risk groups","Reference line","95% CI","Loess"),
       col=c(cols[3],"black",cols[3],"grey"),lty=c(0,2,1,1),
       pch=c(16,NA,NA,NA),bty="n")



############# Survival assessment

test0<-read.dta("baseline_model.dta")
test0$smoke<-factor(test0$smoke)
test0$prev_sinf<-0
test0$prev_tnf<-1
test0$age2<-0
test0$treatfail<-0
test0$steroids2<-0
test0$lowsteroids<-0
for (i in (1:length(test0$studyno))){
  if (test0$age[i] > 60) {test0$age2[i] = 1}
  if (test0$previous_dmards[i] > 5) {test0$treatfail[i] = 1}
  if ((test0$steroids[i] == 1) & 
      (test0$diabetes[i] == 0) & (test0$smoke[i] == 0)){
    test0$steroids2[i] = 1}
  if ((test0$steroids[i] == 1) & (test0$steroids2[i] == 0)){
    test0$lowsteroids[i] = 1}
}

test1<-read.dta("validation_6months.dta")
test1$smoke<-factor(test1$smoke)
test1$prev_tnf<-1
test1$age2<-0
test1$treatfail<-0
test1$steroids2<-0
test1$lowsteroids<-0
for (i in (1:length(test1$studyno))){
  if (test1$age[i] > 60) {test1$age2[i] = 1}
  if (test1$previous_dmards[i] > 5) {test1$treatfail[i] = 1}
  if ((test1$steroids[i] == 1) & 
      (test1$diabetes[i] == 0) & (test1$smoke[i] == 0)){
    test1$steroids2[i] = 1}
  if ((test1$steroids[i] == 1) & (test1$steroids2[i] == 0)){
    test1$lowsteroids[i] = 1}
}

test2<-read.dta("validation_12months.dta")
test2$smoke<-factor(test2$smoke)
test2$prev_tnf<-1
test2$age2<-0
test2$treatfail<-0
test2$steroids2<-0
test2$lowsteroids<-0
for (i in (1:length(test2$studyno))){
  if (test2$age[i] > 60) {test2$age2[i] = 1}
  if (test2$previous_dmards[i] > 5) {test2$treatfail[i] = 1}
  if ((test2$steroids[i] == 1) & 
      (test2$diabetes[i] == 0) & (test2$smoke[i] == 0)){
    test2$steroids2[i] = 1}
  if ((test2$steroids[i] == 1) & (test2$steroids2[i] == 0)){
    test2$lowsteroids[i] = 1}
}

test3<-read.dta("validation_18months.dta")
test3$smoke<-factor(test3$smoke)
test3$prev_tnf<-1
test3$age2<-0
test3$treatfail<-0
test3$steroids2<-0
test3$lowsteroids<-0
for (i in (1:length(test3$studyno))){
  if (test3$age[i] > 60) {test3$age2[i] = 1}
  if (test3$previous_dmards[i] > 5) {test3$treatfail[i] = 1}
  if ((test3$steroids[i] == 1) & 
      (test3$diabetes[i] == 0) & (test3$smoke[i] == 0)){
    test3$steroids2[i] = 1}
  if ((test3$steroids[i] == 1) & (test3$steroids2[i] == 0)){
    test3$lowsteroids[i] = 1}
}

test4<-read.dta("validation_24months.dta")
test4$smoke<-factor(test4$smoke)
test4$prev_tnf<-1
test4$age2<-0
test4$treatfail<-0
test4$steroids2<-0
test4$lowsteroids<-0
for (i in (1:length(test4$studyno))){
  if (test4$age[i] > 60) {test4$age2[i] = 1}
  if (test4$previous_dmards[i] > 5) {test4$treatfail[i] = 1}
  if ((test4$steroids[i] == 1) & 
      (test4$diabetes[i] == 0) & (test4$smoke[i] == 0)){
    test4$steroids2[i] = 1}
  if ((test4$steroids[i] == 1) & (test4$steroids2[i] == 0)){
    test4$lowsteroids[i] = 1}
}

test5<-read.dta("validation_30months.dta")
test5$smoke<-factor(test5$smoke)
test5$prev_tnf<-1
test5$age2<-0
test5$treatfail<-0
test5$steroids2<-0
test5$lowsteroids<-0
for (i in (1:length(test5$studyno))){
  if (test5$age[i] > 60) {test5$age2[i] = 1}
  if (test5$previous_dmards[i] > 5) {test5$treatfail[i] = 1}
  if ((test5$steroids[i] == 1) & 
      (test5$diabetes[i] == 0) & (test5$smoke[i] == 0)){
    test5$steroids2[i] = 1}
  if ((test5$steroids[i] == 1) & (test5$steroids2[i] == 0)){
    test5$lowsteroids[i] = 1}
}

test6<-read.dta("validation_36months.dta")
test6$smoke<-factor(test6$smoke)
test6$prev_tnf<-1
test6$age2<-0
test6$treatfail<-0
test6$steroids2<-0
test6$lowsteroids<-0
for (i in (1:length(test6$studyno))){
  if (test6$age[i] > 60) {test6$age2[i] = 1}
  if (test6$previous_dmards[i] > 5) {test6$treatfail[i] = 1}
  if ((test6$steroids[i] == 1) & 
      (test6$diabetes[i] == 0) & (test6$smoke[i] == 0)){
    test6$steroids2[i] = 1}
  if ((test6$steroids[i] == 1) & (test6$steroids2[i] == 0)){
    test6$lowsteroids[i] = 1 
  }
}


# survival temporal assessment - predictions (1) [all low dose]

### SCENARIO 1

sum10<- (-4.191) + 0.470*test0$age2 + 0.309*(test0$ovmean - 3.16) +
  0.484*test0$lung + 0.415*test0$renal + 0.992*test0$prev_sinf +
  0.397*test0$treatfail + 0.782*test0$steroids + 0.589*test0$prev_tnf
test0$pred<-1-exp(-exp(sum10))

sum11<- (-4.191) + 0.470*test1$age2 + 0.309*(test1$ovmean - 3.16) +
  0.484*test1$lung + 0.415*test1$renal + 0.992*test1$prev_sinf +
  0.397*test1$treatfail + 0.782*test1$steroids + 0.589*test1$prev_tnf
test1$pred<-1-exp(-exp(sum11))

sum12<- (-4.191) + 0.470*test2$age2 + 0.309*(test2$ovmean - 3.16) +
  0.484*test2$lung + 0.415*test2$renal + 0.992*test2$prev_sinf +
  0.397*test2$treatfail + 0.782*test2$steroids + 0.589*test2$prev_tnf
test2$pred<-1-exp(-exp(sum12))

sum13<- (-4.191) + 0.470*test3$age2 + 0.309*(test3$ovmean - 3.16) +
  0.484*test3$lung + 0.415*test3$renal + 0.992*test3$prev_sinf +
  0.397*test3$treatfail + 0.782*test3$steroids + 0.589*test3$prev_tnf
test3$pred<-1-exp(-exp(sum13))

sum14<- (-4.191) + 0.470*test4$age2 + 0.309*(test4$ovmean - 3.16) +
  0.484*test4$lung + 0.415*test4$renal + 0.992*test4$prev_sinf +
  0.397*test4$treatfail + 0.782*test4$steroids + 0.589*test4$prev_tnf
test4$pred<-1-exp(-exp(sum14))

sum15<- (-4.191) + 0.470*test5$age2 + 0.309*(test5$ovmean - 3.16) +
  0.484*test5$lung + 0.415*test5$renal + 0.992*test5$prev_sinf +
  0.397*test5$treatfail + 0.782*test5$steroids + 0.589*test5$prev_tnf
test5$pred<-1-exp(-exp(sum15))

sum16<- (-4.191) + 0.470*test6$age2 + 0.309*(test6$ovmean - 3.16) +
  0.484*test6$lung + 0.415*test6$renal + 0.992*test6$prev_sinf +
  0.397*test6$treatfail + 0.782*test6$steroids + 0.589*test6$prev_tnf
test6$pred<-1-exp(-exp(sum16))

test0$surv_pred0 <- 1-test0$pred
test1$surv_pred0 <- 1-test1$pred
test2$surv_pred0 <- 1-test2$pred
test3$surv_pred0 <- 1-test3$pred
test4$surv_pred0 <- 1-test4$pred
test5$surv_pred0 <- 1-test5$pred
test6$surv_pred0 <- 1-test6$pred

LT <- c(0,0.5,1,1.5,2,2.5,3)
Cstat<-0
CstatL<-0
CstatU<-0
BS<-0

## LT = 0
c0 <- coxph(Surv(timevent,event)~sum10, data=test0)
Cstat[1]<-as.numeric(summary(c0)$concordance[1])
CstatL[1]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[1]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

BS[1]<-sbrier(obj = Surv(test0$timevent,test0$event), 
              test0$surv_pred0, btime = 1)

## LT = 0.5
c0 <- coxph(Surv(timevent,event)~sum11, data=test1)
Cstat[2]<-as.numeric(summary(c0)$concordance[1])
CstatL[2]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[2]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

BS[2]<-sbrier(obj = Surv(test1$timevent,test1$event), 
              test1$surv_pred0, btime = 1.5)

## LT = 1
c0 <- coxph(Surv(timevent,event)~sum12, data=test2)
Cstat[3]<-as.numeric(summary(c0)$concordance[1])
CstatL[3]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[3]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

BS[3]<-sbrier(obj = Surv(test2$timevent,test2$event), 
              test2$surv_pred0, btime = 2)

## LT = 1.5
c0 <- coxph(Surv(timevent,event)~sum13, data=test3)
Cstat[4]<-as.numeric(summary(c0)$concordance[1])
CstatL[4]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[4]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

BS[4]<-sbrier(obj = Surv(test3$timevent,test3$event), 
              test3$surv_pred0, btime = 2.5)

## LT = 2
c0 <- coxph(Surv(timevent,event)~sum14, data=test4)
Cstat[5]<-as.numeric(summary(c0)$concordance[1])
CstatL[5]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[5]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

BS[5]<-sbrier(obj = Surv(test4$timevent,test4$event), 
              test4$surv_pred0, btime = 3)

## LT = 2.5
c0 <- coxph(Surv(timevent,event)~sum15, data=test5)
Cstat[6]<-as.numeric(summary(c0)$concordance[1])
CstatL[6]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[6]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))

# Brier
BS[6]<-sbrier(obj = Surv(test5$timevent,test5$event), 
              test5$surv_pred0, btime = 3.5)

## LT = 3
c0 <- coxph(Surv(timevent,event)~sum16, data=test6)
Cstat[7]<-as.numeric(summary(c0)$concordance[1])
CstatL[7]<-as.numeric(summary(c0)$concordance[1]
                      -(1.96*summary(c0)$concordance[2]))
CstatU[7]<-as.numeric(summary(c0)$concordance[1]
                      +(1.96*summary(c0)$concordance[2]))
BS[7]<-sbrier(obj = Surv(test6$timevent,test6$event), 
              test6$surv_pred0, btime = 4)

RRS.TA_surv<-cbind(LT, Cstat, CstatL, CstatU, BS)
RRS.TA_surv<-as.data.frame(RRS.TA_surv)
write.table(RRS.TA_surv, file = "rrsTA_surv.csv", sep = ",", col.names = NA,
            qmethod = "double")
