frogMFCCs = read.csv(file.choose(), header = T)

library(class)

#Step 1:
frogMFCCs = frogMFCCs[,1:23] #remove extra columns

#classes
CL1 = subset(frogMFCCs, frogMFCCs[,23]=='Dendrobatidae') 
CL2 = subset(frogMFCCs, frogMFCCs[,23]=='Hylidae') 
CL3 = subset(frogMFCCs, frogMFCCs[,23]=='Leptodactylidae') 

#triple Dendrobatidae cases; randomly remove 1420 Leptodactylidae cases
CL1 = rbind(CL1, CL1, CL1) 
CL3 = CL3[-sample(c(1:4420), 1420),] 

#cleaned up data
frogDATA = rbind(CL1,CL2,CL3)

#standardize
m = apply(frogDATA[,-23], 2, mean)
s = apply(frogDATA[,-23], 2, sd)

CL1 = cbind(CL1[23], t(t(sweep(CL1[,-23],2,m))/s))
CL2 = cbind(CL2[23], t(t(sweep(CL2[,-23],2,m))/s))
CL3 = cbind(CL3[23], t(t(sweep(CL3[,-23],2,m))/s))
SDATA = rbind(CL1,CL2,CL3)

#Step 2:
#1.3, discriminating power(DP)
m1 = apply(CL1[,-1], 2, mean); m2 = apply(CL2[,-1], 2, mean); m3 = apply(CL3[,-1], 2, mean)
sd1 = apply(CL1[,-1], 2, sd); sd2 = apply(CL2[,-1], 2, sd); sd3 = apply(CL3[,-1], 2, sd)

dp12 = abs(m1-m2) / sqrt((abs(sd1^2 + sd2^2))/ (nrow(CL1) + nrow(CL2)))
dp13 = abs(m1-m3) / sqrt((abs(sd1^2 + sd3^2))/ (nrow(CL1) + nrow(CL3)))
dp23 = abs(m2-m3) / sqrt((abs(sd2^2 + sd3^2))/ (nrow(CL2) + nrow(CL3)))

#top 3 for each:
sort(dp12)[20:22]
sort(dp13)[20:22]
sort(dp23)[20:22]

#t-test
t.test(CL1[,'MFCCs_.3'], CL2[,'MFCCs_.3'])
t.test(CL1[,'MFCCs_17'], CL3[,'MFCCs_17'])
t.test(CL2[,'MFCCs_19'], CL3[,'MFCCs_19'])

t.test(CL1[,'MFCCs_10'], CL3[,'MFCCs_10'])

#histograms
par(mfrow=c(2,3))
hist(CL1[,'MFCCs_.3']); hist(CL2[,'MFCCs_.3']); hist(CL3[,'MFCCs_.3'])
hist(CL1[,'MFCCs_17']); hist(CL2[,'MFCCs_17']); hist(CL3[,'MFCCs_17'])
hist(CL1[,'MFCCs_19']); hist(CL2[,'MFCCs_19']); hist(CL3[,'MFCCs_19'])
hist(CL1[,'MFCCs_10']); hist(CL2[,'MFCCs_10']); hist(CL3[,'MFCCs_10'])

#Step 3:
#assign selections ~20% of each class
r1 = sort(sample(nrow(CL1),nrow(CL1)*0.2))
r2 = sort(sample(nrow(CL2),nrow(CL2)*0.2))
r3 = sort(sample(nrow(CL3),nrow(CL3)*0.2))
#define test&training sets
testCL1 = CL1[r1,]; trainCL1=CL1[-r1,]
testCL2 = CL2[r2,]; trainCL2=CL2[-r2,]
testCL3 = CL3[r3,]; trainCL3=CL3[-r3,]

TRAINSET=rbind(trainCL1,trainCL2,trainCL3)
TESTSET=rbind(testCL1,testCL2,testCL3)

#looking for best k value
k_try = c(1,5,10,20,30,40,50,100)

trainperfK=NULL
testperfK=NULL

for (i in k_try){
  trainknn = knn(TRAINSET[,-1], TRAINSET[,-1], TRAINSET$Family, i)
  trainperfK = c(trainperfK, mean(TRAINSET$Family==trainknn))
  testknn = knn(TRAINSET[,-1], TESTSET[,-1], TRAINSET$Family, i)
  testperfK = c(testperfK, mean(TESTSET$Family==testknn))
}
overfit = trainperfK - testperfK
plot(k_try, trainperfK, type='l', ylab='performance')
lines(k_try, testperfK, col='blue')
legend("topright",c("TRAIN acc","TEST acc"),col=c("black","blue"), lwd=1, lty=c(1,1))

plot(k_try, overfit, type='l', col='red')

#looks most promising in 1-10 range:
k_try = seq(1,10)

trainperfK=NULL
testperfK=NULL

for (i in k_try){
  trainknn = knn(TRAINSET[,-1], TRAINSET[,-1], TRAINSET$Family, i)
  trainperfK = c(trainperfK, mean(TRAINSET$Family==trainknn))
  testknn = knn(TRAINSET[,-1], TESTSET[,-1], TRAINSET$Family, i)
  testperfK = c(testperfK, mean(TESTSET$Family==testknn))
}
overfit = trainperfK - testperfK

plot(k_try, trainperfK, type='l', ylab='performance')
lines(k_try, testperfK, col='blue')
legend("topright",c("TRAIN acc","TEST acc"),col=c("black","blue"), lwd=1, lty=c(1,1))
plot(overfit)
#k=5 has lowest overfit

#for best k = 5:
train = knn(TRAINSET[,-1], TRAINSET[,-1], TRAINSET$Family, k=5)
test = knn(TRAINSET[,-1], TESTSET[,-1], TRAINSET$Family, k=5)

#confusion matrices 
trainconf = table(TRAINSET$Family,train)
testconf = table(TESTSET$Family,test)
#global acc
sum(diag(trainconf))/sum(trainconf)
sum(diag(testconf))/sum(testconf)

#conf.matrices in %'s
trainconf = (trainconf/apply(trainconf,1,sum))[-1,-1]
testconf = (testconf/apply(testconf,1,sum))[-1,-1]

#check confidence intervals
#training set:
ptrain=as.numeric(diag(trainconf))
sigmaCL1train=sqrt(ptrain[1]*(1-ptrain[1])/nrow(trainCL1))
sigmaCL2train=sqrt(ptrain[2]*(1-ptrain[2])/nrow(trainCL2))
sigmaCL3train=sqrt(ptrain[3]*(1-ptrain[3])/nrow(trainCL3))
#90% confidence interval
intervalCL1train=c(ptrain[1]-sigmaCL1train*qnorm(1-0.1/2), ptrain[1]+sigmaCL1train*qnorm(1-0.1/2))
intervalCL2train=c(ptrain[2]-sigmaCL1train*qnorm(1-0.1/2), ptrain[2]+sigmaCL2train*qnorm(1-0.1/2))
intervalCL3train=c(ptrain[3]-sigmaCL1train*qnorm(1-0.05/2), ptrain[3]+sigmaCL3train*qnorm(1-0.05/2))

#test set:
ptest=as.numeric(diag(testconf))
sigmaCL1test=sqrt(ptest[1]*(1-ptest[1])/nrow(testCL1))
sigmaCL2test=sqrt(ptest[2]*(1-ptest[2])/nrow(testCL2))
sigmaCL3test=sqrt(ptest[3]*(1-ptest[3])/nrow(testCL3))
#90% confidence interval
intervalCL1test=c(ptest[1]-sigmaCL1test*qnorm(1-0.1/2), ptest[1]+sigmaCL1test*qnorm(1-0.1/2))
intervalCL2test=c(ptest[2]-sigmaCL1test*qnorm(1-0.1/2), ptest[2]+sigmaCL2test*qnorm(1-0.1/2))
intervalCL3test=c(ptest[3]-sigmaCL1test*qnorm(1-0.05/2), ptest[3]+sigmaCL3test*qnorm(1-0.05/2))

incorrectTrain = TRAINSET[which(train != TRAINSET$Family),]
correctTrain = TRAINSET[which(train == TRAINSET$Family),]
incorrectTest = TESTSET[which(test != TESTSET$Family),]
correctTest = TESTSET[which(test == TESTSET$Family),]

#Part 4, PCA
#correlation matrix of standardized features
CORR=cor(SDATA[,-1])

#eigenvalues
lamd = eigen(CORR)
ev = lamd$values
W = lamd$vectors

plot(ev, ylab='eigenvalue', xlab='r', pch=20)

#Getting PVE for r=1 to 22
PVE = function(r){
  return (sum(ev[1:r])/22)
}
PVEs = NULL
for(r in 1:22) PVEs=c(PVEs, PVE(r))
pveTable = as.table(setNames((PVEs),seq(1:22)))
#PVE > 95% at r=12

W.t = t(W) #transpose W

Y = as.matrix(SDATA[,-1]) %*%  W.t 
# take only the first 12 components
Y = Y[,1:12]

#KNN with these new features
#(Repeating the previous steps in part3, but for new features):
newDATA = cbind(SDATA[,1], as.data.frame(Y))
names(newDATA)[1] = 'Family'

#classes
CL1 = subset(newDATA, newDATA[,1]=='Dendrobatidae') 
CL2 = subset(newDATA, newDATA[,1]=='Hylidae') 
CL3 = subset(newDATA, newDATA[,1]=='Leptodactylidae') 

#assign selections ~20% of each class
r1 = sort(sample(nrow(CL1),nrow(CL1)*0.2))
r2 = sort(sample(nrow(CL2),nrow(CL2)*0.2))
r3 = sort(sample(nrow(CL3),nrow(CL3)*0.2))
#define test&training sets
testCL1 = CL1[r1,]; trainCL1=CL1[-r1,]
testCL2 = CL2[r2,]; trainCL2=CL2[-r2,]
testCL3 = CL3[r3,]; trainCL3=CL3[-r3,]

TRAINSET=rbind(trainCL1,trainCL2,trainCL3)
TESTSET=rbind(testCL1,testCL2,testCL3)

#looking for best k value
k_try = c(1,5,10,20,30,40,50,100)

trainperfK=NULL
testperfK=NULL

for (i in k_try){
  trainknn = knn(TRAINSET[,-1], TRAINSET[,-1], TRAINSET$Family, i)
  trainperfK = c(trainperfK, mean(TRAINSET$Family==trainknn))
  testknn = knn(TRAINSET[,-1], TESTSET[,-1], TRAINSET$Family, i)
  testperfK = c(testperfK, mean(TESTSET$Family==testknn))
}
overfit = trainperfK - testperfK
plot(k_try, trainperfK, type='l', ylab='performance')
lines(k_try, testperfK, col='blue')
legend("topright",c("TRAIN acc","TEST acc"),col=c("black","blue"), lwd=1, lty=c(1,1))

plot(k_try, overfit, type='l', col='red')

#looks most promising in 1-10 range:
k_try = seq(1,10)

trainperfK=NULL
testperfK=NULL

for (i in k_try){
  trainknn = knn(TRAINSET[,-1], TRAINSET[,-1], TRAINSET$Family, i)
  trainperfK = c(trainperfK, mean(TRAINSET$Family==trainknn))
  testknn = knn(TRAINSET[,-1], TESTSET[,-1], TRAINSET$Family, i)
  testperfK = c(testperfK, mean(TESTSET$Family==testknn))
}
overfit = trainperfK - testperfK

plot(k_try, trainperfK, type='l', ylab='performance')
lines(k_try, testperfK, col='blue')
legend("topright",c("TRAIN acc","TEST acc"),col=c("black","blue"), lwd=1, lty=c(1,1))
plot(overfit)

#k=3 has lowest overfit
train = knn(TRAINSET[,-1], TRAINSET[,-1], TRAINSET$Family, k=3)
test = knn(TRAINSET[,-1], TESTSET[,-1], TRAINSET$Family, k=3)

#confusion matrices 
trainconf = table(TRAINSET$Family,train)
testconf = table(TESTSET$Family,test)
#global acc
sum(diag(trainconf))/sum(trainconf)
sum(diag(testconf))/sum(testconf)

#conf.matrices in %'s
trainconf = (trainconf/apply(trainconf,1,sum))[-1,-1]
testconf = (testconf/apply(testconf,1,sum))[-1,-1]
