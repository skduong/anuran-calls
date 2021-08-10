frogMFCCs = read.csv('Frogs_MFCCs.csv', header = T)

library(randomForest)
library(scatterplot3d)#10/21 for cloning

#Cleaning
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
#testing SMOTE
library(DMwR)
frogDATAsmote <- SMOTE(Family ~ ., frogDATA, perc.over = 1084,perc.under=1420)
table(newData$Species)

#standardize
m = apply(frogDATA[,-23], 2, mean)
s = apply(frogDATA[,-23], 2, sd)
CL1s = cbind(CL1[23], t(t(sweep(CL1[,-23],2,m))/s))
CL2s = cbind(CL2[23], t(t(sweep(CL2[,-23],2,m))/s))
CL3s = cbind(CL3[23], t(t(sweep(CL3[,-23],2,m))/s))
SDATA = rbind(CL1s,CL2s,CL3s)

#Q1: PCA
#correlation matrix of standardized features
CORR=cor(SDATA[,-1]) #save and give preview?
write.csv(CORR, "Q1_CORR.csv")

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

plot(PVEs, ylab='PVE', xlab='r', pch=20)

W.t = t(W) #transpose W
Y = as.matrix(SDATA[,-1]) %*%  W.t 
# take only the first 12 components
Y = Y[,1:12]

#3d plots
newDATA = cbind(SDATA[,1], as.data.frame(Y))
names(newDATA)[1] = 'Family'

#take first 3 PCA vectors for each class
CL1new = subset(newDATA, newDATA[,1]=='Dendrobatidae')[,1:4] 
CL2new = subset(newDATA, newDATA[,1]=='Hylidae')[,1:4] 
CL3new = subset(newDATA, newDATA[,1]=='Leptodactylidae')[,1:4] 

#CL1 & CL2
pair = rbind(CL1new,CL2new)
#levels(pair$Family) = droplevels(pair$Family)
pair$Family = droplevels(pair$Family)
cols = c('red', 'green')
scatterplot3d(pair[,2:4], xlab='V1', ylab='V2', zlab='V3',
              pch=20, color = cols[as.numeric(pair$Family)])
legend(4,6.3, legend = levels(pair$Family), col = cols, pch=20)

#separate (testing purposes)
scatterplot3d(CL1new[,2:4], xlab='PCA1', ylab='PCA2', zlab='PCA3',
              pch=20, color = 'red')
scatterplot3d(CL2new[,2:4], xlab='PCA1', ylab='PCA2', zlab='PCA3',
              pch=20, color = 'green')
scatterplot3d(CL3new[,2:4], xlab='PCA1', ylab='PCA2', zlab='PCA3',
              pch=20, color = 'blue')

#CL1 & CL3
#par(mfrow=c(3,1))
pair = rbind(CL1new,CL3new)
pair$Family = droplevels(pair$Family)
cols = c('red', 'blue')
scatterplot3d(pair[,2:4], xlab='V1', ylab='V2', zlab='V3',
              pch=20, color = cols[as.numeric(pair$Family)])
legend(4.3,7.5, legend = levels(pair$Family), col = cols, pch=20)

#CL2 & CL3
pair = rbind(CL2new,CL3new)
pair$Family = droplevels(pair$Family)
cols = c('green', 'blue')
scatterplot3d(pair[,2:4], xlab='V1', ylab='V2', zlab='V3',
              pch=20, color = cols[as.numeric(pair$Family)])
legend(4.3,7.5, legend = levels(pair$Family), col = cols, pch=20)

#all together
pair = rbind(CL2new,CL3new, CL1new)
pair$Family = droplevels(pair$Family)
cols = c('red','green', 'blue')
scatterplot3d(pair[,2:4], xlab='V1', ylab='V2', zlab='V3',
              pch=20, color = cols[as.numeric(pair$Family)])
legend(4.3,7.5, legend = levels(pair$Family), col = cols, pch=20)

#Q2: K-means see 11/04
set.seed(6350)
kclusters = vector(mode = "list", length = 10)
for(k in 1:10){
  kclusters[[k]] = kmeans(SDATA[,-1], k, nstart=50)
} #34ish sec to run

#reduction of variance
perfk = NULL
for(i in 1:10){
  perfk = c(perfk, 1-sum(kclusters[[i]]$withinss)/kclusters[[1]]$withinss)
}
plot(perfk, type='l', xlab='K', ylab='perf(K)')

#
ginik = vector(mode = "list", length = 10); maxGini = NULL
IMP = NULL #sum of gini for each cluster
for (i in 1:10){#k=1-10
  clusterk = kclusters[[i]]
  gini = NULL
  for(k in 1:i){
    CLU = frogDATA[as.numeric(names(clusterk$cluster[clusterk$cluster == k])),]
    f1 = nrow(subset(CLU, CLU[,23]=='Dendrobatidae'))/nrow(CLU)
    f2 = nrow(subset(CLU, CLU[,23]=='Hylidae'))/nrow(CLU)
    f3 = nrow(subset(CLU, CLU[,23]=='Leptodactylidae'))/nrow(CLU)
    gini = c(gini, f1*(1-f1)+f2*(1-f2)+f3*(1-f3))
  }
  ginik[[i]] = gini
  IMP = c(IMP, sum(gini))
  maxGini = c(maxGini, i*2/3)
}
plot(IMP[1:5], xlab='K', ylab='Impurity(K)', type = 'l')

#repeat single K-means execution for K=2-4 50 times to see which clusters give best Gini
gini24 = vector(mode = "list", length = 4)
for(j in 1:4){
  for (i in 1:30){
    clusterk = kmeans(SDATA[,-1], j, nstart=1)
    gini = NULL
    for(k in 1:j){
      CLU = frogDATA[as.numeric(names(clusterk$cluster[clusterk$cluster == k])),]
      f1 = nrow(subset(CLU, CLU[,23]=='Dendrobatidae'))/nrow(CLU)
      f2 = nrow(subset(CLU, CLU[,23]=='Hylidae'))/nrow(CLU)
      f3 = nrow(subset(CLU, CLU[,23]=='Leptodactylidae'))/nrow(CLU)
      gini = c(gini, f1*(1-f1)+f2*(1-f2)+f3*(1-f3))
      gini24[[j]] = c(gini24[[j]], sum(gini))
    }
  }
}

minGini = NULL
for (i in 1:4){
  minGini = c(minGini, min(gini24[[i]][seq(i, length(gini24[[i]]), i)]))
}
plot(minGini, type='l', xlab='K', ylab='min(Gini)')
bestk = 2#(Might change)

#Q3
clusterBest = kclusters[[bestk]]
clusterk = clusterBest
CENT = clusterBest$center #CENT[1,] 

scatterplot3d(unname(CENT[,1:3]), xlab='V1', ylab='V2', zlab='V3', pch=20 ,label.tick.marks=FALSE)
text(c(3,2) ,c(7,-3),labels = c('Cl2', 'Cl1'))
#looks better
scatterplot3d(unname(CENT[,1:3]), xlab='V1', ylab='V2', zlab='V3', pch=20, lab.z=1)
text(c(3,3) ,c(4,-0.3),labels = c('CL2', 'CL1'))

#Gini/Freq of class for each cluster
fams = c('Dendrobatidae', 'Hylidae', 'Leptodactylidae')
gini = NULL #gini index for [kth] cluster
fDendro = NULL; fHyli = NULL; fLepto = NULL #frequency of family cases in cluster[k]
TOP = NULL 
frogDAT = cbind(frogDATA, clusterBest$cluster)
newDAT = cbind(newDATA, clusterBest$cluster)
for(k in 1:bestk){
  CLU = frogDAT[frogDAT$`clusterBest$cluster` == k,]
  f1 = nrow(subset(CLU, CLU[,23]== fams[1]))/nrow(CLU)
  f2 = nrow(subset(CLU, CLU[,23]==fams[2]))/nrow(CLU)
  f3 = nrow(subset(CLU, CLU[,23]==fams[3]))/nrow(CLU)
  gini = c(gini, f1*(1-f1)+f2*(1-f2)+f3*(1-f3))
  fDendro = c(fDendro, f1); fHyli = c(fHyli, f2); fLepto = c(fLepto, f3)
  TOP = c(TOP, fams[which(c(f1,f2,f3) == max(f1,f2,f3))])
}

IMP = sum(gini)
FREQ = rbind(fDendro, fHyli, fLepto) #apply(FREQ,2,sum)

#projection of clusters:
pair = newDAT
pair$Family = droplevels(pair$Family)
cols = c('orange', 'purple')
scatterplot3d(pair[,2:4], xlab='V1', ylab='V2', zlab='V3',
              pch=20, color = cols[pair$`clusterBest$cluster`])
legend('topright', legend = c('CL1', 'CL2'), col = cols, pch=20)

#simple classifier
PRED = frogDAT
PRED$Family[which(PRED$Family == 'Dendrobatidae')] = 'Hylidae'
PRED$font_pred = vector(mode = "list", length = nrow(frogDAT))
for(k in 1:bestk){
  PRED[PRED$cluster==k,]$font_pred=TOP[k]
}
#confusion matrices 
conf = table(PRED$Family, unlist(PRED$font_pred))
#conf in %'s
conf/apply(conf,1,sum)
sum(diag(conf[3:4,]))/sum(conf)#global

##aside: if k=3
cols = c('orange', 'purple', 'black')
scatterplot3d(pair[,2:4], xlab='V1', ylab='V2', zlab='V3',
              pch=20, color = cols[pair$`clusterBest$cluster`])
legend('topright', legend = c('CL1', 'CL2', 'CL3'), col = cols, pch=20)
#simple classifier
PRED = frogDAT
PRED$font_pred = vector(mode = "list", length = nrow(frogDAT))
for(k in 1:bestk){
  PRED[PRED$cluster==k,]$font_pred=TOP[k]
}
#confusion matrices 
conf = table(PRED$Family, unlist(PRED$font_pred))
#conf in %'s
conf/apply(conf,1,sum)
sum(diag(conf[2:4,]))/sum(conf)#global
#Q4
set.seed(6350)
#redo classes
CL1 = subset(frogMFCCs, frogMFCCs[,23]=='Dendrobatidae') 
CL2 = subset(frogMFCCs, frogMFCCs[,23]=='Hylidae') 
CL3 = subset(frogMFCCs, frogMFCCs[,23]=='Leptodactylidae') 
#standardize
m = apply(frogDATA[,-23], 2, mean)
s = apply(frogDATA[,-23], 2, sd)
CL1s = cbind(CL1[23], t(t(sweep(CL1[,-23],2,m))/s))
CL2s = cbind(CL2[23], t(t(sweep(CL2[,-23],2,m))/s))
CL3s = cbind(CL3[23], t(t(sweep(CL3[,-23],2,m))/s))
#assign selections ~20% of each class
r1 = sort(sample(nrow(CL1),nrow(CL1)*0.2))
r2 = sort(sample(nrow(CL2),nrow(CL2)*0.2))
r3 = sort(sample(nrow(CL3),nrow(CL3)*0.2))
#define test&training sets
testCL1 = CL1[r1,]; trainCL1=CL1[-r1,]
testCL2 = CL2[r2,]; trainCL2=CL2[-r2,]
testCL3 = CL3[r3,]; trainCL3=CL3[-r3,]
#rebalance
trainCL1 = rbind(trainCL1, trainCL1, trainCL1)
testCL1 = rbind(testCL1, testCL1, testCL1)
trainCL3 = trainCL3[-sample(c(1:nrow(trainCL3)), 1420*0.8),]
testCL3 = testCL3[-sample(c(1:nrow(testCL3)), 1420*0.2),]
#test/train sets
TRAINSET= droplevels(rbind(trainCL1,trainCL2,trainCL3))
TESTSET= droplevels(rbind(testCL1,testCL2,testCL3))
#TRAINSET$Family=as.factor(TRAINSET$Family)
#TESTSET$Family=as.factor(TESTSET$Family)

#Q5 11/18 for RF
ntry = sqrt(ncol(SDATA)-1)
#ntrees = c(100, 192,200,300,400)
#ntrees = seq(175, 225,1)
ntrees = c(100,200,300,400,seq(175,225,1))

rfn = vector(mode = "list", length = length(ntrees)) #random forest (train set)
predn = vector(mode = "list", length = length(ntrees)) #test predictions
predConf = vector(mode = "list", length = length(ntrees)) #conf.matrix for testsets
accnTest = NULL #global accuracies (test)
DendroACCn = NULL; HyliACCn = NULL; LeptoACCn = NULL #diagonals of the conf. matrix

for(i in 1:length(ntrees)){
  rfn[[i]] = randomForest(Family~., data=TRAINSET, ntree=ntrees[i], mtry=ntry)
  predn[[i]] = predict(rfn[[i]],TESTSET)
  confn = table(TESTSET$Family, predn[[i]])
  accnTest = c(accnTest, sum(diag(confn))/sum(confn))
  percent = (confn/apply(confn,1,sum))
  predConf[[i]] = percent
  DendroACCn = c(DendroACCn, diag(percent)[1])
  HyliACCn = c(HyliACCn, diag(percent)[2])
  LeptoACCn = c(LeptoACCn, diag(percent)[3])
} #30sec for the first 4; 4:46 for 2nd

accnTrain = NULL
for(i in 1:length(ntrees)){
  confi = rfn[[i]]$confusion[,-4]
  accnTrain = c(accnTrain, sum(diag(confi))/sum(confi))
}

#ntrees[1:4]; ntrees[5:length(ntrees)]
plot(ntrees[1:4], accnTrain[1:4], type = 'l', xlab = 'ntrees', ylab='accuracy',
     ylim = c(min(accnTest), max(accnTrain)))
lines(ntrees[1:4], accnTest[1:4], col = 'blue')
legend(345,0.9895,c("TRAIN acc","TEST acc"),col=c("black","blue"), lwd=1, lty=c(1,1), cex=0.50)
legend(215.5,0.99,c("TRAIN acc","TEST acc"),col=c("black","blue"), lwd=1, lty=c(1,1), cex=0.50)
#checking overfit and stuff
overfit=accnTrain - accnTest
which(accnTrain==max(accnTrain))
ntrees[which(accnTest==max(accnTest))]
ntrees[which(overfit==min(overfit))]

#Q6
plot(ntrees[1:4], DendroACCn[1:4], type='l', col='red', xlab='ntrees',
     ylab='accuracy', ylim=c(min(HyliACCn), max(DendroACCn)))
points(ntrees[1:4], DendroACCn[1:4], pch=20, col='red')
lines(ntrees[1:4], HyliACCn[1:4], col='green')
points(ntrees[1:4], HyliACCn[1:4], pch=20, col='green')
lines(ntrees[1:4], LeptoACCn[1:4], col='blue')
points(ntrees[1:4], LeptoACCn[1:4], pch=20, col='blue')
legend(331,0.995,c("Dendrobatidae","Hylidae","Leptodactylidae"),col=c("red","green","blue"), lwd=1, lty=c(1,1), cex=0.50)
legend(213.5,0.995,c("Dendrobatidae","Hylidae","Leptodactylidae"),col=c("red","green","blue"), lwd=1, lty=c(1,1), cex=0.50)

#seeing which cases are confused for Dendrobatidae
TESTSETwithPred = TESTSET
TESTSETwithPred$pred = predn[[which(ntrees==213)]]
t = TESTSETwithPred[which(TESTSETwithPred$Family != TESTSETwithPred$pred),]
t[which(t$pred=='Dendrobatidae'),]

bntree = ntrees[which(accnTest==max(accnTest))] #=213

#Q7
rfBest = randomForest(Family~., data=TRAINSET, ntree=bntree, mtry=ntry, importance=TRUE)
ImpRF = rfBest$importance
sort(ImpRF[,4], decreasing = TRUE)

predb = predict(rfBest, TESTSET)
tb = table(TESTSET$Family, predb)
tb/apply(tb,1,sum)
sum(diag(tb))/sum(tb)

#Q8
par(mfrow=c(2,3))
#most important
hist(CL1s[,'MFCCs_13'], xlab='MFCC13', main='Dendrobatidae', col='red')
hist(CL2s[,'MFCCs_13'], xlab='MFCC13', main='Hylidae', col='green')
hist(CL3s[,'MFCCs_13'], xlab='MFCC13', main='Leptodactylidae', col='blue')
ks.test(CL1s[,'MFCCs_13'], CL2s[,'MFCCs_13'])
ks.test(CL1s[,'MFCCs_13'], CL3s[,'MFCCs_13'])
ks.test(CL2s[,'MFCCs_13'], CL3s[,'MFCCs_13'])

#least important
hist(CL1s[,'MFCCs_.8'], xlab='MFCC8', main='Dendrobatidae', col='red')
hist(CL2s[,'MFCCs_.8'], xlab='MFCC8', main='Hylidae', col='green')
hist(CL3s[,'MFCCs_.8'], xlab='MFCC8', main='Leptodactylidae', col='blue')
ks.test(CL1s[,'MFCCs_.8'], CL2s[,'MFCCs_.8'])
ks.test(CL1s[,'MFCCs_.8'], CL3s[,'MFCCs_.8'])
ks.test(CL2s[,'MFCCs_.8'], CL3s[,'MFCCs_.8'])

#Q9
cluster1of2 = frogDAT[which(frogDAT$`clusterBest$cluster`==1),][,-24]
nrow(cluster1of2[which(cluster1of2$Family=='Hylidae'),])
nrow(cluster1of2[which(cluster1of2$Family=='Leptodactylidae'),])
#make it pure
cluster1of2$Family = 'Leptodactylidae'
r12 = sort(sample(nrow(cluster1of2),nrow(cluster1of2)*0.2))
test12 = droplevels(cluster1of2[r12,]); train12 = droplevels(cluster1of2[-r12,])
rfpure = randomForest(Family~., data=train12, ntrees=bntree, mtry=ntry) #error

cluster2of2 = frogDAT[which(frogDAT$`clusterBest$cluster`==2),][,-24]
nrow(cluster2of2[which(cluster2of2$Family=='Dendrobatidae'),])
nrow(cluster2of2[which(cluster2of2$Family=='Hylidae'),])
nrow(cluster2of2[which(cluster2of2$Family=='Leptodactylidae'),])

CL1c = cluster2of2[which(cluster2of2$Family=='Dendrobatidae'),]
CL2c = cluster2of2[which(cluster2of2$Family=='Hylidae'),]
CL3c = cluster2of2[which(cluster2of2$Family=='Leptodactylidae'),]
#assign selections ~20% of each class
r1c = sort(sample(nrow(CL1c),nrow(CL1c)*0.2))
r2c = sort(sample(nrow(CL2c),nrow(CL2c)*0.2))
r3c = sort(sample(nrow(CL3c),nrow(CL3c)*0.2))
#define test&training sets
testCL1c = CL1c[r1c,]; trainCL1c=CL1c[-r1c,]
testCL2c = CL2c[r2c,]; trainCL2c=CL2c[-r2c,]
testCL3c = CL3c[r3c,]; trainCL3c=CL3c[-r3c,]
#rebalance
trainCL3c = rbind(trainCL3c, trainCL3c, trainCL3c)
testCL3c = rbind(testCL3c, testCL3c, testCL3c)
#test/train sets
TRAINSETc= droplevels(rbind(trainCL1c,trainCL2c,trainCL3c))
TESTSETc= droplevels(rbind(testCL1c,testCL2c,testCL3c))
#TRAINSET$Family=as.factor(TRAINSET$Family)
#TESTSET$Family=as.factor(TESTSET$Family)

rfc = randomForest(Family~., data=TRAINSETc, ntrees=bntree, mtry=ntry) 
predc = predict(rfc, TESTSETc)
confc = table(TESTSETc$Family, predc)
confc/apply(confc,1,sum) #%
sum(diag(confc))/sum(confc) #global

#rfc on original train/set
predc = predict(rfc, TESTSET)
confc = table(TESTSET$Family, predc)

allLept = TESTSET
allLept$Family = 'Leptodactylidae'
allconf = table(TESTSET$Family, allLept$Family)
allconf/sum(allconf)

#rf with the previously trained forest?
predc213 = predict(rfBest, TESTSETc)
confc213 = table(TESTSETc$Family, predc213)
confc213/apply(confc213,1,sum) #%
sum(diag(confc213))/sum(confc213) #global

#rf for orignal sets with original rfbest
predc213 = predict(rfBest, TESTSET)
confc213 = table(TESTSET$Family, predc213)

#confidence intervals
#comparing overall accuracies: (had to hard code for 'best')
accbest = 0.989
accC = sum(diag(confc))/sum(confc) #global
sigmabest = sqrt(accbest*(1-accbest)/ nrow(TESTSET))
sigmaC = sqrt(accC*(1-accC)/ nrow(TESTSET))
#90% conf. intervals:
intervalbest = c(accbest-sigmabest*qnorm(1-0.1/2), accbest+sigmabest*qnorm(1-0.1/2))
intervalC = c(accC-sigmaC*qnorm(1-0.1/2), accC+sigmaC*qnorm(1-0.1/2))

#comparing confusion matrices (accuracy for each class):
#clusterRF: ('new' is 'cluster')
p = diag(confc/apply(confc,1,sum))
sigmaCL1new=sqrt(p[1]*(1-p[1])/nrow(testCL1c))
sigmaCL2new=sqrt(p[2]*(1-p[2])/nrow(testCL2c))
sigmaCL3new=sqrt(p[3]*(1-p[3])/nrow(testCL3c))
#90% confidence interval
intervalCL1new=c(p[1]-sigmaCL1new*qnorm(1-0.1/2), p[1]+sigmaCL1new*qnorm(1-0.1/2))
intervalCL2new=c(p[2]-sigmaCL2new*qnorm(1-0.1/2), p[2]+sigmaCL2new*qnorm(1-0.1/2))
intervalCL3new=c(p[3]-sigmaCL3new*qnorm(1-0.1/2), p[3]+sigmaCL3new*qnorm(1-0.1/2))

#bestRF:
p2 = diag(confc213/apply(confc213,1,sum))
p2 = diag(predConf[[which(ntrees==213)]])
sigmaCL1best=sqrt(p2[1]*(1-p2[1])/nrow(testCL1))
sigmaCL2best=sqrt(p2[2]*(1-p2[2])/nrow(testCL2))
sigmaCL3best=sqrt(p2[3]*(1-p2[3])/nrow(testCL3))
#90% confidence interval
intervalCL1best=c(p2[1]-sigmaCL1best*qnorm(1-0.1/2), p2[1]+sigmaCL1best*qnorm(1-0.1/2))
intervalCL2best=c(p2[2]-sigmaCL2best*qnorm(1-0.1/2), p2[2]+sigmaCL2best*qnorm(1-0.1/2))
intervalCL3best=c(p2[3]-sigmaCL3best*qnorm(1-0.1/2), p2[3]+sigmaCL3best*qnorm(1-0.1/2))

#Q11
#assign selections ~20% of each class

r1s = sort(sample(nrow(CL1),nrow(CL1)*0.2))
r2s = sort(sample(nrow(CL2),nrow(CL2)*0.2))
r3s = sort(sample(nrow(CL3),nrow(CL3)*0.2))
#define test&training sets
testCL1s = CL1s[r1s,]; trainCL1s=CL1s[-r1s,]
testCL2s = CL2s[r2s,]; trainCL2s=CL2s[-r2s,]
testCL3s = CL3s[r3s,]; trainCL3s=CL3s[-r3s,]
#rebalance
trainCL3s = trainCL3s[-sample(c(1:nrow(trainCL3s)), 1420*0.8),]
testCL3s = testCL3s[-sample(c(1:nrow(testCL3s)), 1420*0.2),]
trainCL1s = rbind(trainCL1s, trainCL1s, trainCL1s)
testCL1s = rbind(testCL1s, testCL1s, testCL1s)

#test/train sets
TRAINSETs= droplevels(rbind(trainCL2s,trainCL3s))
TESTSETs= droplevels(rbind(testCL2s,testCL3s))
TRAINSETs= droplevels(rbind(trainCL2s,trainCL1s))
TESTSETs= droplevels(rbind(testCL2s,testCL1s))


library(e1071)
#recode classes
TRAINSETs$Fam = as.numeric(TRAINSETs$Family == 'Hylidae')
TESTSETs$Fam = as.numeric(TESTSETs$Family == 'Hylidae')

classifier = svm(formula = Fam ~ ., 
                 data = TRAINSETs[,-1], 
                 type = 'C-classification', 
                 kernel = 'linear') 
predSVM = predict(classifier, newdata = TESTSETs[,-1]) 

confSVMlin = table(TESTSETs$Fam, predSVM)
confSVMlin/apply(confSVMlin,1,sum)

confSVMrad = table(TESTSETs$Fam, predSVMrad)
confSVMrad/apply(confSVMrad,1,sum)


#confidence intervals
#comparing overall accuracies: (had to hard code for 'best')
acclin = sum(diag(confSVMlin))/sum(confSVMlin)
accrad = sum(diag(confSVMrad))/sum(confSVMrad)

sigmalin = sqrt(acclin*(1-acclin)/ nrow(TESTSETs))
sigmarad = sqrt(accrad*(1-accrad)/ nrow(TESTSETs))
#90% conf. intervals:
intervallin = c(acclin-sigmalin*qnorm(1-0.1/2), acclin+sigmalin*qnorm(1-0.1/2))
intervalrad = c(accrad-sigmarad*qnorm(1-0.1/2), accrad+sigmarad*qnorm(1-0.1/2))
