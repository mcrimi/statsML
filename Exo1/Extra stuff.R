#Create a training and testing data partition

# README: In this particular excercuse I chose to NOT USE training and testing 
# partitions and use cross validation instead considering the small amount 
# of observations (n=33)


## 75% of the sample size goes for training
trn_size <- floor(0.75 * nrow(data))
#randomly pick the training individuals
train_ind <- sample(seq_len(nrow(data)), size = trn_size)

#create training and testing subdatasets
trainData <- data[train_ind, ]
testData <- data[-train_ind, ]

#Training dataset
Ytrn<-trainData$Y
Xtrn<-subset.data.frame(trainData, select = -c(Y))

#Testing dataset
Ytst<-testData$Y
Xtst<-subset.data.frame(testData, select = -c(Y))


X<-Xtrn
Y<-Ytrn



#Aadd intercept
Xint<-X
Xint['Int']= as.numeric(1)


#Just for me:

#Compare test performance of tree vs both


fitBoth_AIC


YpredAIC<-predict(fitBoth_AIC, newdata = Xtst)
YpredTree<-predict(fitPrunedTreeVS, newdata = Xtst)
YpredGLM<-predict(fitGLM, newdata = Xtst)

RMSE(fitted=YpredAIC, true=Ytst)
RMSE(fitted=YpredTree, true=Ytst)
RMSE(fitted=YpredGLM, true=Ytst)

plot(YpredGLM,Ytst)
abline(0,1)

plot(YpredAIC,Ytst)
abline(0,1)

plot(YpredTree,Ytst)
abline(0,1)

summary(fitBoth_AIC)
boxplot(Xscaled)

gitGLM<-lm(Y~x1+x2+x4+x5+x8,X)



##################################
#VALIDATION USING VALIDATION SET #
##################################
#TESITNG THE FORWARD VARIABLE SELECTION MODEL
fitFWTest<- predict(fitFW, newdata=Xtst)
fitFWTest.R2<-R2(fitted= fitFWTest, true=Ytst)
fitFWTest.RMSE<-RMSE(fitted= fitFWTest, true=Ytst)
fitFWTest.MSE<-MSE(fitted= fitFWTest, true=Ytst)

#TESITNG THE BACKWARDS VARIABLE SELECTION MODEL
fitBWTest<- predict(fitBW, newdata=Xtst)
fitBWTest.R2=R2(fitted= fitBWTest, true=Ytst)
fitBWTest.RMSE<-RMSE(fitted= fitBWTest, true=Ytst)
fitBWTest.MSE<-MSE(fitted= fitBWTest, true=Ytst)

#TESITNG THE FORWARDBACKWARDS VARIABLE SELECTION MODEL
fitFWBWTest<- predict(fitFWBW, newdata=Xtst)
fitFWBWTest.R2=R2(fitted= fitFWBWTest, true=Ytst)
fitFWBWTest.RMSE<-RMSE(fitted= fitFWBWTest, true=Ytst)
fitFWBWTest.MSE<-MSE(fitted= fitFWBWTest, true=Ytst)



#STABILIZATION FOR RIDGE

#Matrix containing the iterations 1se lambda and corresponding cross validation error
MSEs_Lambdas_Ridge <- NULL

for (i in 1:100){
  cvRidgeMean <- cv.glmnet(x=as.matrix(Xint),y=Y,alpha=0, intercept=TRUE, nfolds=nfolds)  
  lambda.cvRidgeMean=cvRidgeMean$lambda.1se  #the best value of lambda among the grid with CV and 1 SE rule
  i <- which(cvRidgeMean$lambda == cvRidgeMean$lambda.1se)
  MSEs_Lambdas_Ridge  <- rbind(MSEs_Lambdas_Ridge, c(cvRidgeMean$cvm[i], cvRidgeMean$lambda.1se))
}

#Take the cross validation error that is closer to the mean of the iterations
fitRidgeMean.CVMSE<- MSEs_Lambdas_Ridge[which.min(abs(mean(MSEs_Lambdas_Ridge[,1])-MSEs_Lambdas_Ridge[,1])),1]
#Take the lambda corresponding to that cross validation
tunedLambdaRidge<- MSEs_Lambdas_Ridge[which.min(abs(mean(MSEs_Lambdas_Ridge[,1])-MSEs_Lambdas_Ridge[,1])),2]

fitRidgeMean=glmnet(x=as.matrix(Xint),y=Y,intercept=TRUE, lambda=tunedLambdaRidge,alpha=0)
fitRidgeMean.R2<-fitRidgeMean$dev.ratio


#DELETE
getlasso<- function(X,Y) #Y the respose, X is the matrix with the explanatory variables
{
  #p=dim(X)[2] #number of explanatory variables
  lambda.cv=cv.glmnet(X,Y,intercept=FALSE)$lambda.1se  #the best value of lambda among the grid with CV and 1 SE rule
  bh=glmnet(X,Y,intercept = FALSE,lambda=lambda.cv)$beta
  getlasso=bh
}

getridge<- function(X,Y) #Y the respose, X is the matrix with the explanatory variables
{
  #p=dim(X)[2] #number of explanatory variables
  lambda.cv=cv.glmnet(X,Y,alpha=0,intercept=FALSE)$lambda.1se  #the best value of lambda among the grid with CV and 1 SE rule
  bh=glmnet(X,Y,intercept = FALSE,alpha=0,lambda=lambda.cv)$beta
  getridge=bh
}
#/DELETE


generalizedStep<- function(X,Y, direction, stopping){
  
  if (direction="forward") { 
    if (stopping="R2") { 
      expr1
    } else if (stopping="F") {
      expr2
    } else if  (stopping="AIAC") {
      expr3
    } else {
      print("Unrecognized stopping parameter. Valid values are 'R2', 'F' and 'AIAC'")
    }
  } else if (direction="backward") {
    fitBW <- lm(Y~.,X)
    if (stopping="R2") { 
      expr1
    } else if (stopping="F") {
      expr2
    } else if  (stopping="AIAC") {
      expr3
    } else {
      print("Unrecognized stopping parameter. Valid values are 'R2', 'F' and 'AIAC'")
    }
  } else if  (direction="both") {
    if (stopping="R2") { 
      expr1
    } else if (stopping="F") {
      expr2
    } else if  (stopping="AIAC") {
      expr3
    } else {
      print("Unrecognized stopping parameter. Valid values are 'R2', 'F' and 'AIAC'")
    }
  } else {
    print("Unrecognized direction parameter. Valid values are 'forward', 'backward' and 'both'")
  }
  
}







#manual backwards strategy using R2 adjusted as the metric





cvfitBW<-CVlm(data=X,form.lm=fitBW, m=nfolds, plotit="Observed", seed=132)

#Store Overall cross validation MSE and R2adj of the model
#for comparison with other models
fitBW.CVMSE<-attributes(cvfitBW)$ms
fitBW.RdAdj<-summary(fitBW)$adj.r.squared


#Get better results if I remove x10, presumed to be noisy


#MANUAL FORWARD using R2 adjusted as the selection metric




cvfitFW<-CVlm(data=X,form.lm=fitFW, m=nfolds, plotit="Observed", seed=132)

#Store Overall cross validation MSE and R2adj of the model
#for comparison with other models
fitFW.CVMSE<-attributes(cvfitFW)$ms
fitFW.RdAdj<-summary(fitFW)$adj.r.squared




#Both directions, using a library



# Stepwise regression model, using Akaike information criterion

fitFWBW <- stepAIC(fitGLM, direction = "both", trace = FALSE) 

#Cross validation
cvfitFWBW<-CVlm(data=X,form.lm=fitFWBW, m=nfolds, plotit="Observed", seed=132)

#Store Overall cross validation MSE and R2adj of the model
#for comparison with other models
fitFWBW.CVMSE<-attributes(cvfitFWBW)$ms
fitFWBW.RdAdj<-summary(fitFWBW)$adj.r.squared

summary(fitFWBW)
summary(fitFWBW)$adj.r.squared


##########################
#MANUAL PCR AND TESTING  #
##########################

#LM using all of the components
projectedTrain <- as.data.frame(predict(pca.scaled,X),
                                stringsAsFactors = FALSE)

projectedTrain$y <- Y
vars = colnames(projectedTrain)
varexpr = paste(vars, collapse="+")
fmla = paste("y ~", varexpr)

pcrManualAll <- lm(fmla,data=projectedTrain)
summary(pcrManualAll)



#PCR VALIDATION USING VALIDATION SET                      #
pcr.testPred <- predict(pcr, newdata=Xtst, scale=TRUE)
maxR2=0

# loop over the fitted values by number of components used
for (i in c(10:1))
{
  currentR2<-R2(fitted= pcr.testPred[,"Y",i], true=Ytst)
  if (currentR2>maxR2) {
    maxR2<-currentR2
    ncomp<-i
  }
  
}

pcr.testPredR2 <-maxR2
pcr.testPredR2

#TODO:ncomp=7 gives us the lowest MSE so far, it remains to be seen why.
#This is aligned with the CV and adjCV error we got on 9 comps via summary(pcr)






