#dependencies
#install.packages("pls")
#install.packages("Matrix")
#install.packages("ggfortify")
#install.packages("ggplot2")
#install.packages('glmnet')
#install.packages('DAAG')
#install.packages("rpart.plot")
#install.packages("VSURF")

library(ggplot2) #for ggbiplot
library(ggfortify)
library(pls)
library(Matrix) #for matrix/df conversions
library(glmnet) #for glmnet
library(DAAG) #for cv.lm
library(MASS) #for stepAIC
library(tidyverse)
library(caret)
library(leaps)
library(rpart) #for cart
library(rpart.plot) #for cool trees
library(VSURF) #for random forest variable selection

## set the seed to make randomizations reproducible
seed=sample(5)
set.seed(seed)


#Some functions
#stepR2adj 
#This function attempts to do variable selection
#using R2Adjusted as the metric for selection
#supporting multiple directions
stepR2adj <- function (X,Y, direction = c("both", "backward", 
                                        "forward")){
  #Number of features
  p= length(X)
  if (direction=="forward"){
    #Fit model with only the first explanatory variable
    fit<- lm(Y~.,as.data.frame(X[,1]))
    
    for (j in c(2:p)){
      #Add a new explanatory variable
      stepFw <- lm(Y~.,as.data.frame(X[,c(1:j)]))
      #If the model with new variable performs worst than the previous one
      #then exit the process.
      if (summary(stepFw)$adj.r.squared<summary(fit)$adj.r.squared){
      break
      }
      #If not, store the the model and move on
      else {fit <- stepFw}
    }
    return(fit)
  } else if (direction=="backward"){
    #Fit full model
    fit <- lm(Y~.,X)
    
    for (i in c(p:1)){
      #Remove one variable from the model
      stepBack <- lm(Y~.,as.data.frame(X[,-c(i:10)]))
      #If the model with one variable less performs worst than the previous one
      #then exit the process. 
      if (summary(stepBack)$adj.r.squared<summary(fit)$adj.r.squared) {
        break
      }
      #If not, store the the model and move on 
      else {fit <- stepBack}
    }
    return(fit)
    
  } else if (direction=="both"){
    stop("Not implemented")
  } else {
    stop("Direction parameter is not valid")
  }
}

stepFisher <- function (X,Y, direction = c("both", "backward", 
                                          "forward")){
  #Number of features
  p<-length(X)
  
  #Suggested alpha values for class to prevent early
  #stop on the forward case
  alphaBW<-0.05
  alphaFW<-0.1
 
  if (direction=="forward"){
    #Fit model with only the intercept
    fit<-lm(Y~1)
    #variables to evaluate
    remainingVars<-c(1:p)
    #variables already included in the model
    includedVars<-NULL
    #Loop all the possible lenghts of the formula
    for (j in c(1:p)){
      #find the best new variable to add for a model with length j
      min_p_value<-1
      v<-1
      for (v in c(1:length(remainingVars))){
        
        stepFw <- lm(Y~.,subset(X, select=c(includedVars,remainingVars[v])))
        #get the jp-value of the jth coefficient
        p_value <- summary(stepFw)$coefficients[j+1,4]
        #browser()
        if (p_value<min_p_value){
          min_p_value<-p_value
          min_p_value_var<-remainingVars[v]
        }
      }
 
      if (min_p_value > alphaFW){
        break
      }else{
        includedVars[j]<-min_p_value_var
        remainingVars<-remainingVars[-match(min_p_value_var,remainingVars)]
      }

    }
    fit <- lm(Y~.,subset(X, select=includedVars))
    summary(fit)
    return(fit)
    
    
  } else if (direction=="backward"){
    stop("Not implemented")
  } else if (direction=="both"){
    stop("Not implemented")
  } else {
    stop("Direction parameter is not valid")
  }
}

#Root mean squared error
RMSE <- function(fitted, true){
  sqrt(mean((fitted - true)^2))
}

#Squared error
SE <- function(fitted, true){
  sum((fitted - true)^2)
}

#Mean squared error
MSE <- function(fitted, true){
  mean((fitted - true)^2)
}

#Mean absolute error
MAE <- function(fitted, true){
  mean(abs(fitted - true))
}

#R squared
R2 <- function(fitted, true){
  1 - (sum((true - fitted)^2)/sum((true - mean(true))^2))
}

#Manual cross validation for GLM
CVMSE_Custom <- function(data,seed){
  
  #MANUAL K-FOLD crossvalidation
  #Randomly shuffle the data
  randomizedData<-data[sample(nrow(data)),]
  
  #Create nfolds equally size folds
  folds <- cut(seq(1,nrow(randomizedData)),breaks=nfolds,labels=FALSE) 
  #initialize crossvalidation error to 0
  CVerrGLM <- 0
  #Perform nfolds fold cross validation
  for(i in 1:nfolds){
    #Segement your data by fold using the which() function 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- randomizedData[testIndexes, ]
    trainData <- randomizedData[-testIndexes, ]
    
    #Training data for iteration ith
    trainY<-trainData$Y
    trainX<-subset.data.frame(trainData, select = -c(Y))
    
    #Test data for iteration ith
    testY<-testData$Y
    testX<-subset.data.frame(testData, select = -c(Y))
    
    #Fit linear model with training data
    fit<-lm(trainY~.,trainX)
    #Predict using testing data
    prediction<- predict(fit, newdata=testX)
    #Sum the prediction errors
    CVerrGLM= CVerrGLM + nrow(testX)/n * MSE(fitted=prediction, true=testY)

  }
  #return the cross validation error
  return (CVerrGLM)
}


#####################
#1- DATA LOADING    #
#####################
#load dataset
dsfull <- read.delim("C:/Users/MCRIMI/Google Drive/Grad school/DSTI/Stats/Advanced/Submission/Exo1/procespin.txt")

#Extract response variables
Yraw<-dsfull$y
#Extract explanatory data
Xraw<-subset.data.frame(dsfull, select = -c(y))


#######################
#2- BASIC INSPECTION  #
#######################

#Check for missing data
is.na.data.frame(dsfull)

#No misssing data, good news!

#explore
summary(dsfull)
#looks like x1 has a lot of variance
boxplot(dsfull)
#x1 has a lot of variance compared to the other variables
#and it also seems to be in a different scale

#Scaled version of the distributions
boxplot(scale(dsfull))

#A peek at the correlations of the explenatory variables
cor(Xraw)


rk<-rankMatrix(data.matrix(Xraw))
rk[1]
#matrix is full rank

#######################################
#3- DEFINE TRAINING AND LEARNING DATA #
#######################################

#Apply logarithm to Y, as indicated in the assgignment
Y<- log(Yraw)
#Convert to data frame
X<- as.data.frame(Xraw)
#Scaled version of X
Xscaled<- scale(X)
#Bind treated data dataset
data<- cbind(X,Y)
#n has the number of observations
n= nrow(data)
#P has the number of features
p= length(X)


#######################
#3- PCA Analyisis     #
#######################

#A first naive attempt with no scaling in X

pca <- prcomp(X, center = TRUE)
summary(pca)
plot(pca)
autoplot(pca, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)

# Would seem like PC1 explains 99.27% of the dataset variance, x1 being the main contributor
# PC2 would explain 0.35% of the data, with x3, and x2 as the main influencers

#Now, let's look at the scaled version
pca.scaled <- prcomp(X, center = TRUE, scale= TRUE)
summary(pca.scaled)
plot(pca.scaled)
autoplot(pca.scaled, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)
ggbiplot(pca.scaled)

#Things look quite different in in the scaled version. It is going to be important
#to use scaling in the PCR process


##########################
# MULTIPLE LINEAR MODELS #
##########################

#Baseline generlized linear model, no variable selection
fitGLM<-lm(Y~.,X)
summary(fitGLM)
#Intertingly we see that the intercept has a lot of significance

plot(fitGLM$fitted.values,Y)
abline(0,1)
#Fit looks ok

#Testing the gaussian assumption on the error

hist(fitGLM$residuals)
#Residuals look standard normalish, with expectation 0.   

plot(fitGLM)
#Q/Q plots also demonstrates adequate fit

shapiro.test(fitGLM$residuals)
#p-value = >0.05, we cannot reject the hypothesis that 
#the residuals comes from a populatinon which noise is  normally distributed


#CROSS VALIDATION #


#Executing CV function from DAAS package

#Perform cross validation using the baseline model. 
#Three seems like a natural fit for our case where n=33 so we would have #Ln=11
nfolds=5
cvGLM<-CVlm(data=X,form.lm=fitGLM, m=nfolds, plotit="Observed", seed=seed)

#Store Overall cross validation MSE and R2adj of the model
#for comparison with other models
fitGLM.CVMSE<-attributes(cvGLM)$ms
fitGLM.R2Adj<-summary(fitGLM)$adj.r.squared

#Just for fun a manual crossvalidation function I've created.
CVMSE_Custom(data=data, seed=seed)
#Difference with attributes(cvGLM)$ms can be attrubuted to different fold randomization

paste( "Minimum CVMSE of ",  
       fitGLM.CVMSE, 
       " was produced with ", 
       paste(names(fitGLM$coefficients),collapse=' '),
       " as explanatory variables")



####################################################
##PCR (Principal component regression)             #
####################################################


#Fit PCR with all 10 components
fitPcr <-pcr(Y~., data = X, scale = TRUE, validation = "CV") #using crossvalidation and scaling
summary(fitPcr)
plot(fitPcr)
abline(0,1)

#Plot the cross validation MSE
validationplot (fitPcr, val.type="MSEP")

#Looking at the plot it seems that ncomp=9 would gives us the minimum squared error

#Replace model by fitting PCR with 9 components instead of 10
fitPcr <-pcr(Y~., data = X, scale = TRUE, #using crossvalidation, scaling
             validation = "CV", ncomp=9) 

summary(fitPcr)
plot(fitPcr)
abline(0,1)

#Plot the CVMSE of the new model with 9 components
validationplot (fitPcr, val.type="MSEP")

#Extract and store
#minumum cross validation MSE for model comparison
fitPcr.CVMSE = MSEP(fitPcr)$val[1,1,which.min(MSEP(fitPcr)$val[1,1,])]

paste( "Minimum CVMSE of ",  
       fitPcr.CVMSE, 
       " was produced with ", 
       sub(" comps","", names(which.min(MSEP(fitPcr)$val[1,1, ] ))), 
       " components")



##########################
# VARIABLE SELECTION     #
##########################

#USING R2 as Metric
#--------------------

#BACKWARDS
fitBW_R2 <-stepR2adj(X,Y,direction = "backward")
#Cross validation
cvfitBW_R2<-CVlm(data=X,form.lm=fitBW_R2, m=nfolds, plotit="Observed", seed=seed)
#Store cross validation MSE and R2adj of the model for comparison with other models
fitBW_R2.CVMSE<-attributes(cvfitBW_R2)$ms
fitBW_R2.R2Adj<-summary(fitBW_R2)$adj.r.squared

paste( "Minimum CVMSE of ",  
       fitBW_R2.CVMSE, 
       " was produced with ", 
       paste(names(fitBW_R2$coefficients),collapse=' '),
       " as explanatory variables")


#FORWARDS
fitFW_R2 <- stepR2adj(X,Y,direction = "forward")
#Cross validation
cvfitFW_R2<-CVlm(data=X,form.lm=fitFW_R2, m=nfolds, plotit="Observed", seed=seed)
#Store cross validation MSE and R2adj of the model for comparison with other models
fitFW_R2.CVMSE<-attributes(cvfitFW_R2)$ms
fitFW_R2.R2Adj<-summary(fitFW_R2)$adj.r.squared

paste( "Minimum CVMSE of ",  
       fitFW_R2.CVMSE, 
       " was produced with ", 
       paste(names(fitFW_R2$coefficients),collapse=' '),
       " as explanatory variables")


#USING Fisher as Metric
#-----------------------
fitFW_F <-stepFisher(X,Y,direction = "forward")
#Cross validation
cvfitFW_F<-CVlm(data=X,form.lm=fitFW_F, m=nfolds, plotit="Observed", seed=seed)
#Store cross validation MSE and R2adj of the model for comparison with other models
fitFW_F.CVMSE<-attributes(cvfitFW_F)$ms
fitFW_F.R2Adj<-summary(fitFW_F)$adj.r.squared

paste( "Minimum CVMSE of ",  
       fitFW_F.CVMSE, 
       " was produced with ", 
       paste(names(fitFW_F$coefficients),collapse=' '),
       " as explanatory variables")


#USING AIC as Metric
#--------------------
#FORWARD
fitFW_AIC <- stepAIC(fitGLM, direction = "forward", trace = FALSE) 
#Cross validation
cvfitFW_AIC<-CVlm(data=X,form.lm=fitFW_AIC, m=nfolds, plotit="Observed", seed=seed)
#Store Overall cross validation MSE and R2adj of the model for comparison with other models
fitFW_AIC.CVMSE<-attributes(cvfitFW_AIC)$ms
fitFW_AIC.R2Adj<-summary(fitFW_AIC)$adj.r.squared

paste( "Minimum CVMSE of ",  
       fitFW_AIC.CVMSE, 
       " was produced with ", 
       paste(names(fitFW_AIC$coefficients),collapse=' '),
       " as explanatory variables")


#BACKWARDS
fitBW_AIC <- stepAIC(fitGLM, direction = "backward", trace = FALSE) 
#Cross validation
cvfitBW_AIC<-CVlm(data=X,form.lm=fitBW_AIC, m=nfolds, plotit="Observed", seed=seed)
#Store Overall cross validation MSE and R2adj of the model for comparison with other models
fitBW_AIC.CVMSE<-attributes(cvfitBW_AIC)$ms
fitBW_AIC.R2Adj<-summary(fitBW_AIC)$adj.r.squared

paste( "Minimum CVMSE of ",  
       fitBW_AIC.CVMSE, 
       " was produced with ", 
       paste(names(fitBW_AIC$coefficients),collapse=' '),
       " as explanatory variables")


#BOTH
fitBoth_AIC <- stepAIC(fitGLM, direction = "both", trace = FALSE) 
#Cross validation
cvfitBoth_AIC<-CVlm(data=X,form.lm=fitBoth_AIC, m=nfolds, plotit="Observed", seed=seed)
#Store Overall cross validation MSE and R2adj of the model for comparison with other models
fitBoth_AIC.CVMSE<-attributes(cvfitBoth_AIC)$ms
fitBoth_AIC.R2Adj<-summary(fitBoth_AIC)$adj.r.squared

paste( "Minimum CVMSE of ",  
       fitBoth_AIC.CVMSE, 
       " was produced with ", 
       paste(names(fitBoth_AIC$coefficients),collapse=' '),
       " as explanatory variables")


#####################################
#LASSO/RIDGE VARIABLE SELECTION     #
#####################################

#LASSO PENALIZED REGRESSION

#Run cross validation. Using intercept as it seems significant to the model.
#Given the intercept significance reported in the GLSM, we include in the model
#alpha=1 is Lasso
set.seed(seed)
fitLasso<-cv.glmnet(x=as.matrix(X),y=Y,alpha=1, intercept=TRUE, nfolds=nfolds)

#We see the CV MSE plotted against the log(lambda)
plot(fitLasso)
#See the values marked in the graphed
log(fitLasso$lambda.1se)
log(fitLasso$lambda.min)
#We can see the lambda that gives the min cross validation errr and th
#one form 1sd rule

#Extract coefficients for the 1se lambda 
fitLasso.beta <- coef.glmnet(fitLasso, s=fitLasso$lambda.1se)
#Extract associated CV MSE
fitLasso.CVMSE<-fitLasso$cvm[which(fitLasso$lambda == fitLasso$lambda.1se)]

paste( "Minimum CVMSE of ",  
       fitLasso.CVMSE, 
       " was produced with ", 
       paste(fitLasso.beta,collapse=' '),
       " as explanatory variables")


#looking at the coefficients/lambda progression
plot(glmnet(x=as.matrix(X),y=Y,alpha=1, intercept=TRUE),label=TRUE, "lambda")
#plot 1se lambda
abline(v = log(fitLasso$lambda.1se), col="red", lwd=1, lty=2)

############################
#RIDGE PENALIZED REGRESSION#
############################

#Run Cross validation
#Given the intercept significance reported in the GLSM, we include in the model
#alpha=0 is Ridge
set.seed(seed)
fitRidge<-cv.glmnet(x=as.matrix(X),y=Y,alpha=0, intercept=TRUE, nfolds=nfolds)

#We see the CV MSE plotted against the log(lambda)
plot(fitRidge)
#See the values marked in the graphed
log(fitRidge$lambda.1se)
log(fitRidge$lambda.min)
#Extract coefficients
fitRidge.beta <- coef.glmnet(fitRidge, s=fitRidge$lambda.1se)
#Extract CV MSE
fitRidge.CVMSE<-fitRidge$cvm[fitRidge.CVMSE<-fitRidge$cvm[which(fitRidge$lambda == fitRidge$lambda.1se)]]

#Show coefficients and CVMSE
fitRidge.beta
fitRidge.CVMSE

#looking at the coefficients/lambda progression
plot(glmnet(x=as.matrix(X),y=Y,alpha=0, intercept=TRUE),label=TRUE, "lambda")
#plot 1se lambda
abline(v = log(fitRidge$lambda.1se), col="red", lwd=1, lty=2)

#Seems like Ridge doesn't provide any feature selection at all as it 
#returns all coefficients.
#It remains to be seen why the CVMSE is so big compared to the linear
#model with all variables



########################################################################
#EXPERIMENTAL STABILIZATION OF RIDGE AND LASSO PENALTY REGRESSIONS    #
########################################################################
#Given that folds are created randomly, the suggested models can substantially change
#Im trying to stabilize by doing 100 cross validations and taking CV error that is
#closer to the overall mean and then taking the corresponding lambda

#Matrix to store the 100 minimum CV mean square error the lambdas associated to them
MSEs_Lambdas_Lasso <- NULL

for (i in 1:100){
  cvLassoMean <- cv.glmnet(x=as.matrix(X),y=Y,alpha=1, intercept=TRUE, nfolds=nfolds)  
  #Taking the minimum lambda this time, hoping that the stabilization provides a similar
  #effect than taking the 1se one
  lambda.cvLassoMean=cvLassoMean$lambda.min  #the best value of lambda among the grid with CV and 1 SE rule
  #Append minimum CV mean squared erre
  idx<-which(cvLassoMean$lambda == cvLassoMean$lambda.min)
  MSEs_Lambdas_Lasso  <- rbind(MSEs_Lambdas_Lasso, 
                         c(cvLassoMean$cvm[idx],
                         cvLassoMean$lambda.min))
}
#Errors are normally distributed, mean is a good expectation E
hist(MSEs_Lambdas_Lasso[,1])

#Take the lambda of the CV MSE that is closer to the mean
tunedLambdaLasso<- MSEs_Lambdas_Lasso[which.min(abs(mean(MSEs_Lambdas_Lasso[,1])-MSEs_Lambdas_Lasso[,1])),2]

#Create regression with the tuned lambda
fitLassoMean=glmnet(x=as.matrix(X),y=Y,intercept=TRUE, lambda=tunedLambdaLasso,alpha=1)

#Take the CV MSE  that is closer to the mean
fitLassoMean.CVMSE<- MSEs_Lambdas_Lasso[which.min(abs(mean(MSEs_Lambdas_Lasso[,1])-MSEs_Lambdas_Lasso[,1])),1]

fitLassoMean.beta<-coef(fitLassoMean, s=tunedLambdaLasso)


#Look  at the (scaled coefficients vs lambda progression to see where tuned lambda falls )
plot(glmnet(scale(as.matrix(X)),as.matrix(Y),alpha=1),label=TRUE, "lambda")
abline(v = log(tunedLambdaLasso), col="red", lwd=1, lty=2)
#Show coefficients
fitLassoMean.beta



############################
#CART AND RANDOM FOREST    #
############################


#Construct maximal tree
set.seed(seed)
fitMaxTreeFull=rpart(Y~.,data=X, control=rpart.control(minsplit=2,cp=0))
plot(predict(fitMaxTreeFull),Y)
abline(0,1)


#Use Random Forest for variable selection
#VSURF
fitVSURF<-VSURF(x=X,y=Y,mtry = 10)
plot(fitVSURF)
fitVSURF$varselect.pred


#Construct maximal tree using variable selection coming from VSURF
set.seed(seed)
fitMaxTreeVS=rpart(Y~.,data=X[c(fitVSURF$varselect.pred)], control=rpart.control(minsplit=2,cp=0))
plot(predict(fitMaxTreeVS),Y)
abline(0,1)

#Prune Maximal Tree with all the variables
fitMaxTreeFull.CP <- printcp(fitMaxTreeFull)
fitMaxTreeFull.CV <- fitMaxTreeFull.CP[,4]
mincverr= as.numeric(which.min(fitMaxTreeFull.CV))
s=fitMaxTreeFull.CP[mincverr,4]+fitMaxTreeFull.CP[mincverr,5]
s=min(s)
B=1*(fitMaxTreeFull.CV<=s)
idx=min(which(B==1))
cp=fitMaxTreeFull.CP[idx,1]
fitPrunedTreeFull=prune(fitMaxTreeFull,cp=cp)
# 1SE cross validation mean square error for model comparison
fitPrunedTreeFull.CVMSE= fitMaxTreeFull.CP[idx,4]
#Summary and plot
summary(fitPrunedTreeFull)
rpart.plot(fitPrunedTreeFull)

#Prune Maximal Tree with only the selected variables using RF
fitMaxTreeVS.CP <- printcp(fitMaxTreeVS)
fitMaxTreeVS.CV <- fitMaxTreeVS.CP[,4]
#Find minimum CV error
mincverr= as.numeric(which.min(fitMaxTreeVS.CV))

#Find new threshold using the 1-SE rule
s=fitMaxTreeVS.CP[mincverr,4]+fitMaxTreeVS.CP[mincverr,5]
B=1*(fitMaxTreeVS.CV<=s)
idx=min(which(B==1))
cp=fitMaxTreeVS.CP[idx,1]
fitPrunedTreeVS=prune(fitMaxTreeVS,cp=cp)

# 1SE cross validation mean square error for model comparison
fitPrunedTreeVS.CVMSE= fitMaxTreeVS.CP[idx,4]

#Summary and plot
summary(fitPrunedTreeVS)
rpart.plot(fitPrunedTreeVS)


#############
# CONCLUSION#
#############

#Model comparisson

#Create a list with the models

models<-list("name"=deparse(substitute(fitGLM)), "CVMSE"=fitGLM.CVMSE)
models$name[[2]]<-deparse(substitute(fitPcr))
models$CVMSE[[2]]<-fitGLM.CVMSE
models$name[[3]]<-deparse(substitute(fitBW_AIC))
models$CVMSE[[3]]<-fitBW_AIC.CVMSE
models$name[[4]]<-deparse(substitute(fitFW_AIC))
models$CVMSE[[4]]<-fitFW_AIC.CVMSE
models$name[[5]]<-deparse(substitute(fitBoth_AIC))
models$CVMSE[[5]]<-fitBoth_AIC.CVMSE
models$name[[6]]<-deparse(substitute(fitFW_R2))
models$CVMSE[[6]]<-fitFW_R2.CVMSE
models$name[[7]]<-deparse(substitute(fitBW_R2))
models$CVMSE[[7]]<-fitBW_R2.CVMSE
models$name[[8]]<-deparse(substitute(fitFW_F))
models$CVMSE[[8]]<-fitFW_F.CVMSE
models$name[[9]]<-deparse(substitute(fitRidge))
models$CVMSE[[9]]<-fitRidge.CVMSE
models$name[[10]]<-deparse(substitute(fitLasso))
models$CVMSE[[10]]<-fitLasso.CVMSE
models$name[[11]]<-deparse(substitute(fitLassoMean))
models$CVMSE[[11]]<-fitLassoMean.CVMSE
models$name[[12]]<-deparse(substitute(fitPrunedTreeVS))
models$CVMSE[[12]]<-fitPrunedTreeVS.CVMSE
models$name[[13]]<-deparse(substitute(fitPrunedTreeFull))
models$CVMSE[[13]]<-fitPrunedTreeFull.CVMSE


#Print best performing model
paste( "Minimum CVMSE of ",  
       models$CVMSE[which.min(models$CVMSE)], 
       " was produced by the model ", 
       models$name[which.min(models$CVMSE)], 
       ".")

############
#CONCLUSION#
############
#Using this particular seed, the best performing model seems to be X, with model 
#Y following close with a similar performance and even outperforming X in some runs 
#with different seeds. 
#A closer comparison between the two models show no common denominator in the selected 
#explanatory variables, which is very intriguing. 
#Being forced to conclude I'd say, that based on the extreme simplicity of the model created 
#by random forest variable selection posterior CART tree pruning (single node, x8 as a single explanatory variable) 
#we would incline to choose model X, but it remains to be seen why the 
#intercept x1,x2,x4 and x5 perform so well in the linear context while x8 seems to be irrelevant.
AICWins<-o
CARTWins<-0
nstab=100
for (i in 1:nstab){
  
  fitMaxTreeVS=rpart(Y~.,data=X[c(fitVSURF$varselect.pred)], control=rpart.control(minsplit=2,cp=0))
  
  nfolds=roof(runif(1,2,5))
  set.seed(runif(1,-10000,10000))
  if (fitBoth_AIC.CVMSE<fitPrunedTreeVS.CVMSE){
    AICWins=AICWins+1
  }else{
    CARTWins=CARTWins+1 
  }
}

