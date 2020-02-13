#install.packages("synbreedData")
#install.packages("synbreed")
#install.packages("genetics")
#install.packages("SNPassoc")

library(synbreedData) #for the mice dataset
library(genetics) #for allelic counting
library(SNPassoc) #for Genome wide associaton testing
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

#Function to remove columns with no variance
zeroVar <- function(data, useNA = 'ifany') {
  out <- apply(data, 2, function(x) {length(table(x, useNA = useNA))})
  which(out==1)
}

# Load mice datasets from synbreedData
# Dataset description: https://cran.r-project.org/web/packages/synbreedData/synbreedData.pdf
data(mice)

#Structure of the datasets
str(mice)

#With these datasets, I'd like to see if we can find genetic markers that are associated with weight
#In the lack of a BMI, we'll filter covariates to get weithg comparable individuals 
#Given that weight is body weight at age of 6 weeks gs we don't need to control for birth time


#Filter the covariates (sex,caige density, etc) of the individuals that have been genotyped AND phenotyped
covarRaw<-subset(mice$covar, (phenotyped== TRUE) & (genotyped==TRUE))
#Add individuals as rownames instead of a column for easier inner joins with other datasets
covar <- covarRaw[,-1]
rownames(covar) <- covarRaw[,1]

#Merge the covariates with the phenotipe with an INNER JOIN
rawPhenoCovar<-merge (covar,mice$pheno, by="row.names", all=FALSE)

#remove cases with N/A
rawPhenoCovar<-rawPhenoCovar[complete.cases(rawPhenoCovar),]
#Add individuals as rownames instead of a column for easier inner joins with other datasets
phenoCovar <- rawPhenoCovar[,-1]
rownames(phenoCovar) <- rawPhenoCovar[,1]


# Given that the genotypic dataet is quite big:
dim(mice$geno)

# I narrowed the search to the chromosomes involve ind body weight and adiposity regulation
# Looking at the bibliography we see the candidate chromosomes which might be impactful
# Bibliography:
# https://www.ncbi.nlm.nih.gov/pubmed/9268627
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1435867/
# https://www.researchgate.net/publication/7091851_Chromosome_2_locus_Nidd5_has_a_potent_effect_on_adiposity_in_the_TSOD_mouse
  
#We go for Chr2

#Slice of the map to only conider chromosome 2  
ch2Map<-subset(mice$map, (chr == 2))

#SNP Marker names for chromosome 2
ch2MarkerNames<-c(row.names(ch2Map))

#Get the SNP data for the markers 2nd chromosome only
markerdataRaw<-mice$geno[,ch2MarkerNames]
markerdata<-markerdataRaw


# Encode for codominance (additive)
# https://en.wikipedia.org/wiki/Dominance_(genetics)#Co-dominance
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3846273/

# -1= minor frequence allele of the SNP is not present on the marker
#  0= heterozygous for minor frequency allele
#  1= homozygous for minor frequency allele

for (i in c(1:ncol(markerdataRaw))){

  #genotypes for an individual
  genotypes<-genotype(markerdataRaw[,i], sep = "/")

  #Minor allele
  idxMinor<-which.min(summary(genotypes)$allele.freq[,2])
  if (summary(genotypes)$allele.freq[idxMinor,2] != 1){
    minor<-c(names(summary(genotypes)$allele.freq[,2])[idxMinor])
  } else {
    minor<-0
  }
  #Minor allele enconding, counting the number of minor alleles present of the individual
  # -1 to get -1,0,1 
  markerdata[,i]<-(as.numeric(c(allele.count(genotypes, allele.name=minor)))-1)
    
}


dataset1<-merge (markerdata, phenoCovar, by="row.names", all=FALSE)
#remove N/A from Y
dataset2<-na.omit(dataset1)
#remove markers with no variance
fullData<-dataset2[,-zeroVar(dataset2)]

dataset1Raw<-merge (markerdataRaw,phenoCovar, by="row.names", all=FALSE)
#remove N/A from Y
dataset2Raw<-na.omit(dataset1Raw)
#remove markers with no variance
fullDataRaw<-dataset2Raw[,-zeroVar(dataset2Raw)]

#Extract weight as the response variable
Y<-fullData$weight.1

#Extract explanatory date. We'll leave only the SNP markers in the inner
#joined data with complete cases

X<-subset.data.frame(fullData, select = -c(weight.1))
X<-subset.data.frame(X, select = -c(Row.names))
X<-subset.data.frame(X, select = -c(growth.slope.1))
X<-subset.data.frame(X, select = -c(Litter))
X<-subset.data.frame(X, select = -c(CageDensity))
X<-subset.data.frame(X, select = -c(birthMonth))
X<-subset.data.frame(X, select = -c(birthYear))
X<-subset.data.frame(X, select = -c(sex))
X<-subset.data.frame(X, select = -c(CoatColour))
X<-as.data.frame(X)

#Matrix transformations needed to run glmnet
Xm<-as.matrix.data.frame(X)
Xmnum<-apply(Xm, 2, as.numeric)

#LASSO
###############

#Methodology:
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2732298

set.seed(123)

fitLasso<-cv.glmnet(x=Xmnum,y=Y,alpha=1, intercept=TRUE)
#Extract coefficients for the 1se lambda 
fitLasso.beta <- coef.glmnet(fitLasso, s=fitLasso$lambda.1se)
#Select significant coefficients only
LassoCoef <- data.frame(name = fitLasso.beta@Dimnames[[1]][fitLasso.beta@i + 1], coefficient = fitLasso.beta@x)
#These are the significant markes found by lasso regression
LassoCoef


##################
#SNP ASSOCIATION##
#################


myDat<-setupSNP(fullDataRaw,2:(ncol(fullDataRaw)-8),sep="/")

#Performs fisher test for the association of the SNPs with the response variable
#We adjust by Sex and CaigeDensity including them as covariates to the model
fitAssoc<-WGassociation(weight.1~sex+CageDensity, data=myDat, model = c("codominant"), level=0.05, quantitative = TRUE)
fitAssoc
summary(fitAssoc)

#pvalues
pvalues(fitAssoc)
plot(fitAssoc, alpha = 0.05, print.label.SNPs=FALSE)

#apply bonferroni correction to alpha
bonfCoef<-Bonferroni.sig(fitAssoc, model=c("codominant"), include.all.SNPs=TRUE,alpha=0.05)


#Inner join of the significant coefficient found by lass and the genome wide association
#to find the most significant SNPs
resultsRaw<-merge(x = bonfCoef, y = LassoCoef, by.x = "row.names", by.y="name", all = FALSE)

names(resultsRaw)[1]<-"SNP"
names(resultsRaw)[3]<-"p-value"

#Inner join with the map positions to know where these SNPS are located
results<-merge(x = resultsRaw, y = ch2Map, by.x = "SNP", by.y="row.names", all = FALSE)
results

plot(fitAssoc[c(results$SNP),], alpha = 0.05, print.label.SNPs=TRUE)


#We see a very interesting locus in 31-32 cm.
#Looking in the mice gene db filtering by genetic map position in chromosome 12:
#http://www.informatics.jax.org/marker/summary?nomen=&chromosome=2&cm=31-32.9&coordinate=&coordUnit=bp&startMarker=&endMarker=&go=&goVocab=goFunctionTerm&goVocab=goProcessTerm&goVocab=goComponentTerm&interpro=&phenotype=weight
#DISCLAMER: map position comparison between studies (genetic maps) are not
#a super reliable way of matching markers, as they don't map true physical distance: https://en.wikipedia.org/wiki/Centimorgan
#we asume that the Genetic Map referenc ind informatics.jax.org is somewhat comparable
#to the one we have in mice$map

#We find:
#A gene related with weight: Gpd2, glycerol phosphate dehydrogenase 2, mitochondrial
#A QTL (quantitative trait locus) associates with weigth:  Nidd5, non-insulin-dependent diabetes mellitus 5
#http://www.informatics.jax.org/marker/key/63419
#https://www.ncbi.nlm.nih.gov/pubmed/16688528?dopt=Abstract

#Conclusion: So we might conlcude that there is at least one significant locus in Ch2 that is associated with weight

########################
#FURTHER INVESTIGATION #
########################

#We want to run CART filtering by the significant markers + sex for adjustment
n<-c(results$SNP)
idxCart<-append (n,c("sex"))

#Run CART algorithm
fitCART<-rpart(Y~.,data=fullData[,idxCart], control=rpart.control(minsplit=2,cp=0))
fitCART.CP <- printcp(fitCART)
fitCART.CV <- fitCART.CP[,4]

#Find minimum CV error
mincverr= as.numeric(which.min(fitCART.CV))

#Find new threshold using the minumum cross validation error
s=fitCART.CP[mincverr,4]
B=1*(fitCART.CV<=s)
idx=min(which(B==1))
cp=fitCART.CP[idx,1]
#Prune tree
fitCART_pruned=prune(fitCART,cp=cp)
#Plot
rpart.plot(fitCART_pruned)

#Interestingly the CART algorithm gives us M1241 as a first node for male cases
#which is where we have most of the weight variability

#Frequency plot for males groupd by genotype on M1241
ggplot(fullData[fullData$sex=="M",], aes(weight.1, stat(density), colour = M1241)) +
  geom_freqpoly(binwidth = 1)

#Looks like carriers of 1 (minor allele homoziguos tend to have a higher weight )

#stats for males
#Statistics for males with 0 "G" alleles for M1241
summary(subset(fullData, (fullData[,c("M1241","weight.1")]$M1241==-1) & (sex=="M"))$weight.1)
#Statistics for males with 1 "G" alleles for M1241
summary(subset(fullData, (fullData[,c("M1241","weight.1")]$M1241==0) & (sex=="M"))$weight.1)
#Statistics for males with 0 "G" alleles for M1241
summary(subset(fullData, (fullData[,c("M1241","weight.1")]$M1241==1) & (sex=="M"))$weight.1)

#We look at the allele frequencies to identify the minor allele
M1241genotypes<-genotype(mice$geno[,c("M1241")], sep = "/")
summary(M1241genotypes)
plot(M1241genotypes)

#Second conclusion: We can conclude that allele G migth constitutes a risk allele in SNP M1241 
#for higher weight. Significance testing is needed.
