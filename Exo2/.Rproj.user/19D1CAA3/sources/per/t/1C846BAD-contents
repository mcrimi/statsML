#Maybe R part?
fitMaxTreeVS<-rpart(Y~.,data=as.data.frame(X), control=rpart.control(minsplit=2,cp=0, maxdepth = 5))
summary(fitMaxTreeVS)
rpart.plot(fitTree)

fitTree<-rpart(Y~M1241+M1242+M1249+M1250+M1251+M1252+M1497,data=as.data.frame(X), control=rpart.control(minsplit=2,cp=0, maxdepth = 5))
rpart.plot(fitTree)


plot(fitAssoc, alpha = 0.05, print.label.SNPs=FALSE)










#Continuous data
fullContData<-subset.data.frame(fullData, select = -c(CoatColour))
fullContData<-subset.data.frame(fullContData, select = -c(Litter))
fullContData<-subset.data.frame(fullContData, select = -c(sex))
fullContData<-subset.data.frame(fullContData, select = -c(Row.names))



#M8357       
M8306genotypes<-genotype(mice$geno[,c("M8624")], sep = "/")

#Look at M8305    

M8305genotypes<-genotype(mice$geno[,c("M8305")], sep = "/")
# Map positions
mice$map[c("M8305","M8357"),] #same position, and the regression didn't know. Encouraging
expectedGenotypes(M8305genotypes)
allele.names(M8305genotypes)
plot(M8305genotypes)
summary(M8305genotypes)
#might be rs4854342 snp? # A being the risk allele=minor

str(fullData[,"M8305"])

table(fullData[,"M8305"])

snpresponse<-fullData[,c("M8305","weight.1")]

hist(fullData$weight.1)
points(fullData[c(1:90),]$weight.1,cex = 1, col="red")
points(fullData[c(90:150),]$weight.1,cex = 1, col="blue")
hist(snpresponse)           

snpresponseDF<-as.data.frame(snpresponse)

#minor allele
hist(subset(snpresponse, M8305==1))
#No minor allele
hist(subset(snpresponse, M8305==0))  

hist(snpresponse$weight.1[snpresponse$M8305==0],freq = FALSE)
hist(snpresponse$weight.1[snpresponse$M8305==1],freq = FALSE)



#Difference between sexes
###########################

#In males
ggplot(fullData[fullData$sex=="M",], aes(weight.1, stat(density), colour = M1249 )) +
  geom_freqpoly(binwidth = 1)

#Females
#In males

ggplot(fullData[fullData$sex=="F",], aes(weight.1, stat(density), colour = M1249)) +
  geom_freqpoly(binwidth = 1)



#means males
summary(subset(fullData, (snpresponse$M1249==0) & (sex=="M"))$weight.1)
summary(subset(fullData, (snpresponse$M1249==1) & (sex=="M"))$weight.1)

#means females
summary(subset(fullData, (snpresponse$M8651==0) & (sex=="F"))$weight.1)
summary(subset(fullData, (snpresponse$M8651==1) & (sex=="F"))$weight.1)





WGstats(fitAssoc)$M8648
WGstats(fitAssoc)$M8649
qqpval(pvalues(fitAssoc))

#We see a very interesting locus in the range of 50cm to 

#M8357       
M8624genotypes<-genotype(mice$geno[,c("M8648")], sep = "/")

ggplot(fullData, aes(weight.1, stat(density), colour = M8649)) +
  geom_freqpoly()

ggplot(fullData[fullData$sex=="F",], aes(weight.1, stat(density), colour = M8648)) +
  geom_freqpoly(binwidth = 1)

ggplot(fullData[fullData$sex=="M",], aes(weight.1, stat(density), colour = M8648)) +
  geom_freqpoly(binwidth = 1)

#means males
summary(subset(fullData, (snpresponse$M8649==0) & (sex=="M"))$weight.1)
summary(subset(fullData, (snpresponse$M8651==1) & (sex=="M"))$weight.1)





WGstats(fitAssoc)$M8648




getSignificantSNPs(fitAssoc, model="codominant")

which.max(fitAssoc$`log-additive`)
fitAssoc[50]

snpresponse

dev.new(width = 1024, height = 768, unit = "px")
plot(fitAssoc, alpha=0.05, cutPval = c(0, 0.05, 1))

str(fitAssoc)
getSignificantSNPs(fitAssoc,model = c("domin"))





#association(weight.1~M8305+sex, data=myDat, model=c("codominant"), quantitative=TRUE, level-05)
#association(weight.1~M8357+sex, data=myDat, model=c("codominant"), quantitative=TRUE, level=05)
#association(weight.1~M8651+sex, data=myDat, model=c("codominant"), quantitative=TRUE, level=05)



#For proper GWAS (Control-case)


#Tesing (explore)
HWE.test(M8305genotypes)
HWE.chisq(M8305genotypes)
## Compare with individual calculations:
diseq(genotypesPop2)
diseq.ci(genotypesPop2)
HWE.chisq(genotypesPop2)
HWE.exact(genotypesPop2)


























############################
#EXTRA
#########################


#######################################
#GLM
###############

#fitGLM<-lm(Y~.,X)
#summary(fitGLM)
#idxSignificant <- summary(fitGLM)$coeff[-1,4] < 0.05
# select sig. variables
#betas <- names(idxSignificant)[idxSignificant == TRUE] 
#betas



#CART
fitMaxTreeVS<-rpart(Y~.,data=as.data.frame(X), control=rpart.control(minsplit=2,cp=0, maxdepth = 5))
summary(fitMaxTreeVS)
rpart.plot(fitMaxTreeVS)

fitMaxTreeVS.CP <- printcp(fitMaxTreeVS)
fitMaxTreeVS.CV <- fitMaxTreeVS.CP[,4]
#Find minimum CV error
mincverr= as.numeric(which.min(fitMaxTreeVS.CV))

#Find new threshold using the 1-SE rule
s=fitMaxTreeVS.CP[mincverr,4]#+fitMaxTreeVS.CP[mincverr,5]
B=1*(fitMaxTreeVS.CV<=s)
idx=min(which(B==1))
cp=fitMaxTreeVS.CP[idx,1]
fitPrunedTreeVS=prune(fitMaxTreeVS,cp=cp)

rpart.plot(fitPrunedTreeVS)

lm(Sepal.Length ~ Sepal.Width + Species,
   data=iris[iris$Species == "setosa", ])


help(synbreedData)

data(mice)

View(mice$phenoCovar)
View(mice$geno)

#trying to merge
a<-merge(mice$geno, subset(mice$map, (chr == 1)))
data<-merge(mice$phenoCovar[,1,1], mice$covar$CageDensity)



mice$phenoCovar


dim(mice$map) # remove NAs


is.na.data.frame(mice$geno)
summary(mice)
mice$info

mice$covar$id

dim(mice$geno)
mice$geno[c("A063777532"),]
head(mice$map)
dim(subset(mice$map, (chr == 1)))

c(row.names(subset(mice$map, (chr == 1))))




#genotypes for an individual
genodata<-mice$geno[c("A063777532"),c(row.names(subset(mice$map, (chr == 1))))]
#genotypes for a marker in a population

genodataPop<-mice$geno[,c("M949")]

#genotypes for an individual
genotypes<-makeGenotypes(data=genodata, sep = "/", method=as.genotype)

#genotypes for a marker in a population
genotypesPop<-makeGenotypes(data=mice$geno[,c("M249")], sep = "/", method=as.genotype)

summary(genotypes)
plot(genotypes)
summary(genotypesPop)

genotypesPop<-makeGenotypes(data=mice$geno[,c("M229")], sep = "/", method=as.genotype)
genotypesPop2<-genotype(mice$geno[,c("M229")], sep = "/")



plot(genotypesPop)
summary(genotypesPop)


fullData[which(fullData$M8305==1),]
which(fullData$M8305==1)

fullData[c(19,39,41),]$weight.1
fullData[which(fullData$M8305==1),]$weight.1



#SNP ASSOC
######################

fullData[,90:200]



data(SNPs)
# first, we create an object of class 'setupSNP'
datSNP<-setupSNP(SNPs,6:40,sep="")
# case-control study, crude analysis
association(blood.pre~snp10001, data=datSNP)

str(datSNP)
str(myDat)

data(HapMap)
head(resHapMap)
# resHapMap contains the results for a log-additive genetic model
# to get the significant SNPs for chromosome 12
getSignificantSNPs(resHapMap,chromosome=12)
# to get the significant SNPs for chromosome 5
getSignificantSNPs(resHapMap,5)
# to get the significant SNPs for chromosome X at level 1e-8
getSignificantSNPs(resHapMap,5,sig=1e-8)



example.data <- c("D/D","D/I","D/D","I/I","D/D",
                  "D/D","D/D","D/D","I/I","")
g1 <- genotype(example.data)
g1
HWE.chisq(g1)
# compare with
HWE.exact(g1)
# and
HWE.test(g1)
three.data <- c(rep("A/A",8),
                rep("C/A",20),
                rep("C/T",20),
                rep("C/C",10),
                rep("T/T",3))
g3 <- genotype(three.data)
g3
HWE.chisq(g3, B=10000)

