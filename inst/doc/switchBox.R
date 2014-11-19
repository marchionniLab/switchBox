### R code from vignette source 'switchBox.Rnw'

###################################################
### code chunk number 1: start
###################################################
options(width=85)
options(continue=" ")
rm(list=ls())


###################################################
### code chunk number 2: switchBox.Rnw:181-183 (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("switchBox")


###################################################
### code chunk number 3: switchBox.Rnw:187-188
###################################################
require(switchBox)


###################################################
### code chunk number 4: trainData
###################################################
### Load the example data for the TRAINING set
data(trainingData)


###################################################
### code chunk number 5: switchBox.Rnw:214-217
###################################################
class(matTraining)
dim(matTraining)
str(matTraining)


###################################################
### code chunk number 6: switchBox.Rnw:223-225
###################################################
### Show group variable for the TRAINING set
table(trainingGroup)


###################################################
### code chunk number 7: testData
###################################################
### Load the example data for the TEST set
data(testingData)


###################################################
### code chunk number 8: switchBox.Rnw:247-250
###################################################
class(matTesting)
dim(matTesting)
str(matTesting)


###################################################
### code chunk number 9: switchBox.Rnw:256-258
###################################################
### Show group variable for the TEST set
table(testingGroup)


###################################################
### code chunk number 10: switchBox.Rnw:281-289
###################################################
### The arguments to the "SWAP.KTSP.Train" function
args(SWAP.KTSP.Train)
### Train a classifier using default filtering function based on the Wilcoxon test
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup, krange=c(3:15))
### Show the classifier
classifier
### Extract the TSP from the classifier
classifier$TSPs


###################################################
### code chunk number 11: switchBox.Rnw:300-305
###################################################
### The arguments to the "SWAP.KTSP.Train" function
args(SWAP.Filter.Wilcoxon)
### Retrieve the top best 4 genes using default Wilcoxon filtering
### Note that there are ties
SWAP.Filter.Wilcoxon(trainingGroup, matTraining, featureNo=4)


###################################################
### code chunk number 12: switchBox.Rnw:312-318
###################################################
### Train a classifier from the top 4 best genes 
### according to Wilcoxon filtering function
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      FilterFunc=SWAP.Filter.Wilcoxon, featureNo=4)
### Show the classifier
classifier


###################################################
### code chunk number 13: switchBox.Rnw:324-328
###################################################
### To use all features "FilterFunc" must be set to NULL
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup, FilterFunc=NULL)
### Show the classifier
classifier


###################################################
### code chunk number 14: switchBox.Rnw:347-350
###################################################
### An alternative filtering function selecting 20 random features
random10 <- function(situation, data) { sample(rownames(data), 10) }
random10(trainingGroup, matTraining)


###################################################
### code chunk number 15: switchBox.Rnw:358-365
###################################################
### An alternative filtering function based on a t-test
topRttest <- function(situation, data, quant = 0.75) {
	out <- apply(data, 1, function(x, ...) t.test(x ~ situation)$statistic )
	names(out[ abs(out) > quantile(abs(out), quant) ])
}
### Show the top 5% features using the newly defined filtering function
topRttest(trainingGroup, matTraining, quant=0.95)


###################################################
### code chunk number 16: switchBox.Rnw:373-378
###################################################
### Train with t-test and krange
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      FilterFunc = topRttest, quant = 0.9, krange=c(15:30) )
### Show the classifier
classifier


###################################################
### code chunk number 17: switchBox.Rnw:401-405
###################################################
set.seed(123)
somePairs <- matrix(sample(rownames(matTraining), 6^2, replace=FALSE), ncol=2)
head(somePairs)
dim(somePairs)


###################################################
### code chunk number 18: switchBox.Rnw:412-417
###################################################
### Train
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      RestrictedPairs = somePairs, krange=3:16)
### Show the classifier
classifier


###################################################
### code chunk number 19: switchBox.Rnw:425-432
###################################################
### Train
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      RestrictedPairs = somePairs,
			      FilterFunc = topRttest, quant = 0.3,
			      krange=c(3:10) )
### Show the classifier
classifier


###################################################
### code chunk number 20: switchBox.Rnw:455-468
###################################################
### Train a classifier
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      FilterFunc = NULL, krange=8)
### Compute the statistics using the default parameters:
### counting the signed TSP votes
ktspStatDefault <- SWAP.KTSP.Statistics(inputMat = matTraining,
    classifier = classifier)
### Show the components in the output
names(ktspStatDefault)
### Show some of the votes
head(ktspStatDefault$comparisons[ , 1:2])
### Show default statistics
head(ktspStatDefault$statistics)


###################################################
### code chunk number 21: switchBox.Rnw:474-479
###################################################
### Compute
ktspStatSum <- SWAP.KTSP.Statistics(inputMat = matTraining,
    classifier = classifier, CombineFunc=sum)
### Show statistics obtained using the sum
head(ktspStatSum$statistics)


###################################################
### code chunk number 22: switchBox.Rnw:485-490
###################################################
### Compute
ktspStatThreshold <- SWAP.KTSP.Statistics(inputMat = matTraining,
    classifier = classifier,  CombineFunc = function(x) sum(x) > 2 )
### Show statistics obtained using the threshold
head(ktspStatThreshold$statistics)


###################################################
### code chunk number 23: switchBox.Rnw:498-503 (eval = FALSE)
###################################################
## ### Make a heatmap showing the individual TSPs votes
## colorForRows <- as.character(1+as.numeric(trainingGroup))
## heatmap(1*ktspStatThreshold$comparisons, scale="none",
##     margins = c(10, 5), cexCol=0.5, cexRow=0.5,
##     labRow=trainingGroup, RowSideColors=colorForRows)


###################################################
### code chunk number 24: fig1
###################################################
### Make a heatmap showing the individual TSPs votes
colorForRows <- as.character(1+as.numeric(trainingGroup))
heatmap(1*ktspStatThreshold$comparisons, scale="none",
    margins = c(10, 5), cexCol=0.85, cexRow=1,
    labRow=trainingGroup, RowSideColors=colorForRows)


###################################################
### code chunk number 25: switchBox.Rnw:538-546
###################################################
### Show the classifier
classifier
### Apply the classifier to the TRAINING set
trainingPrediction <- SWAP.KTSP.Classify(matTraining, classifier)
### Show
str(trainingPrediction)
### Resubstitution performance in the TRAINING set
table(trainingPrediction, trainingGroup)


###################################################
### code chunk number 26: switchBox.Rnw:558-565
###################################################
### Usr a CombineFunc based on  sum(x) > 5.5
trainingPrediction <- SWAP.KTSP.Classify(matTraining, classifier,
					 DecisionFunc = function(x) sum(x) > 5.5 )
### Show
str(trainingPrediction)
### Resubstitution performance in the TRAINING set
table(trainingPrediction, trainingGroup)


###################################################
### code chunk number 27: switchBox.Rnw:574-578
###################################################
### Classify one sample
testPrediction <- SWAP.KTSP.Classify(matTesting[ , 1, drop=FALSE], classifier)
### Show
testPrediction


###################################################
### code chunk number 28: switchBox.Rnw:586-592
###################################################
### Apply the classifier to the complete TEST set
testPrediction <- SWAP.KTSP.Classify(matTesting, classifier)
### Show
table(testPrediction)
### Resubstitution performance in the TEST set
table(testPrediction, testingGroup)


###################################################
### code chunk number 29: switchBox.Rnw:601-606
###################################################
### APlly the classifier using sum(x)  > 5.5
testPrediction <- SWAP.KTSP.Classify(matTesting, classifier,
				     DecisionFunc = function(x) sum(x) > 5.5 )
### Resubstitution performance in the TEST set
table(testPrediction, testingGroup)


###################################################
### code chunk number 30: switchBox.Rnw:623-628
###################################################
### Compute the scores using all features for all possible pairs
scores <- SWAP.CalculateSignedScore(matTraining,  trainingGroup, FilterFunc=NULL)
### Show scores
class(scores)
dim(scores$score)


###################################################
### code chunk number 31: switchBox.Rnw:635-639
###################################################
### Get the scores
scoresOfInterest <- diag(scores$score[ classifier$TSPs[,1] , classifier$TSPs[,2] ])
### Their absolute value should corresponf to the scores returned by SWAP.KTSP.Train
all(classifier$score == abs(scoresOfInterest))


###################################################
### code chunk number 32: switchBox.Rnw:649-660
###################################################
### Compute the scores with default filtering function
scores <- SWAP.CalculateSignedScore(matTraining, trainingGroup, featureNo=20 )
### Show scores
dim(scores$score)
### Compute the scores without the default filtering function
### and using restricted pairs
scores <- SWAP.CalculateSignedScore(matTraining, trainingGroup,
				    FilterFunc = NULL, RestrictedPairs = somePairs )
### Show scores
class(scores$score)
length(scores$score)


###################################################
### code chunk number 33: switchBox.Rnw:667-668 (eval = FALSE)
###################################################
## hist(scores$score, col="salmon", main="TSP scores")


###################################################
### code chunk number 34: fig2
###################################################
hist(scores$score, col="salmon", main="TSP scores")


###################################################
### code chunk number 35: switchBox.Rnw:701-708
###################################################
### Phenotypic group variable for the 78 samples
table(trainingGroup)
levels(trainingGroup)
### Turn into a numeric vector with values equal to 0 and 1
trainingGroupNum <- as.numeric(trainingGroup) - 1
### Show group variable for the TRAINING set
table(trainingGroupNum)


###################################################
### code chunk number 36: switchBox.Rnw:713-717
###################################################
### Train a classifier using default filtering function based on the Wilcoxon test
classifier <- KTSP.Train(matTraining, trainingGroupNum, n=8)
### Show the classifier
classifier


###################################################
### code chunk number 37: switchBox.Rnw:723-729
###################################################
### Apply the classifier to one sample of the TEST set using
### sum of votes less  than 2.5
trainPrediction <- KTSP.Classify(matTraining, classifier,
				 combineFunc = function(x) sum(x) < 2.5)
### Contingency table
table(trainPrediction, trainingGroupNum)


###################################################
### code chunk number 38: switchBox.Rnw:737-744
###################################################
### Phenotypic group variable for the 307 samples
table(testingGroup)
levels(testingGroup)
### Turn into a numeric vector with values equal to 0 and 1
testingGroupNum <- as.numeric(testingGroup) - 1
### Show group variable for the TEST set
table(testingGroupNum)


###################################################
### code chunk number 39: switchBox.Rnw:751-759
###################################################
### Apply the classifier to one sample of the TEST set using
### sum of votes less than 2.5
testPrediction <- KTSP.Classify(matTesting, classifier,
     combineFunc = function(x) sum(x) < 2.5)
### Show prediction
table(testPrediction)
### Contingency table
table(testPrediction, testingGroupNum)


###################################################
### code chunk number 40: sessioInfo
###################################################
toLatex(sessionInfo())


