\name{KTSP.Classifiy}

\alias{KTSP.Classify}

\title{Function to classify samples using a KTSP classifier.}

\description{
  \code{KTSP.Classify} classifies new test samples
  using KTSP coming out of the function \code{\link{KTSP.Train}}.
  This function was used in Marchionni et al, 2013, BMC Genomics,
  and it is maintained only for backward compatibility.
  It has been replaced by \code{\link{SWAP.KTSP.Classify}}.
}

\usage{
  KTSP.Classify(data, classifier, combineFunc)
}

\arguments{
  \item{data}{the test data: a matrix in which the rows
    represent the genes and the columns the samples.}
  \item{classifier}{The output of  \code{\link{KTSP.Train}},
    a KTSP classifier.}
  \item{combineFunc}{A user defined function to combine
    the predictions of the individual K TSPs.
    If missing the consensus classification among the
    majority of the TSPs will be used.}
}


\seealso{
  \code{\link{KTSP.Train}},
  \code{\link{SWAP.KTSP.Classify}},
}


\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Luigi Marchionni \email{marchion@jhu.edu}
}


\references{
See \link{switchBox} for the references.
}


\examples{

##################################################
### Load gene expression data for the training set
data(trainingData)

### Turn into a numeric vector with values equal to 0 and 1
trainingGroupNum <- as.numeric(trainingGroup) - 1

### Show group variable for the TRAINING set
table(trainingGroupNum)


##################################################
### Train a classifier using default filtering function based on the Wilcoxon test
classifier <- KTSP.Train(matTraining, trainingGroupNum, n=8)

### Show the classifier
classifier


##################################################
### Testing on new data

### Load the example data for the TEST set
data(testingData)

### Turn into a numeric vector with values equal to 0 and 1
testingGroupNum <- as.numeric(testingGroup) - 1

### Show group variable for the TEST set
table(testingGroupNum)

### Apply the classifier to one sample of the TEST set using
### sum of votes grearter than 2
testPrediction <- KTSP.Classify(matTesting, classifier,
     combineFunc = function(x) sum(x) < 2.5)

### Show prediction
table(testPrediction, testingGroupNum)

}

\keyword{KTSP, classification, prediction}
