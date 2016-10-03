\name{SWAP.GetKTSP.TrainTestResults}

\alias{SWAP.GetKTSP.TrainTestResults}

\title{Trains a kTSP on given training data and provides
performance on testing data.}

\description{
  Trains a kTSP on given training data and provides 
  \code{getkTSPResult} output for both training and
  testing data.
}

\usage{
SWAP.GetKTSP.TrainTestResults(trainMat, trainGroup, testMat, 
  testGroup, classes=NULL, predictions=FALSE, 
  decision_values=FALSE, ...)
}

\arguments{
  \item{trainMat}{is a matrix of data for training with rows being the
  features (such as gene names, if the matrix if gene expression
  data) and columns being the samples.}
  \item{trainGroup}{is a factor or a vector providing the phenotype class
  each training sample belongs to. It should correspond to the order of samples
  given by the columns of \code{trainMat}.}
  \item{testMat}{is a matrix of data for testing with rows being the
  features (such as gene names, if the matrix if gene expression
  data) and columns being the samples.}
  \item{testGroup}{is a factor or a vector providing the phenotype class
  each testing sample belongs to. It should correspond to the order of samples
  given by the columns of \code{testMat}.}
  \item{classes}{is a vetor of length 2 providing the two phenotype or
  class labels.}
  \item{predictions}{is a logical indicating whether to return the 
  predictions or not.}
  \item{decision_values}{is a logical indicating whether to return the 
  decision values or not.}
  \item{...}{any further arguments to be passed on for k-TSP training.}
}

\value{
  A list with items:
  \item{classifier}{The trained k-TSP classifier.}
  \item{train}{Training performance.}
  \item{train}{Testing performance.}
}

\seealso{
  \code{\link{SWAP.GetKTSP.Result}}
}


\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Luigi Marchionni \email{marchion@jhu.edu},
  Wikum Dinalankara \email{wdinala1@jhmi.edu}
}


\references{
See \link{switchBox} for the references.
}


\examples{

### Load gene expression data
data(trainingData)
data(testingData)

require(pROC)

### perform training and testing
result = SWAP.GetKTSP.TrainTestResults(matTraining, trainingGroup, 
  matTesting, testingGroup, featureNo=100)

### view results
result$train

result$test

}

\keyword{KTSP, classification, prediction}
