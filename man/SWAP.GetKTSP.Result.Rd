\name{SWAP.GetKTSP.Result}

\alias{SWAP.GetKTSP.Result}

\title{Function for prediction followed by computing 
various performance measures related to prediction.}

\description{
  Given a kTSP classifier and data matrix and class labels, calculates the predictions and 
  vote sums and then applies \code{SWAP.GetKTSP.PredictionStats}.
}

\usage{
SWAP.GetKTSP.Result(classifier, inputMat, Groups, 
  classes=NULL, predictions=FALSE, decision_values=FALSE)
}

\arguments{
  \item{classifier}{a k-TSP classifier computed using
  \code{SWAP.KTSP.Train} or \code{SWAP.Train.1TSP}.}
  \item{inputMat}{is a matrix of data with rows being the
  features (such as gene names, if the matrix if gene expression
  data) and columns being the samples.}
  \item{Groups}{is a factor or a vector providing the phenotype class
  each sample belongs to. It should correspond to the order of samples
  given by the columns of \code{inputMat}.}
  \item{classes}{is a vetor of length 2 providing the two phenotype or
  class labels of \code{Groups}.}
  \item{predictions}{is a logical indicating whether to return the 
  predictions or not.}
  \item{decision_values}{is a logical indicating whether to return the 
  decision values or not.}
}

\value{
  A list with items:
  \item{stats}{A vector providing accuracy, sensitivity, specificity, balanced accuracy,
  and AUC.}
  \item{roc}{An ROC curve object produced by the \code{pROC} package.}
}

\seealso{
  \code{\link{SWAP.GetKTSP.PredictionStats}}
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

### train 1-TSP
classifier = SWAP.Train.1TSP(matTraining, trainingGroup)

### get performance results
SWAP.GetKTSP.Result(classifier, matTesting, testingGroup)$stats


}

\keyword{KTSP, classification, prediction}
