\name{SWAP.KTSP.LOO}

\alias{SWAP.KTSP.LOO}

\title{Performs leave one out cross validation.}

\description{
  Performs leave one out cross validation; then it combines 
  prediction votes by dividing the vote sums by the number
  of TSPs in each fold to produce an overall 
  cross-validation result.
}

\usage{
SWAP.KTSP.LOO(inputMat, Groups, classes = NULL, ...)
}

\arguments{
  \item{inputMat}{is a matrix of data with rows being the
  features (such as gene names, if the matrix if gene expression
  data) and columns being the samples.}
  \item{Groups}{is a factor or a vector providing the phenotype class
  each sample belongs to. It should correspond to the order of samples
  given by the columns of \code{inputMat}.}
  \item{classes}{is a vetor of length 2 providing the two phenotype or
  class labels.}
  \item{...}{any further arguments to be passed on for k-TSP training.}
}

\value{
  A list with items:
  \item{loo}{A list containing the classifier, training
  performance and testing performance for each fold.}
  \item{decision_values}{Decision values obtained for each left-out
  sample.}
  \item{predictions}{Predicted classes for each left-out
  sample.}
  \item{stats}{Overall peformance results.}
  \item{roc}{ROC curve object for overall performance.}
}

\seealso{
  \code{\link{SWAP.KTSP.CV}}
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

### perform leave one out cross-validation
result = SWAP.KTSP.LOO(matTraining, trainingGroup, featureNo=100)

### print results
result$stats

}

\keyword{KTSP, classification, prediction}
