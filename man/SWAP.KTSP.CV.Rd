\name{SWAP.KTSP.CV}

\alias{SWAP.KTSP.CV}

\title{Performs k-fold cross validation.}

\description{
  Partitions the data into k folds and applies 
  \code{SWAP.GetKTSP.TrainTestResults} for each fold.
  Then it combines prediction votes by dividing the vote 
  sums by the number of TSPs in each fold to produce an
  overall cross-validation result.
}

\usage{
SWAP.KTSP.CV(inputMat, Groups, classes = NULL, k = 4, 
  folds = NULL, randomize = TRUE, ...)
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
  \item{k}{an integer giving the number of folds to use.}
  \item{folds}{a list containing the samples to be used in each fold;
  if \code{NULL}, the data will be split into \code{k} folds
  maintaning the proportions between the classes.}
  \item{randomize}{is a logical indicating whether to randomize the sample
  order before diving into \code{k} folds.}
  \item{...}{any further arguments to be passed on for k-TSP training.}
}

\value{
  A list with items:
  \item{folds}{A list containing the sample indices used in each fold.}
  \item{cv}{A list containing the classifier, training
  performance and testing performance for each fold.}
  \item{stats}{Overall cross-validation performance.}
  \item{roc}{ROC curve object for overall cross-validation performance.}
}

\seealso{
  \code{\link{SWAP.KTSP.LOO}}
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
result = SWAP.KTSP.CV(matTraining, trainingGroup, featureNo=100)

### print results
result$stats


}

\keyword{KTSP, classification, prediction}
