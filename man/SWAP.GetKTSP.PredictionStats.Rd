\name{SWAP.GetKTSP.PredictionStats}

\alias{SWAP.GetKTSP.PredictionStats}

\title{Function for computing various performance
measures related to prediction.}

\description{
  Given a list of predicted labels and true labels, provides accuracy, sensitivity, specificity, balanced
  accuracy (i.e. (sensitivity+specificity)/2 ), and AUC if decision values are given.
}

\usage{
SWAP.GetKTSP.PredictionStats(predictions, truth, classes=NULL, 
  decision_values=NULL)
}

\arguments{
  \item{predictions}{is a vector or factor of predicted 
  classes.}
  \item{truth}{is a vector or factor of the true class
  labels.}
  \item{classes}{ is a character vector of length 2
    providing the phenotype class labels (case followed 
    by control). If NULL, the levels of phenoGroup 
    will be taken as the labels.
  }
  \item{decision_values}{is a vector providing the decision values
  (such as sum of votes from a k-TSP classifier). Will be used
  to compute AUC if provided.}
}

\value{
  A vector providing accuracy, sensitivity, specificity, and balanced accuracy
  , and if \code{decision_values} is prodvided, area under the ROC curve (AUC).
}

\seealso{
  \code{\link{SWAP.KTSP.Classify}}
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

### train 1-TSP
classifier = SWAP.Train.1TSP(matTraining, trainingGroup)
predictions = SWAP.KTSP.Classify(matTesting, classifier)

### get performance results
SWAP.GetKTSP.PredictionStats(predictions, testingGroup)

}

\keyword{KTSP, classification, prediction}
