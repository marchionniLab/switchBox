\name{SWAP.PlotKTSP.TrainTestROC}

\alias{SWAP.PlotKTSP.TrainTestROC}

\title{Plots an ROC curve for training and testing results.}

\description{
  Given the output from SWAP.GetKTSP.TrainTestResults(), 
  plots the training and testing ROC curves.
}

\usage{
SWAP.PlotKTSP.TrainTestROC(result, colors=c(), legends=c(), ...)
}

\arguments{
  \item{result}{is either the output from \code{SWAP.GetkTSPTrainTestResults},
  or if manually prepared, a list with items \code{trainroc} 
  and \code{testroc} items, where each is an ROC object 
  produced by the \code{pROC} library.}
  \item{colors}{is a character vector indicating the color to be
  used for each curve.}
  \item{legends}{is a character vector providing any additional
  information to be appended to each curve label in 
  the legend.}
}

\value{
  Produces a plot with two ROC curves
  corresponding to training results and 
  testing/validation results.
}

\seealso{
  \code{\link{SWAP.GetkTSPTrainTestResults}}
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

### plot ROC curves
SWAP.PlotKTSP.TrainTestROC(result)

}

\keyword{KTSP, classification, prediction}
