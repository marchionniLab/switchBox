\name{SWAP.PlotKTSP.Genes}

\alias{SWAP.PlotKTSP.Genes}

\title{Plot features seperated by phenotype}

\description{
  Makes line plots of one or more features seperated by phenotype.
}

\usage{
SWAP.PlotKTSP.Genes(inputMat, Groups, classes, genes, 
  colors=c(), legends=c(), ...)
}

\arguments{
  \item{inputMat}{is a matrix of data with rows being the
  features (such as gene names, if the matrix if gene expression
  data) and columns being the samples.}
  \item{Groups}{is a factor or a vector providing the phenotype class
  each sample belongs to. It should correspond to the order of samples
  given by the columns of \code{inputMat}.}
  \item{classes}{is a vetor of length 2 providing the two phenotype or
  class labels of \code{Groups}.}
  \item{genes}{is a vector of length one or more providing the names
  (from the rownames of \code{inputMat}) of the features to be plotted.}
}

\value{
  Produces a plot containing a line for each feature plotted,
  the x-axis being the ordering of samples and the y-axis being
  the measured value (such as gene expression).
}

\seealso{
  \code{\link{SWAP.PlotKTSP.GenePairScatter}}
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

### train 1-TSP
classifier = SWAP.Train.1TSP(matTraining, trainingGroup)

### plot top pair
SWAP.PlotKTSP.Genes(matTraining, trainingGroup, levels(trainingGroup), classifier$TSPs)

}

\keyword{KTSP, classification, prediction}
