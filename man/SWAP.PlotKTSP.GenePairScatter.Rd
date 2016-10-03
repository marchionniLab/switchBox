\name{SWAP.PlotKTSP.GenePairScatter}

\alias{SWAP.PlotKTSP.GenePairScatter}

\title{Make a scatter plot of two features.}

\description{
  Makes a scatter plot of a pair of features/genes.
}

\usage{
SWAP.PlotKTSP.GenePairScatter(inputMat, Groups, 
  classes, genes, colors=c(), legends=c(), ...)
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
  \item{colors}{is a character vector indicating the color to be
  used for each phenotype.}
  \item{legends}{is a character vector providing any additional
  information to be appended to the phenotype label in 
  the legend.}
}

\value{
  Produces a scatter plot containing points for each
  sample colored by the phenotype, with two axes being
  the measurements for the given two features.
}

\seealso{
  \code{\link{SWAP.PlotKTSP.Genes}}
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
SWAP.PlotKTSP.GenePairScatter(matTraining, trainingGroup, levels(trainingGroup), classifier$TSPs)

}

\keyword{KTSP, classification, prediction}
