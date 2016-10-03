\name{SWAP.PlotKTSP.Votes}

\alias{SWAP.PlotKTSP.Votes}

\title{Plots a heatmap of k-TSP votes.}

\description{
  Given a k-TSP classifer and a matrix of data,
  plots a heatmap of the votes of the pairs 
  computed on the given data.
}

\usage{
SWAP.PlotKTSP.Votes(classifier, inputMat, 
  Groups=NULL, CombineFunc, ...)
}

\arguments{
  \item{classifier}{is a k-TSP classifier produced by
  \code{SWAP.KTSP.Train}.}
  \item{inputMat}{is a matrix of data with rows being the
  features (such as gene names, if the matrix if gene expression
  data) and columns being the samples.}
  \item{Groups}{is a factor or a vector providing the phenotype class
  each sample belongs to. It should correspond to the order of samples
  given by the columns of \code{inputMat}. These phenotype
  labels will be added to the x-axis of the heatmap.}
  \item{CombineFunc}{is a function corresponding to the
  \code{CombineFunc} argument of the \code{SWAP.KTSP.Classify}
  function.}
}

\value{
  Produces a heatmap where the color indicates
  a vote of 1 or 0 for a given sample by a
  top scoring pair.
}

\seealso{
  \code{\link{SWAP.KTSP.Train}}
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

}

\keyword{KTSP, classification, prediction}
