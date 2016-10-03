\name{SWAP.PlotKTSP.GenePairBoxplot}

\alias{SWAP.PlotKTSP.GenePairBoxplot}

\title{Plots a feature pair as boxplots.}

\description{
  Plots two genes or features as a pair 
  of boxplots; optionally, individual samples 
  can be plotted on top of the boxplots as points; 
  for this points can be colored by either gene, or 
  class, or whether first gene < second gene is TRUE
  or FALSE for each sample.
}

\usage{
SWAP.PlotKTSP.GenePairBoxplot(genes, inputMat, Groups=NULL, 
  classes=NULL, points=FALSE, point_coloring="byGene", 
  colors=c(), point_colors=c(), ...)
}

\arguments{
  \item{genes}{is a vector of length two providing the pair
  (from the rownames of \code{inputMat}) of features to be plotted.}
  \item{inputMat}{is a matrix of data with rows being the
  features (such as gene names, if the matrix if gene expression
  data) and columns being the samples.}
  \item{Groups}{is a factor or a vector providing the phenotype class
  each sample belongs to. It should correspond to the order of samples
  given by the columns of \code{inputMat}.}
  \item{classes}{is a vetor of length 2 providing the two phenotype or
  class labels of \code{Groups}.}
  \item{points}{is a logical value indicating whether to overlay
  the boxplot with points for individual samples or not.}
  \item{point_coloring}{can be either 'byGene' or 'byClass' indicating
  whether to color the points by gene/feature or by phenotype. A third
  option is 'byDirection' indicating to color the points by whether 
  the first gene is less than the second gene.}
  \item{colors}{is a character vector indicating the color to be
  used for each boxplot.}
  \item{point_colors}{is a character vector indicating the color to be
  used for the points.}
  \item{...}{any further arguments are supplied to the \code{boxplot}
  function.}
}

\value{
  Produces a pair of boxplots indicating the distribution of 
  the measured values for the pair of features/genes.
}

\seealso{
  \code{\link{SWAP.PlotKTSP.GenePairClassesBoxplot}}
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
SWAP.PlotKTSP.GenePairBoxplot(classifier$TSPs, matTraining, 
  points=TRUE, point_coloring="byGene")

}

\keyword{KTSP, classification, prediction}
