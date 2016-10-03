\name{SWAP.PlotKTSP.GenePairClassesBoxplot}

\alias{SWAP.PlotKTSP.GenePairClassesBoxplot}

\title{Plots a feature pair as seperated by 
class as boxplots.}

\description{
  Plots two genes or features, each as a pair 
  of boxplots seperated to two classes or phenotypes.
}

\usage{
SWAP.PlotKTSP.GenePairClassesBoxplot(genes, inputMat, Groups, 
  classes=NULL, points=FALSE, ordering="byGene", 
  colors=c(), point_colors=c(), point_directions=FALSE, ...)
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
  \item{ordering}{can be either 'byGene' or 'byClass' respectively
  indicating whether to plot two adjacent boxplots for each class/phenotype
  or two adjacent boxplots for each gene/features.}
  \item{colors}{is a character vector indicating the color to be
  used for each class or gene boxplots.}
  \item{point_colors}{is a character vector indicating the color to be
  used for the points.}
  \item{point_directions}{is a logical indicating whether to color
  the points by whether the first gene is less than the second gene.}
  \item{...}{any further arguments are supplied to the \code{boxplot}
  function.}
}

\value{
  Produces a pair of boxplots indicating the distribution of 
  the measured values for the pair of features/genes.
}

\seealso{
  \code{\link{SWAP.PlotKTSP.GenePairBoxplot}}
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
SWAP.PlotKTSP.GenePairClassesBoxplot(classifier$TSPs, matTraining,
  trainingGroup, levels(trainingGroup),
  points=TRUE, ordering="byGene")


}

\keyword{KTSP, classification, prediction}
