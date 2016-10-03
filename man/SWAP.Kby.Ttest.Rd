\name{SWAP.Kby.Ttest}

\alias{SWAP.Kby.Ttest}

\title{K selection for a kTSP classifier.}

\description{
  \code{SWAP.Kby.Ttest} can be supplied to a
  kTSP classifier training function to select
  an optimal k via performing t-tests.
}

\usage{
SWAP.Kby.Ttest(inputMat, phenoGroup, 
  scoreTable, classes, krange, 
  k_opts=list())
}

\arguments{
  \item{inputMat}{is a numerical matrix containing the
    measurements (\emph{e.g.}, gene expression data)
    to be used to build the K-TSP classifier.}
  \item{phenoGroup}{is a factor with two levels containing
    the phenotype information used to train the K-TSP classifier.}
  \item{scoreTable}{a data frame output of \code{SWAP.MakeTSPTable}
  containing TSPs and the accuracy of individuals pairs over
  the training data.}
  \item{classes}{ is a character vector of length 2
    providing the phenotype class labels (case followed 
    by control). If NULL, the levels of phenoGroup 
    will be taken as the labels.
  }
  \item{krange}{an integer (or a vector of integers)
    defining the candidate number of  
    Top Scoring Pairs (TSPs) from which the algorithm chooses
    to build the final classifier.}
  \item{k_opts}{is not used and is left for
    conforming to the arguments of \code{k_selection_fn}.
  }
}

\value{
  A vector of indices of length k indicating
  which pairs from \code{scoreTable} should be
  selected.
}

\seealso{
  \code{\link{SWAP.Kby.Measurement}},
  \code{\link{SWAP.MakeTSPTable}}
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
