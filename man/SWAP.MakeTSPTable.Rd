\name{SWAP.MakeTSPTable}

\alias{SWAP.MakeTSPTable}

\title{Make a table of TSPs in order of TSP score.}

\description{
  Given the output from \code{SWAP.CalculateScores} and a number \code{maxk}, 
  makes a table of the top maxk pairs. The output of this function
  can be provided to a k-selection function such as \code{SWAP.Kby.Ttest} 
  or \code{SWAP.Kby.Measurement} to test out different k-selection methods.
}

\usage{
SWAP.MakeTSPTable(Scores, maxk, disjoint = TRUE)
}

\arguments{
  \item{Scores}{is the output of a scoring function
  such as \code{SWAP.CalculateScores} containing a vector of TSP scores.}
  \item{maxk}{is an integer: the number of pairs to select.}
  \item{disjoint}{a logical indicating whether only disjoint pairs
  should be selected or not.}
}

\value{
  A data frame of \code{maxk} pairs, their score and \code{tieVote}.
}

\seealso{
  \code{\link{SWAP.Kby.Ttest}},
  \code{\link{SWAP.Kby.Measurement}}
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

### load gene expression data
data(trainingData)

### calculate scores
scores = SWAP.CalculateScores(matTraining, trainingGroup, featureNo=5)

### make top 5 pair table
SWAP.MakeTSPTable(scores, 5, FALSE)

}

\keyword{KTSP, classification, prediction}
