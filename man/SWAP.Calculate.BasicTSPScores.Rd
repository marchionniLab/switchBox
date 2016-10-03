\name{SWAP.Calculate.BasicTSPScores}

\alias{SWAP.Calculate.BasicTSPScores}

\title{Function to calculate basic TSP scores.}

\description{
  \code{SWAP.Calculate.BasicTSPScores} calculates basic TSP 
  scores. 
}

\usage{
SWAP.Calculate.BasicTSPScores(phenoGroup, inputMat1, 
  inputMat2 = NULL, classes = NULL, RestrictedPairs = NULL, 
    handleTies = FALSE, verbose = FALSE, score_opts=list())
}

\arguments{
  \item{phenoGroup}{ is a factor containing the training
    phenotypes with two levels.}
  \item{inputMat1}{is a numerical matrix containing the
    measurements (\emph{e.g.}, gene expression data)
    for choosing the first item of a top scoring pair.}
  \item{inputMat2}{is a numerical matrix containing the
    measurements for choosing the second item of a 
    top scoring pair. If \code{NULL}, \code{inputMat1}
    will be used for this.}
  \item{classes}{ is a character vector of length 2
    providing the phenotype class labels (case followed 
    by control). If NULL, the levels of phenoGroup 
    will be taken as the labels.
  }
  \item{RestrictedPairs}{is a character matrix with two columns
    containing the feature pairs to be considered for score calculations.
}
  \item{handleTies}{is a logical value indicating whether tie handling
   should be enabled or not. FALSE by default.
    }
  \item{verbose}{is a logical value indicating whether status
    messages will be printed or not throughout the function. 
    FALSE by default.
    }
  \item{score_opts}{is a list of additional variables that
    will be passed on to the scoring function.
    }
}

\value{
The output is a list containing the following items:
\item{labels}{the levels (phenotypes) in \code{phenoGroup}.}
\item{score}{is a vector containing the pair-wise scores.}
\item{tieVote}{is a vector indicating the class the pair
would vote for in the case of a tie.}
}

\seealso{
  See \code{\link{SWAP.Calculate.SignedTSPScores}}
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

### Load gene expression data for the training set
data(trainingData)

### Show group variable for the TRAINING set
table(trainingGroup)

### Compute the scores
scores = SWAP.Calculate.BasicTSPScores(trainingGroup, matTraining[1:3, ])

# View the scores
scores$score

}


\keyword{Pairwise score}
