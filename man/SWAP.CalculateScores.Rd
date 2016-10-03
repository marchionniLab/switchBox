\name{SWAP.CalculateScores}

\alias{SWAP.CalculateScores}

\title{Function to calculate the pair-wise scores with any given 
  score function.}

\description{
  \code{SWAP.CalculateScores} calculates the pair-wise scores
  between features pairs. The user may pass a filtering function
  to reduce the number of starting features, or provide a restricted
  set of pairs to limit the reported scores to this list. The user 
  can also pass a score-calculating function by either passing one of the
  scoring functions available in the package(i.e. 
  \code{SWAP.Calculate.SignedTSPScores} and \code{SWAP.Calculate.BasicTSPScores}) 
  or a custom function.
}

\usage{
SWAP.CalculateScores(inputMat, phenoGroup, classes = NULL, FilterFunc = SWAP.Filter.Wilcoxon, 
  RestrictedPairs = NULL, handleTies = FALSE, verbose = FALSE, 
  score_fn = signedTSPScores, score_opts = list(), ...)
}

\arguments{
  \item{inputMat}{is a numerical matrix containing the
    measurements (\emph{e.g.}, gene expression data)
    to be used to build the K-TSP classifier.
    The columns represent samples and the
    rows represent the features (\emph{e.g.}, genes).
    The number of columns must agree
    with the length of \code{phenoGroup}.
    Note that \code{rownames(inputMat)} will be construed
    as feature names (\emph{e.g.}, gene symbols)
    in all subsequent analyses.}
  \item{phenoGroup}{ is a factor containing the training
    phenotypes with two levels.}
  \item{classes}{ is a character vector of length 2
    providing the phenotype class labels (case followed 
    by control). If NULL, the levels of phenoGroup 
    will be taken as the labels.
  }
  \item{FilterFunc}{is a filtering function to reduce the
    starting number of features to be used to identify the
    Top Scoring Pairs (TSPs). The default filter is based on
    the Wilcoxon rank-sum test
    and alternative filtering functions can be passed too
    (see \code{SWAP.Filter.Wilcoxon} for details). Note the filtering 
	function must return feature names, i.e. a subset of \code{rownames(inputMat)}. }
  \item{RestrictedPairs}{is a character matrix with two columns
    containing the feature pairs to be considered for score calculations.
    Each row should contain a pair of feature names matching the
    \code{rownames(inputMat)}.
    If RestrictedPairs is missing all available feature pairs
    will be considered.}
  \item{handleTies}{is a logical value indicating whether tie handling
   should be enabled or not. FALSE by default.
    }
  \item{verbose}{is a logical value indicating whether status
    messages will be printed or not throughout the function. 
    FALSE by default.
    }
  \item{score_fn}{is a function for calculating TSP scores. 
    By default, the signed TSP scores as calculated by 
    \code{SWAP.Calculate.SignedTSPScores} will be used. The user can 
    also provide \code{SWAP.Calculate.BasicTSPScores} to obtain basic
    TSP scores. The output of any custom function should correspond
    to the same strucure as the output from these two functions.}
  \item{score_opts}{is a list of additional variables that
    will be passed on to the scoring function as the \code{score_opts}
    argument.}
  \item{...}{Additional argument passed to the filtering
  function \code{FilterFunc}.}
}

\value{
The output is a list containing the following items:
\item{labels}{the levels (phenotypes) in \code{phenoGroup}.}
\item{score}{is a vector containing the pair-wise scores.}
}

\seealso{
  See \code{\link{SWAP.KTSP.Train}},
  \code{\link{SWAP.Filter.Wilcoxon}},
  \code{\link{SWAP.Calculate.BasicTSPScores}},
  \code{\link{SWAP.Calculate.SignedTSPScores}},
  and \code{\link{SWAP.KTSP.Statistics}}.
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

### Compute the scores using all features (a matrix will be returned)
scores <- SWAP.CalculateScores(matTraining,  trainingGroup, FilterFunc=NULL)


}


\keyword{Pairwise score}
