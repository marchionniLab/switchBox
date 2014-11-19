\name{SWAP.CalculateSignedScore}

\alias{SWAP.CalculateSignedScore}

\title{Function to calculate the pair-wise scores.}

\description{
  \code{SWAP.CalculateSignedScore} calculates the pair-wise scores
  between features pairs. The user may pass a filtering function
  to reduce the number of starting features, or provide a restricted
  set of pairs to limit the reported scores to this list.
}

\usage{
SWAP.CalculateSignedScore(inputMat, phenoGroup,
  FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs, ...)
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
  \item{...}{Additional argument passed to the filtering
  function \code{FilterFunc}.}
}

\value{
The output is a list containing the following items:
\item{labels}{the levels (phenotypes) in \code{phenoGroup}.}
\item{P}{a matrix or a vector containing the probability of
  comparisons in samples with phenotype equal to \code{label[1]}.
  In case \code{RestrictedPairs} is not given,
  \code{P} is a matrix and
  \code{P[i,j] = P(InputMat1[i,] < InputMat2[j,] |
    phenoGroup == $label[1])}.
  In case \code{RestrictedPairs} is given,
  \code{P[k] = P(InputMat1[RestrictedPairs[k,1],] <
    InputMat2[RestrictedPairs[k,2],] | phenoGroup == $label[1])}.}
\item{Q}{a matrix or a vector containing the probability
  of comparisons in samples with phenotype \code{label[2]}.}
\item{score}{a matrix or a vector containing the pair-wise scores.
  Basically, \code{score = P - Q + C}.
  The \code{C} term is the tie breaker and proportion to
  the secondary score to avoid the ties.}
Note that the \code{P}, \code{Q}, and \code{score}
list elements are matrices when scores are computed
for all possible feature pairs, while they are vectors
when scores are computed for restricted pairs
defined by  \code{RestrictedPairs}.
}

\seealso{
  See \code{\link{SWAP.KTSP.Train}},
  \code{\link{SWAP.Filter.Wilcoxon}},
  and \code{\link{SWAP.KTSP.Statistics}}.
}


\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Luigi Marchionni \email{marchion@jhu.edu}
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
scores <- SWAP.CalculateSignedScore(matTraining,  trainingGroup, FilterFunc=NULL, )

### Show scores
class(scores)
dim(scores$score)

### Get the scores for a couple of features
diag(scores$score[ 1:3 , 5:7 ])

### Compute the scores using the default filtering function for 20 features
scores <- SWAP.CalculateSignedScore(matTraining, trainingGroup, featureNo=20)

### Show scores
dim(scores$score)

### Creating some random pairs
set.seed(123)
somePairs <- matrix(sample(rownames(matTraining), 25, replace=FALSE), ncol=2)

### Compute the scores for restricted pairs (a vector will be returned)
scores <- SWAP.CalculateSignedScore(matTraining, trainingGroup,
				    FilterFunc = NULL, RestrictedPairs = somePairs )

### Show scores
class(scores$score)
length(scores$score)

}


\keyword{Pairwise score}
