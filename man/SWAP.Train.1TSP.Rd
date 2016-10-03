\name{SWAP.Train.1TSP}

\alias{SWAP.Train.1TSP}

\title{Function for training the 1-TSP classifier.}

\description{
  \code{SWAP.Train.1TSP} trains a binary TSP classifier 
  with a single top scoring pair. The classifiers resulting 
  from using this function can be
  passed to \code{\link{SWAP.KTSP.Classify}}
  for samples classification.
}

\usage{
SWAP.Train.1TSP(inputMat, phenoGroup, classes = NULL, 
  FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = NULL,
  handleTies = FALSE, disjoint = TRUE,
  score_fn = signedTSPScores, score_opts = NULL, 
  verbose = FALSE, ...)
}

\arguments{
  \item{inputMat}{is a numerical matrix containing the
    measurements (\emph{e.g.}, gene expression data)
    to be used to build the K-TSP classifier.
    The columns represent samples and the
    rows represent the features (\emph{e.g.}, genes).
    The number of columns must agree
    with the length of \code{phenoGroup}.
    Note that \code{rownames(inputMat)} will be used
    as the feature names (\emph{e.g.}, gene symbols)
    in all subsequent analyses.}
  \item{phenoGroup}{is a factor with two levels containing
    the phenotype information used to train the K-TSP classifier.
    In order to identify the best TSP to be included in the classifier,
    the features contained in \code{inputMat} will be compared
    between the two groups defined by this factor.
    Levels from \code{phenoGroup} will be also used
    to reorder the features in each TSP such as the first feature
    is larger than the second one in the group corresponding
    to first level, and \emph{vice-versa}.}
  \item{classes}{ is a character vector of length 2
    providing the phenotype class labels (case followed 
    by control). If NULL, the levels of phenoGroup 
    will be taken as the labels.
  }
  \item{FilterFunc}{is a filtering function to reduce the
    starting number of features to be used to identify the
    Top Scoring Pairs (TSP). The default filter is differential 
    expression test based on the Wilcoxon rank-sum test
    and alternative filtering functions can be passed too
    (see \code{SWAP.Filter.Wilcoxon} for details). The output of the function 
    must be subset of \code{rownames(inputMat)}}
  \item{RestrictedPairs}{is a character matrix with two columns
    containing the feature pairs to be considered for score calculations.
    Each row should contain a pair of feature names matching the
    rownames of \code{inputMat}.
    If RestrictedPairs is missing all available feature pairs
    will be considered.}
  \item{handleTies}{is a logical value indicating whether tie handling
   should be enabled or not. FALSE by default.
    }
  \item{disjoint}{is a logical value indicating whether only disjoint
  pairs should be considered in the final set of selected pairs; i.e. all
  features occur only once among the set of TSPs.}
  \item{score_fn}{is a function for calculating TSP scores. 
    By default, the signed TSP scores as calculated by 
    \code{SWAP.Calculate.SignedTSPScores} will be used. The user can 
    also provide \code{SWAP.Calculate.BasicTSPScores} to obtain basic
    TSP scores. The output of any custom function should correspond
    to the same strucure as the output from these two functions.}
  \item{score_opts}{is a list of additional variables that
    will be passed on to the scoring function as the \code{score_opts}
    argument.}
  \item{verbose}{is a logical value indicating whether status
    messages will be printed or not throughout the function. 
    FALSE by default.
    }
  \item{...}{Additional argument passed to the filtering
  function \code{FilterFunc}.}
}

\value{
  The TSP classifier, in the form of a list,
  which contains the following components:
  \item{name}{The classifier name.}
  \item{TSPs}{A 1 by 2 matrix, containing
    the feature names for the selected TSP. These names
    correspond to the \code{rownames(inputData)}.}
  \item{score}{scores TSP for the top TSP.}
  \item{label}{the class labels. These labels correspond to
    the \code{phenoGroup} factor lelves and will be used
    lable any new sample classified by the \code{SWAP.KTSP.Classify}
    function.}
  \item{tieVote}{indicates which class the pair would vote
  for in case of a tie.
  }
}

\seealso{
  \code{\link{SWAP.KTSP.Classify}},
  \code{\link{SWAP.Filter.Wilcoxon}},
  \code{\link{SWAP.CalculateSignedScore}}
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

##################################################
### Load gene expression data for the training set
data(trainingData)


### Show group variable for the TRAINING set
table(trainingGroup)


##################################################
### Train a classifier using default filtering function based on the Wilcoxon test
classifier <- SWAP.Train.1TSP(matTraining, trainingGroup)

### Show the classifier
classifier

}

\keyword{KTSP, classification, prediction}
