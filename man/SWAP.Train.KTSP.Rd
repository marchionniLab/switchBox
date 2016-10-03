\name{SWAP.Train.KTSP}

\alias{SWAP.Train.KTSP}

\title{Function for training the K-TSP classifier.}

\description{
  \code{SWAP.Train.KTSP} trains a binary K-TSP classifier.
  The classifiers resulting from using this function can be
  passed to \code{\link{SWAP.KTSP.Classify}}
  for samples classification.
}

\usage{
SWAP.Train.KTSP(inputMat, phenoGroup, classes = NULL, krange = 2:10,
  FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = NULL, 
  handleTies = FALSE, disjoint = TRUE,
  k_selection_fn = KbyTtest, k_opts = list(), score_fn = signedTSPScores, 
  score_opts = NULL, verbose = FALSE, ...)
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
  \item{krange}{an integer (or a vector of integers)
    defining the candidate number of  
    Top Scoring Pairs (TSPs) from which the algorithm chooses
    to build the final classifier. The algorithm uses the mechanism
    in Afsari et al (AOAS, 2014) to select the number of pairs and pair of features.
    Default is the range from 2 to 10.}
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
  \item{k_selection_fn}{is a function for selecting the optimal k
  once the TSP scores have been calculated for all the candidate pairs.
  This can be either \code{SWAP.Kby.Measurement} or \code{SWAP.Kby.Ttest}(default),
  or a user defined function.
  }
  \item{k_opts}{a list of additional arguments to be passed on to a custom
  k selection function.}
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
  The KTSP classifier, in the form of a list,
  which contains the following components:
  \item{name}{The classifier name.}
  \item{TSPs}{A \code{k} by 2 matrix, containing
    the feature names for each TSP. These names
    correspond to the \code{rownames(inputData)}.
    In this matrix each row corresponds to a specific TSP.
    For each TSP (\emph{i.e.} row in the \code{TSPs} matrix)
    the order of the features is such that the first one
    is on average smaller than the second one in the
    phenotypic group defined by the first levels
    of the \code{phenoGroup} factor
    and \emph{vice-versa}.
    The algorithm uses the mechanism
    in Afsari et al (2014) to select the number
    of pairs and pair of features.}
  \item{score}{scores TSP for the top \code{k} TSPs.}
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
classifier <- SWAP.Train.KTSP(matTraining, trainingGroup)

### Show the classifier
classifier


##################################################
### Train another classifier from the top 4 best features 
### according to the deafault  filtering function
classifier <- SWAP.Train.KTSP(matTraining, trainingGroup,
			      FilterFunc=SWAP.Filter.Wilcoxon, featureNo=4)

### Show the classifier
classifier

}

\keyword{KTSP, classification, prediction}
