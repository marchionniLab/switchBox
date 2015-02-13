\name{SWAP.KTSP.Train}

\alias{SWAP.KTSP.Train}

\title{Function for training the K-TSP classifier.}

\description{
  \code{SWAP.KTSP.Train} trains a binary K-TSP classifier.
  The classifiers resulting from using this function can be
  passed to \code{\link{SWAP.KTSP.Classify}}
  for samples classification.
}

\usage{
SWAP.KTSP.Train(inputMat, phenoGroup, krange = c(3, 5, 7:10),
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
  \item{krange}{an integer (or a vector of integers)
    defining the candidate number of  
    Top Scoring Pairs (TSPs) from which the algorithm chooses
    to build the final classifier. The algorithm uses the mechanism
    in Afsari et al (AOAS, 2014) to select the number of pairs and pair of features. }
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
  \item{$score}{scores TSP for the top \code{k} TSPs.}
  \item{$label}{The class labels. These labels correspond to
    the \code{phenoGroup} factor lelves and will be used
    lable any new sample classified by the \code{SWAP.KTSP.Classify}
    function.}
}

\seealso{
  \code{\link{SWAP.KTSP.Classify}},
  \code{\link{SWAP.Filter.Wilcoxon}},
  \code{\link{SWAP.CalculateSignedScore}}
}


\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Luigi Marchionni \email{marchion@jhu.edu}
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
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup, krange=c(3, 5, 8:15))

### Show the classifier
classifier


##################################################
### Train another classifier from the top 4 best features 
### according to the deafault  filtering function
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      FilterFunc=SWAP.Filter.Wilcoxon, featureNo=4)

### Show the classifier
classifier


##################################################
### To use all features "FilterFunc" must be set to NULL
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup, FilterFunc=NULL)

### Show the classifier
classifier


##################################################
### Train a classifier using and alternative filtering function.
### For instance we can use the a  "t.test"  to selec the features 
### with an absolute t-statistics larger than a specified quantile
topRttest <- function(situation, data, quant = 0.75) {
	out <- apply(data, 1, function(x, ...) t.test(x ~ situation)$statistic )
	names(out[ abs(out) > quantile(abs(out), quant) ])
}

### Show the top features selected 
topRttest(trainingGroup, matTraining, quant=0.95)

### Train a classifier using the alternative filtering function
### and also define the maximum number of TSP using "krange"
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      FilterFunc = topRttest, quant = 0.75, krange=c(15:30) )

### Show the classifier
classifier


##################################################
### Training with restricted pairs

### Define a set of specific pairs to be used for classifier development
### For this example we will a random set of features
### In a real example these pairs should  be provided by the user.
set.seed(123)
somePairs <- matrix(sample(rownames(matTraining), 6^2, replace=FALSE), ncol=2)
head(somePairs, n=3)
dim(somePairs)

### Train a classifier using the restricted feature pairs and the default filtering
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      RestrictedPairs = somePairs, krange=3:16)

### Show the classifier
classifier

}

\keyword{KTSP, classification, prediction}
