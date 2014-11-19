\name{SWAP.Filter.Wilcoxon}

\alias{SWAP.Filter.Wilcoxon}

\title{Statistical feature filtering based on Wilcoxon test on the ranks of expressions.}

\description{\code{SWAP.Filter.Wilcoxon} filters the features to top
 differential expressed 
  to be used for KTSP classifier implementation.
}

\usage{
SWAP.Filter.Wilcoxon(phenoGroup, inputMat, featureNo = 100, UpDown = TRUE)
}

\arguments{
  \item{phenoGroup}{a factor with levels containing training labels
    for the phenotype of interest.}
  \item{inputMat}{a numerical matrix containing  feature
    measurements to be used to implement the classifier
    (\emph{e.g.}, the set of gene expression values).
    The columns of this matrix correspond to samples and
    must correspond to \code{phenoGroup}.
    The rows represent the features and
    \code{rownames(inputMat)} will be used as feature names.}
  \item{featureNo}{an integer specifying the number
    of different features to be returned.}
  \item{UpDown}{logical value specifying whether
    an equal proportion of features displaying opposite
    change across the two phenotypes should be returned
    (\emph{e.g.}an equal number of up- and down-regulated genes).}
}

\value{The names of the features that survived the statistical filtering, i.e. differential expressed features.}

 
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

### Load gene expression data for the training set
data(trainingData)

### Return equal numbers of up- and down- regulated features (default)
SWAP.Filter.Wilcoxon(trainingGroup, matTraining, featureNo=10)

### Return the top 10 differentially expressed features irrispective to
### the direction of change.
### By setting the argument 'UpDown' equal to FALSE the number of
### up- and down- regulated features can be different
SWAP.Filter.Wilcoxon(trainingGroup, matTraining, featureNo=10, UpDown=FALSE)

}
