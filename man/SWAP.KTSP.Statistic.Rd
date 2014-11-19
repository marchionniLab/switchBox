\name{SWAP.KTSP.Statistics}

\alias{SWAP.KTSP.Statistics}

\title{Function computing TSP votes (comparisons) and combine their votes. The default is the kTSP statistics, sum of the votes.}

\description{
  \code{SWAP.KTSP.Statistics} computes the votes in favor
  of one of the classes or the other for each TSP. This function
  also computes the final, combined, consensus of all
  TSP votes based on a specific decision rules. The default is the kTSP statistics, sum of the votes.
  }

\usage{
SWAP.KTSP.Statistics(inputMat, classifier, CombineFunc)
}

\arguments{
  \item{inputMat}{is a numerical matrix containing the
    measurements (\emph{e.g.}, gene expression data)
    to be used to compute the individual TSP votes
    and their consensus.
    like the matrix used for training classifier
    (in \code{\link{SWAP.KTSP.Train}} function), inputMatrix 
	rows represent the features and the columns represent the samples.
   }
  \item{classifier}{the classifier obtained by invoking 
    \code{\link{SWAP.KTSP.Train}}.}
  \item{CombineFunc}{ is the function used to combine
    the votes (i.e., comparisons) of individual TSPs contained in the classifier.
    By default, the consensus is the count of the votes taking into
    account the order of the features in each TSP.
    Using this argument alternative aggregating functions
    can be also passed  to \code{SWAP.KTSP.Statistics}
    as described below (see ``details'').}
}  


\value{
  A list containing the following two components:
  \item{statistics}{ a named vector containing the
    aggregated summary statistics computed by
    \code{CombineFunc}. The names correspond
    to samples and are derived from \code{colnames(inputMat)}.}
  \item{comparisons}{ a logical matrix containing the individual
    TSP votes (\code{TRUE} if the first pair  feature is larger
    then the second one, \code{FALSE} otherwise).
    The columns of this matrix correspond to TSP comparisons
    and are named accordingly using feature names derived
    from \code{rownames(inputMat)}.
    The columns of this matrix correspond to the samples
    and are named accordingly using \code{colnames(inputMat)}.}
}


\details{
  For each TSP in the KTSP classifier,
  \code{SWAP.KTSP.Statistics} computes the vote in favor of 
  one of classes or the other.
  This function also aggregates the individual TSP votes
  and computes a final consensus of all TSP votes
  based on specific combination rules.
  By default, this combination is achieved by counting the
  comparisons (votes) of TSPs as follows:
  If the first feature is larger than the second one,
  the TSP vote is positive, else the TSP vote is negative.
  Different combination rules can also be specified
  by defining an alternative combination function
  and by passing it to \code{SWAP.KTSP.Statistics}
  using the \code{CombineFunc} argument.
  A combination function takes
  as its input a logical vector \code{x} corresponding
  to the sample TSP comparisons
  (\code{TRUE} if the first feature in the pair is
  larger then the second, \code{FALSE}
  in the opposite case). The output of the
  \code{CombineFunction} is a single value
  summarizing the votes of all individual TSPs
  (see examples below).
  Note that \code{CombineFunction} function must operate on a logical
  vector as input and the outcome must be real value number.
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
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup,
			      FilterFunc = NULL, krange=8)

### Show the TSP in the classifier 
classifier$TSPs


##################################################
### Compute the TSP votes and combine them using various methods

### Here we will use the count of the signed TSP votes
ktspStatDefault <- SWAP.KTSP.Statistics(inputMat = matTraining,
    classifier = classifier)

### Here we will use the sum of the TSP votes
ktspStatSum <- SWAP.KTSP.Statistics(inputMat = matTraining,
    classifier = classifier, CombineFunc=sum)

### Here, for instance, we will apply a hard treshold equal to 2
ktspStatThreshold <- SWAP.KTSP.Statistics(inputMat = matTraining,
    classifier = classifier,  CombineFunc = function(x) sum(x) > 2 )

### Show components
names(ktspStatDefault)

### Show some of the votes
head(ktspStatDefault$comparisons[ , 1:2])

### Show default statistics
head(ktspStatDefault$statistics)

### Show statistics obtained using the sum
head(ktspStatSum$statistics)

### Show statistics obtained using the hard threshold
head(ktspStatThreshold)

### Make a heatmap showing the individual TSPs votes
colorForRows <- as.character(1+as.numeric(trainingGroup))
heatmap(1*ktspStatDefault$comparisons, scale="none",
    margins = c(10, 5), cexCol=0.5, cexRow=0.5,
    labRow=trainingGroup, RowSideColors=colorForRows)

}

\keyword{KTSP, classification, prediction}

