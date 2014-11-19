\name{SWAP.KTSP.Classifiy}

\alias{SWAP.KTSP.Classify}

\title{Function to classify samples using a KTSP classifier.}

\description{
  \code{SWAP.KTSP.Classify} classifies new test samples
  using KTSP coming out of the function \code{\link{SWAP.KTSP.Train}}.
}

\usage{
SWAP.KTSP.Classify(inputMat, classifier, DecisionFunc)
}

\arguments{
  \item{inputMat}{is a numerical matrix containing the
    measurements (\emph{e.g.}, gene expression data)
    to be used with a K-TSP classifier to classify the samples
    in a specific class or the other.
    In this numerical matrix the columns represent the samples
    and the rows represent the features (\emph{e.g.}, genes)
    used by the classification rule.
    Note that \code{rownames(inputMat)} will be used
    to select the features (\emph{e.g.}, gene symbols)
    contained in the K-TSP classifier.}
  \item{classifier}{the classifier obtained by invoking
    \code{\link{SWAP.KTSP.Train}}.}
  \item{DecisionFunc}{ is the function used to generate the final 
   classification prediction by combining the comparisons of the TSPs 
    in the classifier. By default each sample is classified
    according to the class voted by the majority of the TSPs
    (``majority wins''  principle).
    Different decision rules can be also specified using
    alternative functions passed \code{DecisionFunc},
    as described below (see ``details'').}
}

\value{
  This function returns the predicted class for each sample
  in the form of a factor.
}


\details{
  The \code{SWAP.KTSP.Classify} classifies new test samples
  based on a specific decision rule. By default, each sample
  is classified based on the the majority voting rule of 
  the comparisons of TSPs in the classifier.
  Alternative rules can be defined by the user and passed
  to \code{SWAP.KTSP.Classify} using the argument
  \code{DecisionFunc}. A decision function takes
  as its input a logical vector \code{x} corresponding
  to the individual decision of each TSP
  (\code{TRUE} if the first feature in the pair is
  larger then the second, \code{FALSE}
  in the opposite case). The output of the
  \code{DecisionFunction} is a single logical
  value summarizing all votes of the individual TSPs
  (see examples below).
}

 
\seealso{
  \code{\link{SWAP.KTSP.Train}},
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

### Apply the classifier to the TRAINING set using default decision rule
trainingPrediction <- SWAP.KTSP.Classify(matTraining, classifier)

### Resubstitution performance in the TRAINING set
### Define a "positive" test result if needed
table(trainingPrediction, trainingGroup)


### Use an alternative DecideFunction to classify each patient
### Here for instance at least two TSPs must agree
trainingPrediction <- SWAP.KTSP.Classify(matTraining, classifier,
					 DecisionFunc = function(x) sum(x) > 5.5 )

### Contingency table for the TRAINING set
table(trainingPrediction, trainingGroup)


##################################################
### Testing on new data

### Load the example data for the TEST set
data(testingData)

### Show group variable for the TEST set
table(testingGroup)

### Apply the classifier to one sample of the TEST set using default decision rule
testPrediction <- SWAP.KTSP.Classify(matTesting[ , 1, drop=FALSE], classifier)

### Show prediction
testPrediction

### Apply the classifier to the complete the TEST set
### using decision rule defined above (agreement of two TSPs)
testPrediction <- SWAP.KTSP.Classify(matTesting,
  classifier, DecisionFunc = function(x) sum(x) > 5.5)

### Show prediction
head(testPrediction, n=10)

### Contingency table for the TEST set
table(testPrediction, testingGroup)

}

\keyword{KTSP, classification, prediction }

