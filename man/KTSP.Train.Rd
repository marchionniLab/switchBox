\name{KTSP.Train}

\alias{KTSP.Train}

\title{Funtion for training the K-TSP classifier.}

\description{
  \code{KTSP.Train} trains a K-TSP classifier
  for the specific phenotype of interest.
  The classifiers resulting from using this function can be
  passed to \code{\link{KTSP.Classify}}
  for samples classification.
  This function was used in Marchionni et al, 2013, BMC Genomics,
  and it is maintained only for backward compatibility.
  It has been replaced by \code{\link{SWAP.KTSP.Train}}.
}

\usage{
KTSP.Train(data, situation, n)
}

\arguments{
  \item{data}{the matrix of the values (usually gene expression)
    to be used to train the classifier. The columns represents samples
    and the rows represents the genes.}
  \item{situation}{an integer vector containing the training labels.
    Its elements should be one or zero.}
  \item{n}{The number of disjoint TSP used for classification.
    If before n pairs, the score drops to zero, the TSP with zero score are ignored.}
}

\value{
  The KTSP classifier, a list containing the following elements:
  \item{TSPs}{a matrix containing TSPs indexes.}
  \item{score}{a vector containing TSPs scores.}
  \item{geneNames}{a matrix containing TSPs feature names.}
  It should be passed to KTSP.Classify for classification of test samples.
}


\seealso{
  \code{\link{KTSP.Classify}},
  \code{\link{SWAP.KTSP.Train}},
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

### Turn into a numeric vector with values equal to 0 and 1
trainingGroupNum <- as.numeric(trainingGroup) - 1

### Show group variable for the TRAINING set 
table(trainingGroupNum)


##################################################
### Train a classifier using default filtering function based on the Wilcoxon test
classifier <- KTSP.Train(matTraining, trainingGroupNum, n=8)

### Show the classifier
classifier

}

\keyword{KTSP, classification, prediction}

