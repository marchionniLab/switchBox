\name{trainingGroup}

\alias{trainingGroup}

\docType{data}

\title{Training set phenotypes}

\description{A factor with two levels
  describing the phenotypes for the training data
  (Glas et al cohort,
  see the \code{\link[mammaPrintData]{mammaPrintData}} package).
}

\usage{
data(trainingData)
}

\format{
  The \code{trainingGroup} factor contains phenotypic information
  for the 78 samples of the training dataset.
}

\details{
This phenotype factor corresponds to the breast cancer patients' cohort
published by Glas and colleagues in BMC Genomics (2006).
The information was obtained from
the \code{mammaPrintData} package as described 
by Marchionni and colleagues in BMC Genomics (2013).
}


\seealso{
  \code{\link{testingGroup}}
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

### Show the class of the ``trainingGroup'' object
class(trainingGroup)

### Show group variable
table(trainingGroup)

}

\keyword{datasets}
