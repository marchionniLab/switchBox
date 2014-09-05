\name{testingGroup}

\alias{testingGroup}

\docType{data}

\title{Testing set phenotypes}

\description{A factor with two levels
  describing the phenotypes for the testing data
  (Buyse et al cohort,
  (see the \code{\link[mammaPrintData]{mammaPrintData}} package).
}

\usage{
data(testingData)
}

\format{
  The \code{matTesting} factor contains phenotypic information
  for the 307 samples of the testing dataset.
}

\details{
This phenotype factor corresponds to the breast cancer patients' cohort
published by Buyse and colleagues in JNCI (2006).
The gene expression matrix was obtained from
the \code{mammaPrintData} package as described 
by Marchionni and colleagues in BMC Genomics (2013).
}


\seealso{
  \code{\link{trainingGroup}}
}

\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Luigi Marchionni \email{marchion@jhu.edu}
}


\references{
See \link{switchBox} for the references.
}


\examples{

### Load gene expression data for the test set
data(testingData)

### Show the class of the ``testingGroup'' object
class(testingGroup)

### Show group variable
table(testingGroup)

}

\keyword{datasets}
