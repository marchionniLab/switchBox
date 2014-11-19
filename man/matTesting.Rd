\name{matTesting}

\alias{matTesting}

\docType{data}

\title{Gene expression matrix for test set data}

\description{A numerical matrix containing gene expression matrix
  for 70 genes and 307 breast cancer patients (test set data)
  from the Buyse et al cohort
  (see the \code{\link[mammaPrintData]{mammaPrintData}} package).
  
}

\usage{
data(testingData)
}

\format{
  The \code{matTesting} matrix contains normalized
  expression values for the 70 gene signature (rows)
  across 307 samples (columns).
  Group information (emph{``bad'' versus ``good'' }prognosis)
  is shown in \code{colnames(matTesting)}. 
}

\details{
This dataset corresponds to the breast cancer patients' cohort
published by Buyse and colleagues in JNCI (2006).
The gene expression matrix was obtained from
the \code{mammaPrintData} package as described 
by Marchionni and colleagues in BMC Genomics (2013).
}


\seealso{
  \code{\link{matTraining}}
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

### Show the class of the ``matTesting'' object
class(matTesting)

### Show the dimentions of the ``matTesting'' matrix
dim(matTesting)

### Show the first  10 sample names of the ``matTest'' matrix
head(colnames(matTesting), n=10)
testingGroup[1:10]

}

\keyword{datasets}
