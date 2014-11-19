\name{matTraining}

\alias{matTraining}

\docType{data}

\title{Gene expression matrix for training set data}

\description{A numerical matrix containing gene expression matrix
  for 70 genes and 78 breast cancer patients (training set data)
  from the Glas et al cohort
  (see the \code{\link[mammaPrintData]{mammaPrintData}} package).
    
}

\usage{
data(trainingData)
}

\format{
  The \code{matTraining} matrix contains normalized
  expression values for the 70 gene signature (rows)
  across 78 samples (columns).
  Group information (emph{``bad'' versus ``good'' }prognosis)
  is shown in \code{colnames(matTraining)}.
}

\details{
This dataset corresponds to the breast cancer patients' cohort
published by Glas and colleagues in BMC Genomics (2006).
The gene expression matrix was obtained from
the \code{mammaPrintData} package as described 
by Marchionni and colleagues in BMC Genomics (2013).
}


\seealso{
  \code{\link{matTesting}}
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

### Show the class of the ``matTraining'' object
class(matTraining)

### Show the dimentions of the ``matTraining'' matrix
dim(matTraining)

### Show the first  10 sample names of the ``matTraining'' matrix
head(colnames(matTraining), n=10)

}

\keyword{datasets}
