##################################################
### Functions for backward compatibility with package used
### in Marchionni et al, 2013, BMC Genomics


### #################################################
### #################################################
### KTSP.tiedrank(x) applies rank to the columns of the matrix x.
### It is used to calculate the secondary score.

KTSP.tiedRank <- function(x) {
	return(apply(x, 2, rank))
}



### #################################################
### #################################################
### void CalculateSignedScoreCore(int *situation, int rowLen, double *data1,int colNo1, double *data2, int colNo2,//inputs
### double *score, double *a, double *M, double *k, double *N, double *P,double *Q)//outputs
### Calculates pairwise score between every gene in data1 and every gene in data2. It is signed because
### score(i,j)=P(X_i<X_j|1)-P(X_i<X_j|0) - P(X_i>X_j|1) +P(X_i>X_j|0) + C
### The third and the forth terms are for avoiding the equlities. The C term is proportion to the secondary score
### to avoid the ties.  C = (E(X_j-X_i|1)-E(X_j-X_i|0))/K where K is big enough to make sure that the secondary score
### does not intervene with the primary score.
### a,M,k,N are usable for fisher exact test. M is the size of the class situation == 0 and N is the total number
### of samples. a_ij = #(X_i<X_j|0) k_ij = #(X_i<X_j)
### P_ij = P(X_i<X_j|0) and Q(X_i<X_j|1)

KTSP.CalculateSignedScore <- function(situation , data1, data2) {

	## Check argument 'situation'
	if (! is.numeric(situation)) {
		stop("The argument 'situation' must be a numeric vector with values equal to '0' or '1'. \n")
	}

	## Check argument 'data1'
	if ( ncol(data1) != length(situation)) {
		stop("The number of 'data1' columns must be equal to the length of 'situation'. \n")
	}
	
	## Check argument 'data2'
	if ( ncol(data2) != length(situation)) {
		stop("The number of 'data2' columns must be equal to the length of 'situation'. \n")
	}

	## Check arguments 'data1 and 'data2'
	if ( ncol(data1) != ncol(data2)) {
		stop("The number columns of 'data1' and 'data2' must be equal. \n")
	}

	
	n <- length(situation)
	m1 <- nrow(data1)
	m2 <- nrow(data2)
	
	d <- .C(
		"CalculateSignedScoreCore",
		as.integer(situation),
		as.integer(n),
		as.double(data1), as.integer(m1),
		as.double(data2), as.integer(m2),
		as.double(matrix(0,m1,m2)),
		as.double(matrix(0,m1,m2)),
		as.double(1),
		as.double(matrix(0,m1,m2)),
		as.double(1),
		as.double(matrix(0,m1,m2)),
		as.double(matrix(0,m1,m2))
		)

	retVal<-list(score=(matrix(d[[7]],nrow=m1)),a=(matrix(d[[7]],nrow=m1)),
		     M=(matrix(d[[8]],nrow=1)),k=(matrix(d[[9]],nrow=m1)),
		     N=(matrix(d[[10]],nrow=1)),
		     P=(matrix(d[[11]],nrow=m1)),Q=(matrix(d[[12]],nrow=m1)))
	return(retVal)

}


### #################################################
### #################################################
### KTSP.Train trains the KTSP classifier. 'data' is the matrix of the expression values 
### whose columns represents samples and the rows represents the genes. 
### 'situation' is a vector containing the training labels. Its elements should be one or zero.
### n is the number of top disjoint pairs. 
### If after picking some pairs, only pairs with score left,
### no more pair is chosen even the 'n' has not been reached. 

KTSP.Train <- function(data, situation, n) {

	## Check argument 'situation'
	if (! is.numeric(situation)) {
		stop("The argument 'situation' must be a numeric vector with values equal to '0' or '1'. \n")
	}

	## Check argument 'data'
	if ( ncol(data) != length(situation)) {
		stop("The number of 'data' columns must be equal to the length of 'situation'. \n")
	}

	## Check argument 'n'
	if (! is.numeric(situation)) {
		stop("The argument 'n' must be a single integer specifying the TSP number used in the classifier. \n")
	}

	
	data <- switchBox:::KTSP.tiedRank( data )
	KTSPout <- switchBox:::KTSP.CalculateSignedScore(situation , data, data)
	TSPs <- matrix(0, n, 2)
	TSPscore <- vector(length=n)
	TSPGenes <- matrix("", n, 2)
	geneNames <- rownames(data)
	score <- KTSPout$score

	for (i in 1:n)	{
		TSPs[i,] <- arrayInd(which.max(score), .dim=dim(score))
		TSPscore[i] <- score[TSPs[i, 1], TSPs[i, 2]]/2
		TSPGenes[i, ] <- geneNames[TSPs[i, ]]
		score[TSPs[i , 1] , ] <- 0
		score[TSPs[i , 2] , ] <- 0
		score[ , TSPs[i , 1]] <- 0
		score[ , TSPs[i , 2]] <- 0

		if (norm(score) < 1e-4) {
			TSPs <- TSPs[ 1:i , ]
			TSPscore <- TSPscore[ 1:i ]
			break
		}
	}

	retVal <- list(TSPs=TSPs, score=TSPscore, geneNames=TSPGenes)
	return(retVal)
}



### #################################################
### #################################################
### New function that checks for gene names and enable user define
### rules for tsp combination
### The classifier for the test data. 'data' is the test data. 

KTSP.Classify <- function(data, classifier, combineFunc) {

	## Check argument 'classifier'
	if (missing(classifier)) {
		stop("A valid 'classifier' is needed to classify new data. \n")
	}

	## Check argument 'combineFunc'
	if (! missing(combineFunc)) {
		if (! is.function(combineFunc)) {
			stop("If provided 'combineFunc' must be a function. \n")
		}
	} else {
		## set combineFunc if missing
		combineFunc <- function(x) {
			mean(x) < 0.5
		}
	}

	##retrieve indexes using gene names: preserve the order!
	tspIndexes <- t(apply(classifier$geneNames, 1,
			      FUN <- function(x, y=data) {
				      a <- which(rownames(y) %in% x[1])
				      b <- which(rownames(y) %in% x[2])
				      c(a,b)
			      }))

	##classify using combineFunc
	if (is.vector(data)) {
		x <- data[tspIndexes[, 1]] < data[tspIndexes[, 2]]
		s <- combineFunc(x)
	} else {
		x <- data[tspIndexes[, 1],] < data[tspIndexes[, 2],]
		s <- apply(x, 2, combineFunc)
	}
	
	return(as.integer(s))

}

