##################################################
##################################################
### An internal function to check the inputs.

SWAP.Check.Input <- function(phenoGroup, inputMat, RestrictedPairs) {

	## Checking "phenoGroup"
	if ( (!is.factor(phenoGroup)) || (length(levels(phenoGroup)) != 2) ) {
		stop("Check input:  'phenoGroup' must be a factor with exactly two levels.")
	}

	## Checking "inputMat"
	if ( ! missing(inputMat) ) {
		if ( is.null(rownames(inputMat)) || (length(phenoGroup) != ncol(inputMat)) ||
		    (!is.matrix(inputMat)) || (!is.numeric(inputMat)) ) {
			stop(paste("Check input: 'inputMat' must be a numeric matrix with rownames",
				   "and a number of columns equal to the length of phenoGroup."))
		}
	}

	## Checking "RestrictedPairs", if available
	if ( ! missing(RestrictedPairs) ) {
		if ( (!is.matrix(RestrictedPairs)) || (!is.character(RestrictedPairs)) ||
		    (ncol(RestrictedPairs) != 2) ) {
			stop(paste("Check input: 'RestrictedPairs' must be a 2 columns matrix",
				   "containing strings matching rownames in inputMat."))
		}

		## Checking "RestrictedPairs" and inputMat consistency
		if ( all( ! RestrictedPairs[ , 1] %in% rownames(inputMat) ) ||
		    all( ! RestrictedPairs[ , 2] %in% rownames(inputMat) ) ) {
			stop("None of the 'RestrictedPairs' matches 'rownames(inputMat)'!" )
		}
	}
}



##################################################
##################################################
### SWAP.format.Input prepares input for lower level training functions

SWAP.Format.Input <- function(phenoGroup, inputMat, FilterFunc, RestrictedPairs, ...) {

	## Apply filter function
	if ( is.null(FilterFunc) ) {
		Filter <- rownames(inputMat)
		cat("No feature filtering procedure will be used...\n")
	} else {
		cat("Applying filtering function to 'inputMat'...\n")
		Filter <- FilterFunc(phenoGroup, inputMat, ...)
	}

	## Prepare inputMat using filtered features
	inputMat <- inputMat[Filter , ]

	## Prepare data to compute the "unrestricted" classifier
	if ( missing(RestrictedPairs) ) {
		## Check if there are at least 4 features left to build the TSPS
		if ( length(rownames(inputMat)) < 4 ) {
			stop("Not enough features left after feature filtering!")
		} else {
			## Return input matrix
			out <- list(inputMat = inputMat)
		}
	} else {
	## Filter the inputMat for computing the classifier with "RestrictedPairs"
	## Check the intersection between filtered features and restriceted pairs
		moreFilter <- intersect(as.vector(RestrictedPairs), rownames(inputMat))
		## CHECK HERE
		## print(paste("Check1:", length(moreFilter)))
		if ( length(moreFilter) < 4 ) {
			stop("Not enough features left after filtering and restiction of pairs!")
		} else {
			cat("Restricting the analysis to the provided candidate TSPs\n")
			## Subset to available feature set
			## CHECK HERE
			## print(paste("Check2:", dim(inputMat)))
			inputMat <- inputMat[ moreFilter , ]
			## CHECK HERE
			## print(paste("Check3:", dim(inputMat)))
			## Get the features
			features <- rownames(inputMat)
			## First find RestrictedPairs for which both features are available
			keepPairs <- which( (RestrictedPairs[,1] %in% features) &
					   (RestrictedPairs[,2] %in% features) )
			## Is any of the RestrictedPairs left to go ahead?
			if ( length(keepPairs) < 2 ) {
				stop("Not enough features left after filtering and restiction of pairs!")
			} else {
				## Retain the available pairs
				FilteredPairs <- RestrictedPairs[ keepPairs , ]
			}
			## Apply filter and return input matrix
			out <- list(inputMat = inputMat, FilteredPairs = FilteredPairs)
		}
	}
	## Return inout data after filtering and/or pair restriction
	return(out)
}



##################################################
##################################################
### The internal  function computing the score using the C code

calculateSignedScore <- function(phenoGroup , inputMat1, inputMat2,
				      RestrictedPairs) {
	
	## Setting variables for C procedure
	n <- length(phenoGroup)
	m1 <- nrow(inputMat1)
	m2 <- nrow(inputMat2)

	##  Get ranks and ties
	datatied <- apply(rbind(inputMat1, inputMat2), 2, rank)
	data1tied <- datatied[1:m1 , ]
	data2tied <- datatied[(m1+1) : (m1+m2) , ]

	### Get labels and create a vector from the group factor
	labels <- levels(phenoGroup)
	situationV <- as.vector(phenoGroup)

	## Launch the C procedure to compute the signed score
	## Non-restricted version
	if ( missing(RestrictedPairs) ) {

		## How many fearures are used
		cat(paste("Computing scores for ", m1, " features.\n",
			  "This will require enough memory for ",
			  formatC(m1*(m1-1)/2), " pairs.\n", sep=""))
		
		## Launch the C code
		d <- .C(
			"CalculateSignedScoreCore",
			as.integer(as.numeric(situationV == labels[1])),
			as.integer(n),
			as.double(data1tied), as.integer(m1),
			as.double(data2tied), as.integer(m2),
			as.double(matrix(0, m1, m2)),
			as.double(matrix(0, m1, m2)),
			as.double(1),
			as.double(matrix(0, m1, m2)),
			as.double(1),
			as.double(matrix(0, m1, m2)),
			as.double(matrix(0, m1, m2))
			)

		## Extract components from the output
		score <- (matrix(d[[7]], nrow = m1))
		P <- (matrix(d[[12]], nrow = m1))
		Q <- (matrix(d[[13]], nrow = m1))

		## Process feature names
		names1 <- rownames(inputMat1)
		names2 <- rownames(inputMat2)

		## Assign row and columns names to scores
		rownames(score) <- names1
		colnames(score) <- names2

		## Assign row and columns names to P
		rownames(P) <- names1
		colnames(P) <- names2

		## Assign row and columns names to Q
		rownames(Q) <- names1
		colnames(Q) <- names2	

		## Prepare output object
		retVal <- list(score=score/2, P=P, Q=Q, labels=labels)
		
	## Launch the C procedure to compute the signed score
	## Restricted case
	} else {

		## Prepare the inputMat
		pairsno <- nrow(RestrictedPairs)
		pairsind <- matrix(0, pairsno, 2)
		pairsind[ , 1] <- match(RestrictedPairs[ , 1], rownames(inputMat1))
		pairsind[ , 2] <- match(RestrictedPairs[ , 2], rownames(inputMat2))

		## How many fearures are used
		cat(paste("Computing scores for ", pairsno, " available restricted pairs.\n",
			  "This will require enough memory for ", pairsno, " pairs.\n", sep=""))

		## Run the C code
		d <- .C(
			"CalculateSignedScoreRestrictedPairsCore",
			as.integer(as.numeric(situationV == labels[1])),
			as.integer(n),
			as.double(data1tied), as.integer(m1),
			as.double(data2tied), as.integer(m2),
			as.integer(pairsind[ , 1]-1),
			as.integer(pairsind [ , 2]-1),
			as.integer(pairsno),
			as.double(matrix(0, pairsno, 1)),
			as.double(matrix(0, pairsno, 1)),
			as.double(matrix(0, pairsno, 1))
			)

		## Extract scores, P, and Q from the output
		score <- d[[10]]
		names(score) <- sprintf("%s,%s", RestrictedPairs[,1], RestrictedPairs[,2])
		P <- d[[11]]
		Q <- d[[12]]

		## Prepare the output object
		retVal <- list(score=score/2, P=P, Q=Q, labels=labels)
	}

	## Return output
	return(retVal)

}



#################################################
#################################################
### An internal function for finding kTSP without restricting RestrictedPairs

SWAP.KTSP.Train.Plain <- function(inputMat, phenoGroup, maxK) {

	## Get class lables
	labels <- levels(phenoGroup)
	
	## Obtain the TSP scores
	KTSPout <- calculateSignedScore(phenoGroup , inputMat1 = inputMat,
					     inputMat2 = inputMat)

	## Process the TSP score
	TSPscore <- vector(length = maxK)
	score <- KTSPout$score
	absscore <- abs(score)
	TSPsInd <- matrix(0, maxK , 2)

	## Select the best TSPs
	for (i in 1:maxK )	{
		nextPair <- arrayInd(which.max(absscore), .dim = dim(score));

		## If the score of the pairs are negligible do not go down the list
		if ( absscore[ nextPair[1] , nextPair[2] ] < 1e-4 ) {
			TSPsInd <- TSPsInd[ 1:i , ]
			TSPscore <- TSPscore[ 1:i ]
			break
		}	

		## If the score of the pairs is negative
		if ( score[ nextPair[1] , nextPair[2] ] < 0 ) {
			TSPsInd[ i , ] <- nextPair
		} else {
		## If it is positive and not negligible
			TSPsInd[ i , ] <- c( nextPair[2] , nextPair[1] )
		}

		### Calculate the TSP score
		TSPscore[i] <- absscore[ TSPsInd[ i , 1 ] , TSPsInd[ i , 2 ] ]

		## Assign 0 value
		absscore[ TSPsInd[ i , 1:2 ] , ] <- 0
		absscore[ , TSPsInd[ i , 1:2 ] ] <- 0
	}

	## PREPARE TSP
	features <- rownames(inputMat)
	TSPs <- cbind( features[TSPsInd[ , 1 ]] , features[TSPsInd[ , 2]] )

	## Prepare output object as a list
	classifier <- list()

	## Prepare the TSPs component
	maxTSPs <- nrow(TSPs)
	classifier$name <- sprintf('%dTSPs', maxTSPs)
	if ( maxTSPs == 1) {
		classifier$TSPs <- matrix(TSPs[1,], ncol=2, nrow=1)
	} else {
		classifier$TSPs <- TSPs[1:maxTSPs , ]
	}

	## Prepare the score and labels components
	classifier$score <- TSPscore[1:maxTSPs]
	classifier$labels <- labels

	## Return classifier object
	return(classifier)

}


#################################################
#################################################
### An internal function for finding kTSP using restricted pairs

SWAP.KTSP.Train.Restricted <- function(inputMat, phenoGroup, maxK, RestrictedPairs) {

	## Get the lables
	labels <- levels(phenoGroup)

	## Then calculate the scores for these RestrictedPairs
	scores <- calculateSignedScore(phenoGroup,
					    inputMat1 = inputMat,
					    inputMat2 = inputMat,
					    RestrictedPairs = RestrictedPairs  )

	## Process scores
	absOrder <- order(abs(scores$score), decreasing = TRUE)

	## Disjoint TSPs
	pairs <- SWAP.KTSP.Pairs.Disjoint(
		pairs=RestrictedPairs[ absOrder , ],
		scores$score[absOrder],
		min(maxK, length(scores$score)))

	## Prepare output object as a list
	classifier = list();

	## Prepare the TSPs component
	maxTSPs <- length(pairs$scores)
	classifier$TSPs <- pairs$pairs[ 1:maxTSPs , ]

	## Swap gene order in the pair for pairs with negative score
	posPairs <- which(pairs$scores > 0)
	classifier$TSPs[posPairs , ] <- classifier$TSPs[posPairs , 2:1]

	### Flippling the negative scores
	classifier$score <- pairs$scores[ 1:maxTSPs , ]
	classifier$score <- abs(classifier$score)
	
	## Prepare the score and labels components
	classifier$name <- sprintf('%dTSPs', maxTSPs)
	classifier$labels <- labels

	## Remove names (to make output consistent with unrestricted version)
	dimnames(classifier$TSPs) <- NULL
	names(classifier$score) <- NULL

	### Return
	return(classifier)

}



##################################################
##################################################
### An internal function for making k-disjoint-TSP 

SWAP.KTSP.Pairs.Disjoint <- function(pairs, scores, k) {

	### Pick the first pair
	kTSPs <- pairs[1 , ]

	## Get the score for the first pair
	kTSPscores <- scores[1]

	## Set counter "n" equal to 0
	n = 1
	for (i in 2:nrow(pairs)) {
		if ( pairs[i, 1] %in% kTSPs == FALSE &&
		    pairs[i, 2] %in% kTSPs == FALSE ) {
			kTSPs <- rbind(kTSPs, pairs[i , ])
			kTSPscores <- rbind(kTSPscores, scores[i])
			## Update counter "n"
			n = n + 1
			if ( n >= k ) {
				break
			}
		}	
	}

	## Create and return result object
	out <- list(pairs=kTSPs, scores=kTSPscores)
	return(out)

}



##################################################
##################################################
### Function to pick the best K by min t-test

bestKByTtest <-  function(classifier, inputMat, phenoGroup, krange) {

	## Define max TSPs
	maxTSPs <- nrow(classifier$TSPs)

	##Handle levels and lables
	s0 <- which( phenoGroup == levels(phenoGroup)[1] )
	s1 <- which( phenoGroup == levels(phenoGroup)[2] )

	## Here is minimum "t-Test"
	mintt <- -Inf
	classifiertest <- classifier

	## Do the test for all k ranges
	for ( k in 1:length(krange) ) {
		classifiertest$TSPs <- classifier$TSP[ 1:min(krange[k] , maxTSPs) , ]
		stat0 <- SWAP.KTSP.Statistics(inputMat[ , s0], classifiertest)$statistics
		stat1 <- SWAP.KTSP.Statistics(inputMat[ , s1], classifiertest)$statistics

		##current t-test
		tt <- ( abs(mean(stat0) - mean(stat1)) /
		       sqrt(var(stat1) + var(stat0) + 0.000000001) )

		## If there is a tie score, we pick the classifier which shows up first
		if (abs(mintt - tt) > 0.0000001 && mintt < tt) {
			mintt <- tt
			kmin <- k
		}
	}

	### Return kmin
	return(kmin)

}
