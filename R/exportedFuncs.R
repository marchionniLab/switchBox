##################################################
##################################################
### An function used to calculate the pairwise scores.
### We can restrict the pairs for which we calculate the scores.

SWAP.CalculateSignedScore <- function(inputMat, phenoGroup, 
				      FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs, ...) {
	
	## Check inputMat conformity
	SWAP.Check.Input(phenoGroup = phenoGroup, inputMat = inputMat,
			 RestrictedPairs = RestrictedPairs)

	## Filter input and format the data
	formattedInput <- SWAP.Format.Input(phenoGroup =phenoGroup,
					    inputMat = inputMat, RestrictedPairs = RestrictedPairs,
					    FilterFunc = FilterFunc, ...)

	## Compute the "unrestricted" classifier
	if ( missing(RestrictedPairs) ) {
		## Obtain the TSP scores
		out <- calculateSignedScore(phenoGroup = phenoGroup ,
					    inputMat1 = formattedInput$inputMat,
					    inputMat2 = formattedInput$inputMat)
	} else {
		## Then calculate the scores for these RestrictedPairs
		out <- calculateSignedScore(phenoGroup = phenoGroup,
					    inputMat1 = formattedInput$inputMat,
					    inputMat2 = formattedInput$inputMat,
					    RestrictedPairs = formattedInput$FilteredPairs )
	}

	### Return output
	return(out)
}



##################################################
##################################################
### SWAP.KTSP.Train trains a KTSP classifier.
### `inputMat' is a numerical matrix with columns representing samples
### and rows representing features (e.g. genes).
### `phenoGroup' is a factor containing the training labels.
### `krange' is the range of top disjoint pairs used in the KTSP classifier
### If after picking some pairs, only pairs with negligible score remains,
### then no more pair is chosen even the 'n' has not been reached.

SWAP.KTSP.Train <- function(inputMat, phenoGroup, krange = c(3, 5, 7:10),
			    FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs, ...) {
	
	## Check inputMat conformity
	SWAP.Check.Input(phenoGroup = phenoGroup, inputMat = inputMat,
			 RestrictedPairs = RestrictedPairs)

	## Filter input and format the data
	formattedInput <- SWAP.Format.Input(phenoGroup, inputMat = inputMat,
					    RestrictedPairs = RestrictedPairs, FilterFunc = FilterFunc, ...)

	## Compute the "unrestricted" classifier
	if ( missing(RestrictedPairs) ) {
		classifier <- SWAP.KTSP.Train.Plain(inputMat = formattedInput$inputMat,
						    phenoGroup = phenoGroup, maxK = max(krange) )
	} else {
		## Develop restricted classifier
		classifier <- SWAP.KTSP.Train.Restricted(inputMat = formattedInput$inputMat,
							 phenoGroup = phenoGroup, maxK = max(krange),
							 RestrictedPairs = formattedInput$FilteredPairs)
	}
	
	## Use min T-test to select k
	cat("Selecting K...\n")
	kmin <- bestKByTtest(classifier = classifier, inputMat = inputMat,
			     phenoGroup = phenoGroup, krange = krange)

	## Best K
	kstar <- min(krange[kmin], nrow(classifier$TSPs))
	if ( kstar != krange[kmin] ) {
		cat(paste("The required range of k is not available!\n",
			  "The minimum number of available TSP (", kstar,
			  ") will be used instead.\n", sep=""))
	} else {
		cat(paste(kstar, "TSP will be used to build the final classifier.\n"))
	}

	## Prepare final output
	out <- list(name = sprintf("%dTSPs", kstar),
		    TSPs = classifier$TSPs[ 1:kstar , ],
		    score = classifier$score[ 1:kstar ],
		    labels = classifier$labels)

	## Return out
	return(out)
}



##################################################
##################################################
### It calculates the KTSP statistics : \sum_{k} I(X_i_k <X_j_k) - I(X_j_k <X_i_k)
### statistics = SWAP.KTSP.Statitsics(datTraining, classifiers)

SWAP.KTSP.Statistics <- function(inputMat, classifier, CombineFunc) {

	## Checking inputMat
	if ( !is.matrix(inputMat) || !is.numeric(inputMat) ) {
		stop("inputMat must be a numeric matrix.")
	}

	## Check if this a one sample classification situation
	if ( ncol(inputMat) == 1 ) {
		inputMat <- as.matrix(inputMat)
	}

	## Comparisons
	comparisons <- t(inputMat[classifier$TSPs[ , 1 ] , ] > inputMat[classifier$TSPs[ , 2 ] , ])
	## Add row and column names
	colnames(comparisons) <- mapply(x = classifier$TSPs[ , 1 ] ,
					y = classifier$TSPs [ , 2 ] ,
					FUN = function(x,y) sprintf("%s>%s", x, y ) )
	rownames(comparisons) <- colnames(inputMat)
	
	## Compute feature by samples switching
	stats1 <- inputMat[classifier$TSPs[ , 1 ] , ] > inputMat[classifier$TSPs[ , 2 ] , ]
	stats2 <- inputMat[classifier$TSPs[ , 1 ] , ] < inputMat[classifier$TSPs[ , 2 ] , ]
	
	## Set the combine function if missing
	if (missing(CombineFunc)) {
		## Compute standard K statistics
		stats1 <- apply(as.matrix(stats1), 2 , sum)
		stats2 <- apply(as.matrix(stats2), 2 , sum)
		KTSPstat <- stats1 - stats2
		## Add names
		names(KTSPstat) <- colnames(inputMat)
	} else {
		## Apply combined function
		x <- stats1
		KTSPstat <- apply(x, 2, CombineFunc)
	}
	
	## ## Make plot if needed
	## if (show)
	## plot the figure like the figure in the paper.
	##

	## Return statistics
	out <- list(statistics = KTSPstat , comparisons = comparisons)
	return(out)

}

	
##################################################
##################################################
### The classifier for the test inputMat. 'inputMat' is the test inputMat. 
### just threshold the outcome of KTSP Statistics 
### This will be handled by a function inthe future

### CombineFunc must return a logical indicator
### TRUE for Good (2nd level in phenoGroup)
### FALSE for Bad (1st level in phenoGroup)

SWAP.KTSP.Classify <- function(inputMat, classifier, DecisionFunc) {

	## Get the group labels
	mylabels <- classifier$labels

	## If CombineFunction is missing use default (majority wins)
	if (missing(DecisionFunc)) {
		## Compute the KTSP statistics
		ktspStat <- SWAP.KTSP.Statistics(inputMat, classifier)$statistics > 0
	} else {
		## Use CombineFunc to compute the KTSP statistics
		ktspStat <- SWAP.KTSP.Statistics(inputMat, classifier, DecisionFunc)$statistics
	}

	## Prepare and return the results
	out <- factor(ifelse(ktspStat, mylabels[[2]], mylabels[[1]]), levels=mylabels)
	return(out)

}



##################################################
##################################################
### A function to filter the genes, other functions taking similar 
### arguments will be passed to SWAP.KTSP.Train() in the future.
### Here the default  filtering is based on the Wilcoxon rank-sum test:
### 50 top up-regulated and 50 top down-regulated genes

SWAP.Filter.Wilcoxon <- function(phenoGroup, inputMat, featureNo = 100,
				 UpDown = TRUE) {

	## Check inputMat conformity
	SWAP.Check.Input(phenoGroup, inputMat)

	## Get ranks and perform test
	tiedData <- apply(inputMat, 2 , rank)
	tiedDataP <- t(apply(tiedData, 1 , rank))
	n <- sum(phenoGroup == levels(phenoGroup)[1])
	m <- sum(phenoGroup == levels(phenoGroup)[2])
	sumzeros <- apply(tiedDataP[ , which(phenoGroup == levels(phenoGroup)[1])], 1 , sum)
	windex <- (sumzeros - n*(n+m+1)/2) / sqrt(n*m*(n+m+1)/12)
	
	## Retrieve the features
	if (UpDown) {
		s <- order(windex, decreasing = TRUE)
		lens <- length(s)
		## UP
		featuresIndexUp <- s[1:min(c(round(featureNo/2) , lens))]
		## DOWN
		featuresIndexDown <- s[ max(c(lens - round(featureNo/2) , 1)) : lens]
		## Gene indexes
		featuresIndex <- unique(c(featuresIndexUp , featuresIndexDown))
	} else {
	## The top features
		s <- order(abs(windex) , decreasing=TRUE)
		featuresIndex <- s[1:min(c(featureNo , length(s)))]
	}

	## Return results as ronames of inputMat
	out <- rownames(inputMat)[featuresIndex]
	return(out)

}

