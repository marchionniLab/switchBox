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
        if ( ! (missing(RestrictedPairs) || is.null(RestrictedPairs)) ) {
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

SWAP.Format.Input <- function(phenoGroup, inputMat, FilterFunc, RestrictedPairs, verbose = FALSE, ...) {

        ## Apply filter function
        if ( is.null(FilterFunc) ) {
                Filter <- rownames(inputMat)
                if(verbose) cat("No feature filtering procedure will be used...\n")
        } else {
                if(verbose) cat("Applying filtering function to 'inputMat'...\n")
                Filter <- FilterFunc(phenoGroup, inputMat, ...)
        }

        ## Prepare inputMat using filtered features
        inputMat <- inputMat[Filter , ]

        ## Prepare data to compute the "unrestricted" classifier
        if ( missing(RestrictedPairs) || is.null(RestrictedPairs) ) {
                ## Check if there are at least 4 features left to build the TSPS
                if ( length(rownames(inputMat)) < 4 ) {
                        stop("Not enough features left after feature filtering!")
                } else {
                        ## Return input matrix
                        out <- list(inputMat = inputMat, FilteredPairs = NULL)
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
                        if(verbose) cat("Restricting the analysis to the provided candidate TSPs\n")
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

calculateScores <- function(inputMat, phenoGroup, classes = NULL, FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = NULL, 
    handleTies = FALSE, verbose = FALSE, score_fn = signedTSPScores, score_opts = list(), ...){

    #print(score_fn)

    ## Check inputMat conformity
    SWAP.Check.Input(phenoGroup = phenoGroup, inputMat = inputMat,
                     RestrictedPairs = RestrictedPairs)

    ## Filter input and format the data
    formattedInput <- SWAP.Format.Input(phenoGroup =phenoGroup,
                                        inputMat = inputMat, RestrictedPairs = RestrictedPairs,
                                        FilterFunc = FilterFunc, verbose = verbose, ...)

    scores = score_fn(phenoGroup=phenoGroup, inputMat1=formattedInput$inputMat, inputMat2=formattedInput$inputMat, classes,
        RestrictedPairs=formattedInput$FilteredPairs, handleTies=handleTies, verbose=verbose, score_opts = score_opts)

    if(! exists("signed", scores))
        scores$signed = FALSE

    return(scores)

}

signedTSPScores = function(phenoGroup, inputMat1, inputMat2 = NULL, classes = NULL, RestrictedPairs = NULL, 
    handleTies = FALSE, verbose = FALSE, score_opts=list()){

    if(is.null(inputMat2))
        inputMat2 = inputMat1

    scores = calculateSignedScore(phenoGroup=phenoGroup, 
        inputMat1=inputMat1, inputMat2=inputMat2,
        RestrictedPairs=RestrictedPairs, handleTies=handleTies, verbose=verbose, classes = classes)

    # for a non-restricted classifiers, scores are in a matrix; convert to a vector
    if(is.matrix(scores$score)){
        scores$score = score_matrix_to_vector(scores$score)

        # if no tie handling is done, then set tieVote = 0
        if(exists("tieVote", scores))
            scores$tieVote = score_matrix_to_vector(scores$tieVote)
        else
            scores$tieVote = rep(0, length(scores$score))
    }else{
        # scores are already in a vector
        if(! exists("tieVote", scores))
            scores$tieVote = rep(0, length(scores$score))
    }

    scores$signed = TRUE

    return(scores)

}

basicTSPScores = function(phenoGroup, inputMat1, inputMat2 = NULL, classes = NULL, RestrictedPairs = NULL, 
    handleTies = FALSE, verbose = FALSE, score_opts=list()){

    if(is.null(inputMat2))
        inputMat2 = inputMat1

    if(is.null(classes))
        levels = levels(phenoGroup)
    else
        levels = rev(classes)

    scores = list(labels=levels)

    if(is.null(RestrictedPairs)){

        score = t( sapply(rownames(inputMat1), function(xi) sapply(rownames(inputMat2), function(xj) 
            calculate_tsp_instance(inputMat1[xi, ], inputMat2[xj, ], phenoGroup, levels) ) ) )

        rownames(score) = rownames(inputMat1)
        colnames(score) = rownames(inputMat2)

        scores$score = score_matrix_to_vector(score)

    }else{

        score = apply(RestrictedPairs, 1, function(p) calculate_tsp_instance(inputMat1[p[1], ], inputMat2[p[2], ], phenoGroup, levels) )
        names(score) = apply(RestrictedPairs, 1, function(p) sprintf("%s,%s", p[1], p[2]))

        scores$score = score

    }

    if(verbose && handleTies) cat("Tie handling not available for basic TSP score\n")

    scores$tieVote = rep(0, length(scores$score))

    return(scores)

}

calculateSignedScore <- function(phenoGroup , inputMat1, inputMat2, 
                                      RestrictedPairs=NULL, handleTies=FALSE, verbose = FALSE, classes = NULL) {

        ## Setting variables for C procedure
        n <- length(phenoGroup)
        m1 <- nrow(inputMat1)
        m2 <- nrow(inputMat2)

        ##  Get ranks and ties
        datatied <- apply(rbind(inputMat1, inputMat2), 2, rank)
        data1tied <- datatied[1:m1 , ]
        data2tied <- datatied[(m1+1) : (m1+m2) , ]

        ### Get labels and create a vector from the group factor
        if(is.null(classes)){
            labels <- levels(phenoGroup)
        }else{
            labels <- rev(classes)
        }
        
        situationV <- as.vector(phenoGroup)

        if(handleTies){
                if ( missing(RestrictedPairs) || is.null(RestrictedPairs) ) {

                        ## How many features are used
                        if(verbose) cat(paste("Computing scores for ", m1, " features.\n",
                                  "This will require enough memory for ",
                                  formatC(m1*(m1-1)/2), " pairs.\n", sep=""))

                        ## Launch the C code
                        d <- .C(
                                "CalculateSignedScoreCoreTieHandler",
                                as.integer(as.numeric(situationV == labels[1])),
                                as.integer(n),
                                as.double(data1tied), as.integer(m1),
                                as.double(data2tied), as.integer(m2),
                                as.double(matrix(0, m1, m2)),
                                as.double(matrix(0, m1, m2)),
                                as.double(matrix(0, m1, m2)),
                                as.double(matrix(0, m1, m2))
                                )

                        ## Extract components from the output
                        score <- (matrix(d[[7]], nrow = m1))
                        P <- (matrix(d[[8]], nrow = m1))
                        Q <- (matrix(d[[9]], nrow = m1))
                        tieVote <- (matrix(d[[10]], nrow = m1))

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

                        rownames(tieVote) <- names1
                        colnames(tieVote) <- names2


                        ## Prepare output object
                        retVal <- list(score=score, labels=labels, tieVote=tieVote) #P=P, Q=Q, 

                ## Launch the C procedure to compute the signed score
                ## Restricted case
                }else{

                        ## Prepare the inputMat
                        pairsno <- nrow(RestrictedPairs)
                        pairsind <- matrix(0, pairsno, 2)
                        pairsind[ , 1] <- match(RestrictedPairs[ , 1], rownames(inputMat1))
                        pairsind[ , 2] <- match(RestrictedPairs[ , 2], rownames(inputMat2))

                        ## How many fearures are used
                        if(verbose) cat(paste("Computing scores for ", pairsno, " available restricted pairs.\n",
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
                                as.double(matrix(0, pairsno, 1)),
                                as.double(matrix(0, pairsno, 1))
                                )

                        ## Extract scores, P, and Q from the output
                        score <- d[[10]]
                        names(score) <- sprintf("%s,%s", RestrictedPairs[,1], RestrictedPairs[,2])
                        P <- d[[11]]
                        Q <- d[[12]]
                        tieVote<- d[[13]]

                        ## Prepare the output object
                        retVal <- list(score=score/2, labels=labels,tieVote=tieVote) #P=P, Q=Q, 
                }
        }else{
                ## Launch the C procedure to compute the signed score
                ## Non-restricted version
                if ( missing(RestrictedPairs) || is.null(RestrictedPairs) ) {

                        ## How many fearures are used
                        if(verbose) cat(paste("Computing scores for ", m1, " features.\n",
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
                        retVal <- list(score=score/2, labels=labels) #P=P, Q=Q, 

                ## Launch the C procedure to compute the signed score
                ## Restricted case
                } else {

                        ## Prepare the inputMat
                        pairsno <- nrow(RestrictedPairs)
                        pairsind <- matrix(0, pairsno, 2)
                        pairsind[ , 1] <- match(RestrictedPairs[ , 1], rownames(inputMat1))
                        pairsind[ , 2] <- match(RestrictedPairs[ , 2], rownames(inputMat2))

                        ## How many fearures are used
                        if(verbose) cat(paste("Computing scores for ", pairsno, " available restricted pairs.\n",
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
                        retVal <- list(score=score/2, labels=labels) #P=P, Q=Q, 
                }
        }
        ## Return output
        return(retVal)

}



#################################################
#################################################
### An internal function for finding kTSP without restricting RestrictedPairs

SWAP.KTSP.Train.Plain <- function(inputMat, phenoGroup, maxK, verbose = FALSE) {

        ## Get class lables
        labels <- levels(phenoGroup)

        ## Obtain the TSP scores
        KTSPout <- calculateSignedScore(phenoGroup , inputMat1 = inputMat,
                                             inputMat2 = inputMat, verbose = verbose)

        ## Process the TSP score
        TSPscore <- vector(length = maxK)
        score <- KTSPout$score
        absscore <- abs(score)
        TSPsInd <- matrix(0, maxK , 2)

        ## Select the best TSPs
        for (i in 1:maxK )      {
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

SWAP.KTSP.Train.Restricted <- function(inputMat, phenoGroup, maxK, RestrictedPairs, verbose = FALSE) {

        ## Get the lables
        labels <- levels(phenoGroup)

        ## Then calculate the scores for these RestrictedPairs
        scores <- calculateSignedScore(phenoGroup,
                                            inputMat1 = inputMat,
                                            inputMat2 = inputMat,
                                            RestrictedPairs = RestrictedPairs, verbose = verbose )

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

SWAP.KTSP.Pairs.Disjoint <- function(pairs, scores, k, tieVote) {

        if(missing(tieVote)){
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
        }else{
                        ### Pick the first pair
                kTSPs <- pairs[1 , ]

                ## Get the score for the first pair
                kTSPscores <- scores[1]
                kTSPtieVotes <- tieVote[1]

                ## Set counter "n" equal to 0
                n = 1
                for (i in 2:nrow(pairs)) {
                        if ( pairs[i, 1] %in% kTSPs == FALSE &&
                                pairs[i, 2] %in% kTSPs == FALSE ) {
                                kTSPs <- rbind(kTSPs, pairs[i , ])
                                kTSPscores <- rbind(kTSPscores, scores[i])
                                kTSPtieVotes  <- rbind(kTSPtieVotes, tieVote [i])
                                ## Update counter "n"
                                n = n + 1
                                if ( n >= k ) {
                                        break
                                }
                        }
                }

                ## Create and return result object
                out <- list(pairs=kTSPs, scores=kTSPscores, tieVotes = kTSPtieVotes)

        }
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

KbyTtest <-  function(inputMat, phenoGroup, scoreTable, classes, krange, k_opts=list()){

    ## Define max TSPs
    maxTSPs <- nrow(scoreTable)

    ##Handle levels and lables
    s0 <- which( phenoGroup == classes[2] )
    s1 <- which( phenoGroup == classes[1] )

    ## Here is minimum "t-Test"
    mintt <- -Inf
    classifiertest <- list(labels = rev(classes))

    ## Do the test for all k ranges
    for ( k in 1:length(krange) ) {
            
            classifiertest$TSPs <- scoreTable[1:min(krange[k], maxTSPs), c("gene1", "gene2")]
            classifiertest$TSPs <- as.matrix(classifiertest$TSPs, nrow(classifiertest$TSPs), 2)

            classifiertest$tieVote = scoreTable$tieVote[ 1:min(krange[k], maxTSPs) ]

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
    return(1:krange[kmin])

}

#################################################
#################################################
### An internal function for finding kTSP without restricting RestrictedPairs with Handlingties

SWAP.KTSP.Train.Plain.Ties <- function(inputMat, phenoGroup, maxK, verbose = FALSE) {


        ## Get class labels
        labels <- levels(phenoGroup)

        ## If we reverse the pairs to make the score positive, the tieVote 1 and 2 become 2 and 1 respectively
        tieVotechange <- c(0,2,1)

        ## Obtain the TSP scores
        KTSPout <- calculateSignedScore(phenoGroup , inputMat1 = inputMat,
                                             inputMat2 = inputMat, handleTies = TRUE, verbose = verbose)

        ## Process the TSP score
        TSPscore <- vector(length = maxK)
        score <- KTSPout$score
        tieVote <- KTSPout$tieVote
        absscore <- abs(score)
        TSPsInd <- matrix(0, maxK , 2)
        TSPtieVote <- vector(length = maxK)

        ## Select the best TSPs
        for (i in 1:maxK )      {
                nextPair <- arrayInd(which.max(absscore), .dim = dim(score));

                #cat(rownames(inputMat)[nextPair], "\t")
                
                ## If the score of the pairs are negligible do not go down the list
                if ( absscore[ nextPair[1] , nextPair[2] ] < 1e-4 ) {
                        TSPsInd <- TSPsInd[ 1:i , ]
                        TSPscore <- TSPscore[ 1:i ]
                        TSPtieVote <- TSPtieVote[ 1:i ]
                        break
                }

                ## If the score of the pairs is negative
                if ( score[ nextPair[1] , nextPair[2] ] < 0 ) {
                        TSPsInd[ i , ] <- nextPair
                        TSPtieVote[i] <- tieVote[ TSPsInd[ i , 1 ] , TSPsInd[ i , 2 ] ]
                } else {
                ## If it is positive and not negligible
                        TSPsInd[ i , ] <- c( nextPair[2] , nextPair[1] )
                        TSPtieVote[i] <- tieVotechange[ tieVote[ TSPsInd[ i , 2 ] , TSPsInd[ i , 1 ] ] + 1]
                }

                #cat("Yes\n")
                
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
        classifier$tieVote <- TSPtieVote[1:maxTSPs]

        ## Return classifier object
        return(classifier)

}


trainkTSP = function(inputMat, phenoGroup, classes = NULL, krange = 2:10,
                            FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = NULL,
                            handleTies = FALSE, disjoint = TRUE,
                            k_selection_fn = KbyTtest, k_opts = list(), score_fn = signedTSPScores, score_opts = NULL, 
                            verbose = FALSE, ...){


    # if no class labels are given, select phenoGroup factor labels as the classes
    if(is.null(classes)){
        classes = rev(levels(phenoGroup))
    }

    if(verbose) cat(sprintf("Selecting %s (case) and %s (control) as classes labels\n", classes[1], classes[2]))

    # calculate pair scores
    S = calculateScores(inputMat, phenoGroup, classes, FilterFunc, RestrictedPairs, handleTies, verbose, score_fn, score_opts, ...)

    maxk = max(krange)

    tspTable = makeTSPTable(S, maxk, disjoint)

    # select final TSPs for classifier
    sel = k_selection_fn(inputMat, phenoGroup, tspTable, classes, krange, k_opts)

    # prepare output
    classifier = list(name=sprintf("%dTSPS", length(sel)))
    classifier$TSPs = as.matrix(tspTable[sel, c("gene1", "gene2")], length(sel), 2)
    classifier$score = tspTable$score[sel]
    classifier$tieVote = tspTable$tieVote[sel]
    classifier$labels = rev(classes)

    classifier$tieVote = factor(classifier$tieVote, levels=0:2, labels=c("both", classifier$labels))

    names(classifier$score) = rownames(classifier$TSPs)
    names(classifier$tieVote) = rownames(classifier$TSPs)

    ## Return classifier object
    return(classifier)

}

makeTSPTable = function(Scores, maxk, disjoint = TRUE){

    ## If we reverse the pairs to make the score positive, the tieVote 1 and 2 become 2 and 1 respectively
    tieVotechange = c(0,2,1)


    # sorting the absolute scores works for both the signed TSP score and the basic TSP score
    absScore = abs(Scores$score)

    # select the top maxk scores
    ix = sort(absScore, decreasing = TRUE, index.return=TRUE)$ix

    # assemble top tsp data frame
    gene1 = c()
    gene2 = c()
    scores = c()
    tieVotes = c()

    for(i in 1:length(ix)){
        j = ix[i]
        pair = unlist(strsplit(names(Scores$score)[j], ","))

        # if the reverse of the pair exists already, skip
        if((pair[1] %in% gene1 && pair[2] %in% gene2) || (pair[1] %in% gene2 && pair[2] %in% gene1)) next

        if(disjoint){
            if( sum(pair %in% unique(c(gene1, gene2))) > 0 ) next
        }

        # if score is positive, reverse the order of the pair
        scores = c(scores, absScore[j])
        gene1 = c(gene1, ifelse(Scores$score[j] > 0 && Scores$signed, pair[2], pair[1]))
        gene2 = c(gene2, ifelse(Scores$score[j] > 0 && Scores$signed, pair[1], pair[2]))
        tieVotes = c(tieVotes, ifelse(Scores$score[j] > 0 && Scores$signed, tieVotechange[ Scores$tieVote[j] + 1 ], Scores$tieVote[j]))

        if(length(gene1) == maxk) break
    }

    tspTable = data.frame(gene1 = gene1, gene2 = gene2, score=scores, tieVote=tieVotes)

    return(tspTable)
}

KbyMeasurement = function(inputMat, phenoGroup, scoreTable, classes, krange, 
    k_opts=list(disjoint=TRUE, measurement="auc")){

    if(exists("disjoint", k_opts))
        disjoint = k_opts$disjoint
    else
        disjoint = TRUE

    if(exists("measurement", k_opts))
        measurement = k_opts$measurement
    else
        measurement = "auc"

    if(! (measurement %in% c("accuracy", "sensitivity", "specificity", "balanced_accuracy", "auc")))
        stop("Measurement must be one of the following: accuracy, sensitivity, specificity, balanced_accuracy, auc")

    # add pairs within krange 
    maxk = max(krange)
    fit = list(TSPs=as.matrix(scoreTable[1, 1:2], 1, 2), tieVote=scoreTable[1, 4], labels=rev(classes))

    ix = c(1)
    stat = c()

    for(i in 1:min(maxk, nrow(scoreTable))){

        if(i > 1){
            # add i-th pair to the current fit
            if(disjoint){
                if( sum(as.matrix(scoreTable[i, 1:2], 1) %in% unique(as.vector(fit$TSPs))) > 0 ) next
            }

            fit$TSPs = as.matrix(rbind(fit$TSPs, scoreTable[i, 1:2]), length(ix)+1, 2)
            fit$tieVote = c(fit$tieVote, scoreTable[i, 4])
            ix = c(ix, i)

        }

        stat = c(stat, getkTSPResult(fit=fit, Mat=inputMat, Groups=phenoGroup, classes=classes)$stats[measurement])
    }

    return(ix[ 1:which.max(stat) ])

}

#################################################
#################################################
### An internal function for finding kTSP using restricted pairs

SWAP.KTSP.Train.Restricted.Ties <- function(inputMat, phenoGroup, maxK, RestrictedPairs, verbose = FALSE) {

        ## Get the lables
        labels <- levels(phenoGroup)

        ## If we reverse the pairs to make the score positive, the tieVote 1 and 2 become 2 and 1 respectively
        tieVotechange <- c(0,2,1)

        ## Then calculate the scores for these RestrictedPairs
        scores <- calculateSignedScore(phenoGroup,
                                            inputMat1 = inputMat,
                                            inputMat2 = inputMat,
                                            RestrictedPairs = RestrictedPairs,
                                                handleTies = TRUE, verbose = verbose)

        ## Process scores
        absOrder <- order(abs(scores$score), decreasing = TRUE)

        ## Disjoint TSPs
        pairs <- SWAP.KTSP.Pairs.Disjoint(
                pairs=RestrictedPairs[ absOrder , ],
                scores$score[absOrder],
                min(maxK, length(scores$score)), scores$tieVote)

        ## Prepare output object as a list
        classifier = list();

        ## Prepare the TSPs component
        maxTSPs <- length(pairs$scores)
        classifier$TSPs <- pairs$pairs[ 1:maxTSPs , ]
        classifier$tieVotes <- pairs$tieVotes[ 1:maxTSPs]

        ## Swap gene order in the pair for pairs with negative score
        posPairs <- which(pairs$scores > 0)
        classifier$TSPs[posPairs , ] <- classifier$TSPs[posPairs , 2:1]
        classifier$tieVotes[posPairs] <- tieVotechange[classifier$tieVotes[posPairs]+1]

        ### Flipping the negative scores
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


###
### Calculate 1-TSP
###

train1TSP = function(inputMat, phenoGroup, classes = NULL, 
                            FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = NULL,
                            handleTies = FALSE, disjoint = TRUE,
                            score_fn = signedTSPScores, score_opts = NULL, 
                            verbose = FALSE, ...){

    # calculate pair scores
    S = calculateScores(inputMat, phenoGroup, classes, FilterFunc, RestrictedPairs, handleTies, verbose, score_fn, score_opts, ...)

    # if no class labels are given, select phenoGroup factor labels as the classes
    if(is.null(classes)){
        classes = rev(levels(phenoGroup))
        if(verbose) cat(sprintf("Selecting %s (class 1) and %s (class 0) as classes labels\n", classes[1], classes[2]))
    }

    j = which.max(abs(S$score))[1]

    pair = unlist(strsplit(names(S$score)[ j ], ","))
    score = S$score[ which.max(abs(S$score))[1] ]
    if(score > 0) pair = rev(pair)

    ## If we reverse the pairs to make the score positive, the tieVote 1 and 2 become 2 and 1 respectively
    tieVotechange = c(0,2,1)

    ties = ifelse(score > 0, tieVotechange[ S$tieVote[j] + 1 ], S$tieVote[j])
    score = abs(score)

    tieVote = factor(ties, levels=0:2, labels=c("both", S$labels))

    classifier = list(name="1TSP", TSPs=matrix(pair, 1), score=score, labels=S$labels, tieVote=tieVote)

    return(classifier)
}

# Old function for 1-TSP

Train.1TSP = function(...){
    scores = SWAP.CalculateSignedScore(...)

    if(is.matrix(scores$score)){
        max_i = which(abs(scores$score) == max(abs(scores$score)), arr.ind=T)[1, 1]
        max_j = which(abs(scores$score) == max(abs(scores$score)), arr.ind=T)[1, 2]
        score = scores$score[ max_i, max_j ]
        pair = c(rownames(scores$score)[max_i], colnames(scores$score)[max_j])
    }else{
        pair_str = names(scores$score)[which.max(abs(scores$score))[1]]
        pair = unlist(strsplit(pair_str, ","))
        score = scores$score[which.max(abs(scores$score))[1]]
    }

    # if the score is positive, flip the pair
    if(score > 0) pair = rev(pair)

    classifier = list(name="1TSP", TSPs=matrix(pair, 1), score=score, labels=scores$labels)

    return(classifier)

}




