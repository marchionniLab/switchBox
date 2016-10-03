##################################################
##################################################
### An function used to calculate the pairwise scores.
### We can restrict the pairs for which we calculate the scores.

SWAP.CalculateSignedScore <- function(inputMat, phenoGroup,
                                      FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs, handleTies = FALSE, verbose = FALSE, ...) {

        ## Check inputMat conformity
        SWAP.Check.Input(phenoGroup = phenoGroup, inputMat = inputMat,
                         RestrictedPairs = RestrictedPairs)

        ## Filter input and format the data
        formattedInput <- SWAP.Format.Input(phenoGroup =phenoGroup,
                                            inputMat = inputMat, RestrictedPairs = RestrictedPairs,
                                            FilterFunc = FilterFunc, verbose = verbose, ...)

        ## Compute the "unrestricted" classifier
        if ( missing(RestrictedPairs) ) {
                ## Obtain the TSP scores
                tryCatch(out <- calculateSignedScore(phenoGroup = phenoGroup ,
                                            inputMat1 = formattedInput$inputMat,
                                            inputMat2 = formattedInput$inputMat,
                                                handleTies = handleTies, verbose = verbose)
                ,error = function(c){ 
                    if(verbose) cat(c$message, "\nYou may consider increasing the memory or filtering features")
                })
        } else {
                ## Then calculate the scores for these RestrictedPairs
                tryCatch(out <- calculateSignedScore(phenoGroup = phenoGroup,
                                            inputMat1 = formattedInput$inputMat,
                                            inputMat2 = formattedInput$inputMat,
                                            RestrictedPairs = formattedInput$FilteredPairs,
                                                handleTies = handleTies, verbose = verbose )
                ,error = function(c){ 
                    if(verbose) cat(c$message, "\nYou may consider increasing the memory or filtering features")
                })
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

SWAP.KTSP.Train <- function(inputMat, phenoGroup, krange = 2:10, #c(3, 5, 7:10),
                            FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs, handleTies = FALSE, verbose = FALSE, ...) {

        ## Check inputMat conformity
        SWAP.Check.Input(phenoGroup = phenoGroup, inputMat = inputMat,
                         RestrictedPairs = RestrictedPairs)

        ## Filter input and format the data
        formattedInput <- SWAP.Format.Input(phenoGroup, inputMat = inputMat,
                                            RestrictedPairs = RestrictedPairs, FilterFunc = FilterFunc, ...)


        ## Prepare final output
        if( handleTies){

            ## Compute the "unrestricted" classifier
            if ( missing(RestrictedPairs) ) {
                    classifier <- SWAP.KTSP.Train.Plain.Ties(inputMat = formattedInput$inputMat,
                                                            phenoGroup = phenoGroup, maxK = max(krange), verbose = verbose )
            } else {
                    ## Develop restricted classifier
                    classifier <- SWAP.KTSP.Train.Restricted.Ties(inputMat = formattedInput$inputMat,
                                                             phenoGroup = phenoGroup, maxK = max(krange),
                                                             RestrictedPairs = formattedInput$FilteredPairs, verbose = verbose)
            }

        }else{

            ## Compute the "unrestricted" classifier
            if ( missing(RestrictedPairs) ) {
                    classifier <- SWAP.KTSP.Train.Plain(inputMat = formattedInput$inputMat,
                                                            phenoGroup = phenoGroup, maxK = max(krange), verbose = verbose )
            } else {
                    ## Develop restricted classifier
                    classifier <- SWAP.KTSP.Train.Restricted(inputMat = formattedInput$inputMat,
                                                             phenoGroup = phenoGroup, maxK = max(krange),
                                                             RestrictedPairs = formattedInput$FilteredPairs, verbose = verbose)
            }
        }

        ## Use min T-test to select k
        if(verbose) cat("Selecting K...\n")
        kmin <- bestKByTtest(classifier = classifier, inputMat = inputMat,
                                 phenoGroup = phenoGroup, krange = krange)

        ## Best K
        kstar <- min(krange[kmin], nrow(classifier$TSPs))
        if ( kstar != krange[kmin] ) {
                if(verbose) cat(paste("The required range of k is not available!\n",
                          "The minimum number of available TSP (", kstar,
                          ") will be used instead.\n", sep=""))
        } else {
                if(verbose) cat(paste(kstar, "TSP will be used to build the final classifier.\n"))
        }

        out <- list(name = sprintf("%dTSPs", kstar),
                        TSPs = classifier$TSPs[ 1:kstar , ],
                        score = classifier$score[ 1:kstar ],
                        labels = classifier$labels)

        if("tieVote" %in% names(classifier))
            out$tieVote = classifier$tieVote[ 1:kstar ]
        
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

        if(exists("tieVote",classifier)){
                ## Comparisons
                #loop over all comparisons
                comparisons <- matrix(nrow = ncol(inputMat), ncol = nrow(classifier$TSPs))
                comparisonnames <- vector(mode="character", length= nrow(classifier$TSPs))
                rownames(comparisons) <- colnames(inputMat)


                for( i in seq_len(nrow(classifier$TSPs))){
                        if(classifier$tieVote[i]=="both"){ #if the equality is not counted in favor of either classes
                                comparisons[,i] <- inputMat[classifier$TSPs[i , 1 ] , ] > inputMat[classifier$TSPs[i , 2 ] , ]
                                                        + 0.5* inputMat[classifier$TSPs[i , 1 ] , ] == inputMat[classifier$TSPs[i , 2 ] , ]
                                comparisonnames[i] = sprintf("%s>%s", classifier$TSPs[ i, 1 ], classifier$TSPs[ i, 2 ] )
                        }else if(classifier$tieVote[i] == classifier$labels[2]){#if the equality is counted in favor of class 1
                                comparisons[,i] <- inputMat[classifier$TSPs[i , 1 ] , ] > inputMat[classifier$TSPs[i , 2 ] , ]
                                comparisonnames[i] = sprintf("%s>=%s", classifier$TSPs[ i, 1 ], classifier$TSPs[ i, 2 ] )
                        }else{#if the equality is counted in favor of class 0
                                comparisons[,i] <- inputMat[classifier$TSPs[i , 1 ] , ] >= inputMat[classifier$TSPs[i , 2 ] , ]
                                comparisonnames[i] = sprintf("%s>>%s", classifier$TSPs[ i, 1 ], classifier$TSPs[ i, 2 ] )
                        }

                }

    ## Add row and column names
                colnames(comparisons) <- comparisonnames


                ## Compute feature by samples switching
                stats1 <- inputMat[classifier$TSPs[ , 1 ] , ] > inputMat[classifier$TSPs[ , 2 ] , ]
                stats2 <- inputMat[classifier$TSPs[ , 1 ] , ] < inputMat[classifier$TSPs[ , 2 ] , ]


                ## Set the combine function if missing
                if (missing(CombineFunc)) {
                        ## Compute standard K statistics
                        KTSPstat <- apply(comparisons, 1 , sum) - nrow(classifier$TSPs)/2
                        ## Add names
                        names(KTSPstat) <- colnames(inputMat)
                } else {
                        ## Apply combined function
                        KTSPstat <- apply(comparisons, 1, CombineFunc)
                }

                ## ## Make plot if needed
                ## if (show)
                ## plot the figure like the figure in the paper.
                ##
        }else{

                ## Comparisons
                comparisons <- t(inputMat[classifier$TSPs[ , 1 ] , ] > inputMat[classifier$TSPs[ , 2 ] , ])
                
                ## For 1-TSP
                if(nrow(classifier$TSPs) == 1)
                  comparisons <- matrix(comparisons, , 1)
                
                ## Add row and column names
                colnames(comparisons) <- mapply(x = classifier$TSPs[ , 1 ] ,
                                                y = classifier$TSPs [ , 2 ] ,
                                                FUN = function(x,y) sprintf("%s>%s", x, y ) )
                rownames(comparisons) <- colnames(inputMat)

                ## Compute feature by samples switching
                stats1 <- inputMat[classifier$TSPs[ , 1 ] , ] > inputMat[classifier$TSPs[ , 2 ] , ]
                stats2 <- inputMat[classifier$TSPs[ , 1 ] , ] < inputMat[classifier$TSPs[ , 2 ] , ]

                ## For 1-TSP
                if(nrow(classifier$TSPs) == 1){
                  stats1 = t(matrix(stats1, , 1))
                  stats2 = t(matrix(stats2, , 1))
                }
                
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
        }

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

### ====================================================================================================================================

###
### Replacement functions
###

SWAP.CalculateScores <- function(inputMat, phenoGroup, classes = NULL, FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = NULL, 
    handleTies = FALSE, verbose = FALSE, score_fn = signedTSPScores, score_opts = list(), ...){
    calculateScores(inputMat, phenoGroup, classes, FilterFunc, RestrictedPairs, handleTies, verbose, score_fn, score_opts, ...)
}

SWAP.Train.KTSP <- function(inputMat, phenoGroup, classes = NULL, krange = 2:10,
                            FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = NULL, 
                            handleTies = FALSE, disjoint = TRUE,
                            k_selection_fn = KbyTtest, k_opts = list(), score_fn = signedTSPScores, score_opts = NULL, 
                            verbose = FALSE, ...){

    classifier = trainkTSP(inputMat, phenoGroup, classes, krange,
                        FilterFunc, RestrictedPairs, handleTies, disjoint,
                        k_selection_fn, k_opts, score_fn, score_opts, 
                        verbose, ...)

    classifier
}

###
### New kTSP functions
###

SWAP.Train.1TSP <- function(inputMat, phenoGroup, classes = NULL, 
                            FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = NULL,
                            handleTies = FALSE, disjoint = TRUE,
                            score_fn = signedTSPScores, score_opts = NULL, 
                            verbose = FALSE, ...){

    train1TSP(inputMat, phenoGroup, classes, FilterFunc, RestrictedPairs, handleTies, disjoint, 
        score_fn, score_opts, verbose, ...)

}

SWAP.Kby.Measurement <- function(inputMat, phenoGroup, scoreTable, classes, krange, 
    k_opts=list(disjoint=TRUE, measurement="auc")){
    KbyMeasurement(inputMat, phenoGroup, scoreTable, classes, krange, k_opts)
}

SWAP.Kby.Ttest <-  function(inputMat, phenoGroup, scoreTable, classes, krange, k_opts=list()){
    KbyTtest(inputMat, phenoGroup, scoreTable, classes, krange, k_opts)
}

SWAP.MakeTSPTable <- function(Scores, maxk, disjoint = TRUE){
    makeTSPTable(Scores, maxk, disjoint)
}

###
### Score functions
###

SWAP.Calculate.SignedTSPScores <- function(phenoGroup, inputMat1, inputMat2 = NULL, classes = NULL, RestrictedPairs = NULL, 
    handleTies = FALSE, verbose = FALSE, score_opts=list()){
    signedTSPScores(phenoGroup, inputMat1, inputMat2, classes, RestrictedPairs, handleTies, verbose, score_opts)
}

SWAP.Calculate.BasicTSPScores <- function(phenoGroup, inputMat1, inputMat2 = NULL, classes = NULL, RestrictedPairs = NULL, 
    handleTies = FALSE, verbose = FALSE, score_opts=list()){
    basicTSPScores(phenoGroup, inputMat1, inputMat2, classes, RestrictedPairs, handleTies, verbose, score_opts)
}

###
### Plotting functions
###

SWAP.PlotKTSP.Genes <- function(inputMat, Groups, classes, genes, colors=c(), legends=c(), ...){
    plotGenes(inputMat, Groups, classes, genes, colors, legends, ...)
}

SWAP.PlotKTSP.GenePairScatter <- function(inputMat, Groups, classes, genes, colors=c(), legends=c(), ...){
    plotGenePairScatter(inputMat, Groups, classes, genes, colors, legends, ...)
}

SWAP.PlotKTSP.TrainTestROC <- function(result, colors=c(), legends=c(), ...){
    plotkTSPTrainTestROC(result, colors, legends, ...)
}

SWAP.PlotKTSP.Votes <- function(classifier, inputMat, Groups=NULL, CombineFunc, ...){
    plotkTSPVotes(classifier, inputMat, Groups, CombineFunc, ...)
}

SWAP.PlotKTSP.GenePairBoxplot <- function(genes, inputMat, Groups=NULL, classes=NULL, points=FALSE, point_coloring="byGene", colors=c(), point_colors=c(), ...){
    plotGenePairBoxplot(genes, inputMat, Groups, classes, points, point_coloring, colors, point_colors, ...)
}

SWAP.PlotKTSP.GenePairClassesBoxplot <- function(genes, inputMat, Groups, classes=NULL, points=FALSE, ordering="byGene", colors=c(), point_colors=c(), 
  point_directions=FALSE, ...){
    plotGenePairClassesBoxplot(genes, inputMat, Groups, classes, points, ordering, colors, point_colors, point_directions, ...)
}

###
### Utility functions
###

SWAP.GetKTSP.PredictionStats <- function(predictions, truth, classes=NULL, decision_values=NULL){
    getPredictionStats(predictions, truth, classes, decision_values)
}

SWAP.GetKTSP.Result <- function(classifier, inputMat, Groups, classes=NULL, predictions=FALSE, decision_values=FALSE){
    getkTSPResult(classifier, inputMat, Groups, classes, predictions, decision_values)
}

SWAP.GetKTSP.TrainTestResults <- function(trainMat, trainGroup, testMat, testGroup, classes=NULL, predictions=FALSE, decision_values=FALSE, ...){
    getkTSPTrainTestResults(trainMat, trainGroup, testMat, testGroup, classes, predictions, decision_values, ...)
}

SWAP.KTSP.CV <- function(inputMat, Groups, classes = NULL, k = 4, folds = NULL, randomize = TRUE, ...){
    kTSPCV(inputMat, Groups, classes, k, folds, randomize, ...)
}

SWAP.KTSP.LOO <- function(inputMat, Groups, classes = NULL, ...){
    kTSPLOO(inputMat, Groups, classes, ...)
}

###
### Miscellaneous functions
###

SWAP.ScoreMatrixToVector <- function(M){
    score_matrix_to_vector(M)
}

SWAP.ScoreVectorToMatrix <- function(V){
    score_vector_to_matrix(V)
}

SWAP.MakeTrainTestData <- function(inputMat, Groups, classes = NULL, p = .5){
    make_train_test_data(inputMat, Groups, classes, p)
}

SWAP.GetKFoldIndices <- function(Groups, k = 4){
    get_kfold_indices(Groups, k)
}




