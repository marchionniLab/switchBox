###
### Utility functions
###

getPredictionStats = function(predictions, truth, classes=NULL, decision_values=NULL){
  
  if(length(predictions) != length(truth))
    stop("Predictions and true label vectors are not of the same length")

  # accuracy
  accuracy = NULL
  # balanced accuracy = (sensitivity + specificity) / 2
  balanced_accuracy = NULL
  # sensitivity = accuracy of class 1 (positive class)
  sensitivity = NULL
  # specificity = accuracy of class 0
  specificity = NULL

  # If classes are not supplied, take the factor levels of Groups argument; 
  # assume the alphanumeric order is class 0 followed by class 1, so reverse the factor levels
  # since we expect the classes argument to be (class 1, class 0)
  if(is.null(classes))
      classes = rev(levels(factor(truth)))
  
  # calculate the stats
  accuracy = sum(predictions == truth)/length(predictions)
  sensitivity = sum(predictions[ truth == classes[1] ] == truth[ truth == classes[1] ])/sum(truth == classes[1])
  specificity = sum(predictions[ truth == classes[2] ] == truth[ truth == classes[2] ])/sum(truth == classes[2])
  balanced_accuracy = (sensitivity + specificity) / 2
  
  if(is.null(decision_values)){
    
    c(accuracy=accuracy, sensitivity=sensitivity, specificity=specificity, balanced_accuracy=balanced_accuracy)
  
  }else{
    
    auc = NULL
    
    tryCatch({
      auc = auc(roc(truth, decision_values, levels=classes, direction=">"))
    }, error = function(e){
      #print(e)
    })

    c(accuracy=accuracy, sensitivity=sensitivity, specificity=specificity, balanced_accuracy=balanced_accuracy, auc=auc)
  }
  
}

getkTSPResult = function(fit, Mat, Groups, classes=NULL, predictions=FALSE, decision_values=FALSE){

  if(is.null(classes))
      classes = rev(levels(factor(Groups)))

  fit_predictions=SWAP.KTSP.Classify(Mat, fit)
  fit_decision_values=SWAP.KTSP.Statistics(Mat, fit)$statistics

  stats = getPredictionStats(predictions=fit_predictions, truth=Groups, classes=classes, decision_values=fit_decision_values)

  # classes Contains class 1 (cases) and class 0 (controls) respectively so we revere it before supplying it to 
  # the roc function; i.e. rev(classes) = (class 0, class 1) 
  # class 1 is expected to have a higher vote median than class 0, so direction is class 0 < class 1
  roc = NULL
  tryCatch({
    roc = roc(Groups, fit_decision_values, plot=FALSE, levels=rev(classes), direction="<")
  }, error = function(e){
    #print(e)
  })

  result = list(stats=stats, roc=roc)
  if(predictions)
    result$predictions = fit_predictions
  if(decision_values)
    result$decision_values = fit_decision_values

  return(result)

}

getkTSPTrainTestResults = function(trainMat, trainGroup, testMat, testGroup, classes=NULL, predictions=FALSE, decision_values=FALSE, ...){

  if(is.null(classes))
      classes = rev(levels(factor(trainGroup)))

  fit = trainkTSP(trainMat, trainGroup, ...)

  train = getkTSPResult(fit, trainMat, trainGroup, classes, predictions, decision_values)
  test = getkTSPResult(fit, testMat, testGroup, classes, predictions, decision_values)

  result = list(classifier=fit, train=train$stats, test=test$stats, trainroc=train$roc, testroc=test$roc)

  if(predictions){
    result$train_predictions = train$predictions
    result$test_predictions = test$predictions
  }
  if(decision_values){
    result$train_decision_values = train$decision_values
    result$test_decision_values = test$decision_values
  }

  return(result)

}

# perform k-fold cross-validation
# 

kTSPCV = function(Mat, Groups, classes = NULL, k = 4, folds = NULL, randomize = TRUE, ...){

  ix = 1:ncol(Mat)
  if(is.null(folds)){
    # randomly permute the order
    if(randomize)
      ix = sample(ix)

    Mat = Mat[, ix]
    Groups = Groups[ix]

    folds = get_kfold_indices(Groups, k)
  }else{
    k = length(folds)
  }

  if(is.null(classes))
      classes = rev(levels(factor(Groups)))

  cv = list()

  for(i in 1:k){
    cv[[ i ]] = getkTSPTrainTestResults(trainMat = Mat[, -folds[[i]] ], trainGroup = Groups[ -folds[[i]] ], 
      testMat = Mat[, folds[[i]] ], testGroup = Groups[ folds[[i]] ], classes = classes, 
      predictions = TRUE, decision_values = TRUE, ...)
  }

  truth = Groups[unlist(folds)]

  # combine fold results
  #predictions = as.factor( unlist(sapply(cv, function(x) as.character(x$test_predictions))) )
  predictions = unlist(sapply(cv, function(x) x$test_predictions))

  # divide no. of votes by no. of tsps in each kTSP fit
  decision_values = unlist(sapply(cv, function(x) x$test_decision_values/nrow(x$classifier$TSPs) ))

  stats = getPredictionStats(predictions = predictions, truth = truth, classes = classes, decision_values = decision_values)

  # make roc
  roc = roc(Groups, decision_values, plot=FALSE, levels=rev(classes), direction="<")

  list(cv=cv, folds=folds, predictions=predictions, decision_values=decision_values, stats=stats, roc=roc, 
    truth=truth, randomized_indices=ix)

}


# perform loo cross-validation

kTSPLOO = function(Mat, Groups, classes = NULL, ...){

  if(is.null(classes))
      classes = rev(levels(factor(Groups)))

  loo = list()

  for(i in 1:ncol(Mat)){
    loo[[ i ]] = getkTSPTrainTestResults(trainMat = Mat[, -i ], trainGroup = Groups[ -i ], 
      testMat = as.matrix(Mat[, i ], , 1), testGroup = Groups[ i ], classes = classes, 
      predictions = TRUE, decision_values = TRUE, ...)
  }

  # combine results
  #predictions = as.factor( unlist(sapply(loo, function(x) as.character(x$test_predictions))) )
  predictions = unlist(sapply(loo, function(x) x$test_predictions))

  # divide no. of votes by no. of tsps in each kTSP fit
  decision_values = unlist(sapply(loo, function(x) x$test_decision_values/nrow(x$classifier$TSPs) ))

  stats = getPredictionStats(predictions = predictions, truth = Groups, classes = classes, decision_values = decision_values)

  # make roc
  roc = roc(Groups, decision_values, plot=FALSE, levels=rev(classes), direction="<")

  list(loo=loo, predictions=predictions, decision_values=decision_values, stats=stats, roc=roc)

}











