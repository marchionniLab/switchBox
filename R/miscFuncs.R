###
### Miscellaneous functions
###

# convert a matrix of scores to a vector with names for each element
score_matrix_to_vector = function(M){

  V = as.vector(M)
  names(V) = as.vector(t(sapply(rownames(M), function(x) sapply(colnames(M), function(y) sprintf("%s,%s", x, y) ))))

  V
}

# convert a vector of scores to a matrix
score_vector_to_matrix = function(V){

  rows = sort(unique(unname(sapply( names(V), function(x) unlist(strsplit(x, ","))[1] ))))
  cols = sort(unique(unname(sapply( names(V), function(x) unlist(strsplit(x, ","))[2] ))))

  M = matrix(NA, length(rows), length(cols), dimnames = list(rows, cols))
  for(i in 1:length(V)){
    pair = unlist(strsplit(names(V)[i], ","))
    M[pair[1], pair[2]] = V[i]
  }

  M
}

make_train_test_data = function(Mat, Groups, classes = NULL, p = .5){

  if(is.null(classes))
    classes = rev(levels(Groups))

  sel0 = sample( 1:sum(Groups == classes[2]), floor( sum(Groups == classes[2]) * p ) )
  sel1 = sample( 1:sum(Groups == classes[1]), floor( sum(Groups == classes[1]) * p ) )

  train_ids = c( which( Groups == classes[2] )[ sel0 ], which( Groups == classes[1] )[ sel1 ] )
  test_ids = c( which( Groups == classes[2] )[ -sel0 ], which( Groups == classes[1] )[ -sel1 ] )

  trainMat = Mat[ , train_ids ]
  trainGroup = Groups[ train_ids ]

  testMat = Mat[ , test_ids ]
  testGroup = Groups[ test_ids ]

  list(trainMat=trainMat, testMat=testMat, trainGroup=trainGroup, testGroup=testGroup, train_ids=train_ids, test_ids=test_ids, classes=classes)

}

calculate_tsp_instance = function(xi, xj, Groups, levels){
    (sum(xi[ Groups == levels[1] ] < xj[ Groups == levels[1] ])/sum(Groups == levels[1])) + 
        (sum(xi[ Groups == levels[2] ] > xj[ Groups == levels[2] ])/sum(Groups == levels[2]))
}

get_kfold_indices = function(y, k=4){

  # list contains samples for each fold
  folds = list()

  # make a list of samples by class label
  for(label_i in lapply(unique(y), function(y_i) sample( which(y == y_i) ) )){
    # for each class label, separate the samples into k groups
    
    t = (1:length(label_i) %% k ) + 1

    for(j in 1:k){
      # append samples into each fold
      if(length(folds) < j)
        folds[[j]] = label_i[ which(t == j)]
      else
        folds[[j]] = c(folds[[j]], label_i[ which(t == j) ])
    }
  }

  folds
}

