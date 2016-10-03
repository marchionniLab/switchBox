### ====== Jitter amount ======
JITTER = 5

### ====== Default colours =======
COLOUR1 = "cyan4" # teal like colour
COLOUR2 = "coral3" # salmon colour
COLOUR3 = "orchid4" # purple colour
COLOUR4 = "darkolivegreen4" # olive green colour

###
### Plotting functions
###

plotGenes = function(Mat, Groups, classes, genes, colors=c(), legends=c(), ...){
  
  # Check if genes are supplied
  if(length(genes) < 1)
    stop("Gene names missing")

  if(! all(genes %in% rownames(Mat)))
    stop("Given gene names missing in data matrix rownames")

  # If colours are not supplied, use 1, 2, 3, ... as the colours.
  if(length(colors) != length(genes))
    colors = 1:length(genes)
  
  # If gene lengeds are not supplied use ""
  if(length(legends) != length(genes))
    legends = rep("", length(genes))

  # If classes are not supplied, take the factor levels of Groups argument; 
  # assume the alphanumeric order is class 0 followed by class 1, so reverse the factor levels
  # since we expect the classes argument to be (class 1, class 0)
  if(missing(classes))
    classes = rev(levels(Groups))
  
  # Re-arrange data in order of classes
  Mat = Mat[, c(which(Groups == classes[2]), which(Groups == classes[1]))]
  
  # Class 1 is drawn in circles and Class 0 is drawn in triangles
  pch = rep(c(19, 17), c(sum(Groups == classes[2]), sum(Groups == classes[1])))

  # Plot 
  for(j in 1:length(genes)){
    if(j == 1){
      y = Mat[genes[j], ]
      plot(y, type="o", col=colors[j], pch=pch, #xlab=xlab, ylab=ylab,
           ylim=c(min(Mat[genes, ]), max(Mat[genes, ])),
           ...)
    }else
      lines(Mat[genes[j], ], type="o", col=colors[j], pch=pch)
  }

  # Add grid
  grid(20, 20, col="darkgray")

  # Add vertical line separating the two classes
  abline(v=sum(Groups == classes[2])+.5, lwd=2)

  # Add legends
  legend("topleft", 
         sapply(1:length(genes), function(j)  ifelse(legends[j]=="", genes[j], paste(genes[j], "(", legends[j], ")", sep="") ) ),
         lwd=4, col=colors, bty="n")
  legend("topright", classes, pch=c(17, 19), col=1, bty="n")
  
}

plotGenePairScatter = function(Mat, Groups, classes, genes, colors=c(), legends=c(), ...){
  
  # Check if genes are supplied
  if(length(genes) < 2)
    stop("Requires two gene names")

  if(! all(genes %in% rownames(Mat)))
    stop("Given gene names missing in data matrix rownames")

  # If more than 2 genes are supplied, we take the first 2
  if(length(genes) > 2)
    genes = genes[1:2]

  # If classes are not supplied, take the factor levels of Groups argument; 
  # assume the alphanumeric order is class 0 followed by class 1, so reverse the factor levels
  # since we expect the classes argument to be (class 1, class 0)
  if(missing(classes))
    classes = rev(levels(Groups))

  # Use default colours if necessary
  if(length(colors) != length(classes))
    colors = c(COLOUR1, COLOUR2)
  
  # If class lengeds are not supplied use ""
  if(length(legends) != length(classes))
    legends = rep("", length(classes))
  
  pch = as.numeric(as.character(factor(Groups, labels=c(17, 19), levels=classes)))
  col = as.character(factor(Groups, labels=colors, levels=classes))

  plot(Mat[genes[1], ], Mat[genes[2], ], pch=pch, col=col, xlab=genes[1], ylab=genes[2], ...)
  abline(a = 0, b = 1, col="black", lty=3, lwd=2)

  grid(20, 20, col="darkgray")

  legend("topleft", 
    sapply(1:length(classes), function(j)  ifelse(legends[j]=="", classes[j], paste(classes[j], "(", legends[j], ")", sep="") ) ),
    pch=c(17, 19), col=colors)
  
}

plotkTSPTrainTestROC = function(result, colors=c(), legends=c(), ...){

  # result is list with a trainroc element and testroc element, and each element is an roc object produced by pROC; 
  # alternatively, the object returned by getkTSPTrainTestResults()

  if( ! all( c("trainroc", "testroc") %in% names(result) ) )
    stop("'result' argument needs to have a trainroc element and a testroc element 
      corresponding to a training ROC curve and a testing ROC curve respectively")

  # Use default colours if necessary
  if(length(colors) != 2)
    colors = c(COLOUR1, COLOUR2)

  # Add empty legends
  if(length(legends) != 2)
    legends = rep("", 2)

  # Plot the roc curves
  plot(result$trainroc, col=colors[1], lwd=2, ...)
  lines(result$testroc, col=colors[2], lwd=2)
  
  legend_str = c(
    paste("Train AUC=", round(auc(result$trainroc), 3), sep=""),
    paste("Test AUC=", round(auc(result$testroc), 3), sep="")
  )

  legend_str = sapply( 1:2, function(i) ifelse(legends[i] == "", legend_str[i], sprintf("%s (%s)", legend_str[i], legends[i]) ) )

  legend("bottomright", legend_str, col=colors, lwd=2)

  grid(20, 20, col="darkgray")

}

plotkTSPVotes = function(fit, Mat, Groups=NULL, CombineFunc, ...){

  # Calculate votes
  votes = t(SWAP.KTSP.Statistics(Mat, fit, CombineFunc)$comparisons * 1)
  
  # If group labels are available, they'll be used instead of sample names
  if(! is.null(Groups))
    colnames(votes) = Groups

  # plot heatmap
  plot.new()

  heatmap.2(votes, Rowv=FALSE, Colv=TRUE, trace="none", dendrogram="column", ...)

}

plotGenePairBoxplot = function(genes, Mat, Groups=NULL, classes=NULL, points=FALSE, point_coloring="byGene", colors=c(), point_colors=c(), ...){

  # Check if genes are supplied
  if(length(genes) != 2)
    stop("Requires two gene names")

  if(! all(genes %in% rownames(Mat)))
    stop("Given gene names missing in data matrix rownames")

  # Use default colours if necessary
  if(length(colors) != 2)
    colors = c(COLOUR1, COLOUR2)
  if(length(point_colors) != 2)
    point_colors = c(COLOUR3, COLOUR4)
     
  l = list(Mat[genes[1], ], Mat[genes[2], ])
  names(l) = genes

  boxplot(l, border=colors, lwd=2, ...)
  if(points){

    if(point_coloring == "byGene"){
      stripchart(l, vertical = TRUE, method="jitter", jitter=JITTER/20, pch=19, add=TRUE, col=point_colors)
    }else if(point_coloring == "byDirection"){
      point_dir_colors = ifelse(l[[1]] < l[[2]], point_colors[1], point_colors[2])
      points(x=jitter(rep(1, length(l[[1]])), factor=JITTER), y=l[[1]], pch=19, col=point_dir_colors)
      points(x=jitter(rep(2, length(l[[1]])), factor=JITTER), y=l[[2]], pch=19, col=point_dir_colors)
      legend("topright", c(sprintf("%s < %s", genes[1], genes[2]), sprintf("%s >= %s", genes[1], genes[2])), pch=19, col=point_colors)
    }else if(point_coloring == "byClass"){
      if(is.null(Groups)){
        warning("Groups not available for plotting individual sample points")
      }else{
        if(is.null(classes))
          classes = rev(levels(Groups))

        point_group_colors = as.character(factor(Groups, labels=point_colors, levels=classes))
        points(x=jitter(rep(1, length(l[[1]])), factor=JITTER), y=l[[1]], pch=19, col=point_group_colors)
        points(x=jitter(rep(2, length(l[[1]])), factor=JITTER), y=l[[2]], pch=19, col=point_group_colors)
        legend("topright", classes, pch=19, col=point_colors)
      }
    }else{
      warning("Invalid point coloring argument")
    }
  }

  grid(20, 20, col="darkgray")

}

plotGenePairClassesBoxplot = function(genes, Mat, Groups, classes=NULL, points=FALSE, ordering="byGene", colors=c(), point_colors=c(), 
  point_directions=FALSE, ...){

  # Check if genes are supplied
  if(length(genes) != 2)
    stop("Requires two gene names")

  if(! all(genes %in% rownames(Mat)))
    stop("Given gene names missing in data matrix rownames")

  if(is.null(classes))
    classes = rev(levels(Groups))

  # Use default colours if necessary
  if(length(colors) != 2)
    colors = c(COLOUR1, COLOUR2)
  if(length(point_colors) != 2)
    point_colors = colors #c("black", "gray")
  
  if(ordering == "byGene"){

    # We draw each gene as two adjacent boxplots - one for each class
    l = list(Mat[genes[1], Groups == classes[2]], Mat[genes[1], Groups == classes[1]], 
      Mat[genes[2], Groups == classes[2]], Mat[genes[2], Groups == classes[1]])

    names(l) = c(sprintf("%s, %s", genes[1], classes[2]), sprintf("%s, %s", genes[1], classes[1]), 
      sprintf("%s, %s", genes[2], classes[2]), sprintf("%s, %s", genes[2], classes[1]))

    # Note: plotting the x-axis separately using the axis() functions means the user cannot alter the behavior of the
    # x-axis labels since the optional arguments (i.e.: '...') are all passed onto the boxplot() function

    boxplot(l, border=c(colors, colors), lwd=2, xaxt="n", ...)
    axis(1, at=c(1.5, 3.5), labels=genes)
    legend("topright", classes, col=colors, lwd=2)
    abline(v=2.5, lwd=2, col="lightgray")
    grid(20, 20, col="darkgray")

    if(points){

      if(point_directions){
            point_group_colors = lapply(rep(rev(classes), 2), function(c) ifelse(
              Mat[genes[1], Groups == c] < Mat[genes[2], Groups == c], point_colors[1], point_colors[2]))
            legend("topleft", c(sprintf("%s < %s", genes[1], genes[2]), sprintf("%s >= %s", genes[1], genes[2])), pch=19, col=point_colors)
      }else{
        point_group_colors = mapply(function(x, y) rep(x, length(y)), rep(point_colors, 2), l)
      }

      for(j in 1:length(l))
        points(x=jitter(rep(j, length(l[[j]])), factor=JITTER), y=l[[j]], pch=19, col=point_group_colors[[j]])

    }

  }else if(ordering == "byClass"){

    # We draw each class as two adjacent boxplots - one for each gene
    l = list(Mat[genes[1], Groups == classes[2]], Mat[genes[2], Groups == classes[2]], 
      Mat[genes[1], Groups == classes[1]], Mat[genes[2], Groups == classes[1]])

    names(l) = c(sprintf("%s, %s", genes[1], classes[2]), sprintf("%s, %s", genes[2], classes[2]), 
      sprintf("%s, %s", genes[1], classes[1]), sprintf("%s, %s", genes[2], classes[1]))

    boxplot(l, border=c(colors, colors), lwd=2, xaxt="n", ...)
    axis(1, at=c(1.5, 3.5), labels=rev(classes))
    legend("topright", genes, col=colors, lwd=2)
    abline(v=2.5, lwd=2, col="lightgray")
    grid(20, 20, col="darkgray")

    if(points){

      if(point_directions){
            point_group_colors = lapply(rep(rev(classes), each=2), function(c) ifelse(
              Mat[genes[1], Groups == c] < Mat[genes[2], Groups == c], point_colors[1], point_colors[2]))
            legend("topleft", c(sprintf("%s < %s", genes[1], genes[2]), sprintf("%s >= %s", genes[1], genes[2])), pch=19, col=point_colors)
      }else{
        point_group_colors = mapply(function(x, y) rep(x, length(y)), rep(point_colors, 2), l)
      }

      for(j in 1:length(l))
        points(x=jitter(rep(j, length(l[[j]])), factor=JITTER), y=l[[j]], pch=19, col=point_group_colors[[j]])
    }

  }else{
    stop("Invalid ordering argument")
  }

}



