pkgs <- c("rpart", "crayon", "caTools", "ada", "caret")

sapply(pkgs, require, character.only = TRUE)

setwd("C:\\Users\\jeon\\Desktop\\lab\\data")

data          <- read.csv('iris.csv', TRUE)


independent   <<- 1:4
dependent     <<- 6


# divide the dataset into input features and output feature
X             <- data[, independent]
y             <- data[, dependent]

####################################################################################################
# predict function
####################################################################################################
predict.adaboost <- function(object, X, type = c("response", "prob"), n_tree = NULL){
  
  # handle args
  type <- match.arg(type)
  
  if(is.null(n_tree)) { tree_seq <- seq_along(object$alphas) } 
  
  else                { tree_seq <- seq(1, n_tree) }
  
  # evaluate score function on sample
  f <- 0
  
  for(i in tree_seq){
    tree       <- object$trees[[i]]
    tree$terms <- object$terms
    pred       <- as.integer(as.character(stats::predict(tree, data.frame(X), type = "class")))
    f          <- f + object$alphas[i] * pred
  }
  
  # handle response type
  if(type == "response")  { sign(f) } 
  
  else if(type == "prob") { 1/(1 + exp(-2 * f)) }
}

#########################################################################################################################
### FlexBoost
#########################################################################################################################
flex.boost <- function(X, y,
                       n_rounds = 100,
                       par.k = 0.1,
                       control = rpart.control(cp = -1, maxdepth = 1)){
  
  # set parameter k
  par.k <- par.k
  
  # set 3 ways of parameter k on exp loss function
  k.list <- c(1/par.k, 1, par.k)
  
  
  # count the number of rows
  n      <- nrow(X)
  
  
  # save parameter k path, globalize path to see after function
  k.path <<- list()
  
  
  # initialize (weight, tree, alpha) on each data
  w      <- list(rep(1/n, n))
  trees  <- list()
  alphas <- list()
  
  
  # save (weight, tree, alpha) of 3 ways
  temp.w <- list()
  temp.trees <- list()
  temp.alphas <- list()
  
  
  # save train accuracy of 3 ways to compare
  temp.result <<- list()
  
  
  # build weak classifiers
  for(i in seq(n_rounds)){
    
    tree <- rpart::rpart(y ~ .,
                         data = data.frame(X), weights = w[[i]],
                         method = "class", control = control,
                         x = FALSE, y = FALSE, model = TRUE)
    
    
    pred <- as.integer(as.character(stats::predict(tree, data.frame(X), type = "class")))
    
    
    # calculate the error of each classifiers
    e    <- sum(w[[i]] * (pred != y))
    
    
    # if error >= 0.5, flip the result
    if(e >= 0.5) { e <- 1 - e }
    
    
    # count number of 3 ways
    n_count = 0
    
    for(i1 in k.list){
      
      # count number of 3 ways
      n_count = n_count + 1
      
      
      # update weak classifier weight
      alpha <- 1/(2*i1) * log((1-e)/e)
      
      
      # update and normalize weight of each data
      # save weight of 3ways
      temp.w[[n_count]]     <- w[[i]] * exp(-alpha*pred*y)
      temp.w[[n_count]]     <- temp.w[[n_count]] / sum(temp.w[[n_count]])
      
      
      # first weak classifier's weight should be 1
      if(i == 1){
        
        alphas[[i]] <- 1
        trees[[i]]  <- tree
        terms       <- tree$terms
        
      }
      
      
      # Remove formulas since they waste memory
      if(i == 1)  { terms       <- tree$terms }
      
      else        { tree$terms  <- NULL }
      
      
      alphas[[i]] <- alpha
      trees[[i]]  <- tree
      
      # save alpha, tree of 3 ways
      temp.alphas[[n_count]] <- alpha
      temp.trees[[n_count]] <-tree
      
      
      result        <- list(terms  = terms,
                            trees  = trees,
                            alphas = unlist(alphas))
      
      class(result) <- "adaboost"
      
      y_hat                   <- stats::predict(result, X)
      result$confusion_matrix <- table(y, y_hat)
      
      # save train accuracy of 3 ways 
      temp.result <- append(temp.result, (sum(diag(result$confusion_matrix)) / sum(result$confusion_matrix)), after = length(temp.result))
      
    }
    
    
    ###### compare 3 train accuracy and update (weight, alphas, tree) for next iteration
    
    # (k > 1) => max
    if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] > temp.result[[3]]){
      
      k.path      <<- append(k.path, 1/par.k, after = length(k.path))
      w[[i+1]]    <- temp.w[[1]]
      alphas[[i]] <- temp.alphas[[1]]
      trees[[i]]  <- temp.trees[[1]]
  
    }
    
    # (k = 1) => max
    else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] > temp.result[[3]]){
      
      k.path <<- append(k.path, 1, after = length(k.path))
      w[[i+1]] <- temp.w[[2]]
      alphas[[i]] <- temp.alphas[[2]]
      trees[[i]] <- temp.trees[[2]]
      
    }
    
    # (k < 1) => max
    else if (temp.result[[3]] > temp.result[[1]] & temp.result[[3]] > temp.result[[2]]){
      
      k.path <<- append(k.path, par.k, after = length(k.path))
      w[[i+1]] <- temp.w[[3]]
      alphas[[i]] <- temp.alphas[[3]]
      trees[[i]] <- temp.trees[[3]]
      
    }
    
    
    # (k > 1, k = 1) => max
    else if (temp.result[[1]] > temp.result[[3]] & temp.result[[1]] == temp.result[[2]]){
      
      k.path <<- append(k.path, 1, after = length(k.path))
      w[[i+1]] <- temp.w[[2]]
      alphas[[i]] <- temp.alphas[[2]]
      trees[[i]] <- temp.trees[[2]]
      
    }
    
    # (k > 1, k < 1) => max
    else if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] == temp.result[[3]]){
      
      k.path <<- append(k.path, par.k, after = length(k.path))
      w[[i+1]] <- temp.w[[3]]
      alphas[[i]] <- temp.alphas[[3]]
      trees[[i]] <- temp.trees[[3]]
      
    }
    
    # (k = 1, k < 1) => max
    else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] == temp.result[[3]]){
      
      k.path <<- append(k.path, 1, after = length(k.path))
      w[[i+1]] <- temp.w[[2]]
      alphas[[i]] <- temp.alphas[[2]]
      trees[[i]] <- temp.trees[[2]]
      
    }
    
    # (k > 1, k = 1,k < 1) => max
    else{
      k.path <<- append(k.path, 1, after = length(k.path))
      w[[i+1]] <- temp.w[[2]]
      alphas[[i]] <- temp.alphas[[2]]
      trees[[i]] <- temp.trees[[2]]
    }
    
    #print train accuracy on each iteration
    cat(red(max(unlist(temp.result))))
    cat("\n")
    
    # initialize each 3 ways value
    temp.w       <- list()
    temp.alphas  <- list()
    temp.trees   <- list()
    temp.result  <- list()
    
    
    
  }
  
  result        <- list(terms  = terms,
                        trees  = trees,
                        alphas = unlist(alphas))
  
  class(result) <- "adaboost"
  
  
  return(result)
}


####################################################################################################
# AdaBoost
####################################################################################################
# base model is decision stump, which splits only once
basic.adaboost <- function(X, y, n_rounds = 100,
                           control = rpart.control(cp = -1, maxdepth = 1)){
  # count the number of rows
  n      <- nrow(X)
  
  # initialize weight on each data, tree and alpha
  w      <- rep(1/n, n)
  trees  <- list()
  alphas <- list()
  
  # build weak classifiers
  for(i in seq(n_rounds)){
    
    tree <- rpart::rpart(y ~ .,
                         data = data.frame(X), weights = w,
                         method = "class", control = control,
                         x = FALSE, y = FALSE, model = TRUE)
    
    pred <- as.integer(as.character(stats::predict(tree, data.frame(X), type = "class")))
    
    # calculate the error of each classifiers
    e    <- sum(w * (pred != y))
    
    # if error >= 0.5, flip the result
    if(e >= 0.5) { e <- 1 - e }
    
    # learning rate(weight) of each classifiers
    alpha <- 1/2 * log((1-e)/e)
    
    # update and normalize weight of each data
    w     <- w * exp(-alpha*pred*y)
    w     <- w / sum(w)
    
    # If classifier's error rate is nearly 0, boosting process ends
    if(abs(e) < 1e-5){
      # if first base classifier predicts data perfectly, boosting process ends
      if(i == 1){
        # first base classifier's weight should be 1
        alphas[[i]] <- 1
        trees[[i]]  <- tree
        terms       <- tree$terms
        
        break
      }
      break
    }
    
    # Remove formulas since they waste memory
    if(i == 1)  { terms       <- tree$terms }
    
    else        { tree$terms  <- NULL }
    
    alphas[[i]] <- alpha
    trees[[i]]  <- tree
  }
  
  result        <- list(terms  = terms,
                        trees  = trees,
                        alphas = unlist(alphas))
  
  class(result) <- "adaboost"
  
  # create confusion matrix for in-sample fits
  y_hat                   <- stats::predict(result, X)
  result$confusion_matrix <- table(y, y_hat)
  
  return(result)
}


flex.boost(X, y, 30, 0.5)
