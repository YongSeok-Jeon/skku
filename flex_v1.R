pkgs <- c("rpart", "ada", "caTools", "crayon", "ggplot2", "gplots", "PairedData", "car", "tidyr", "reshape")
sapply(pkgs, require, character.only = T)


load("Dataset.RData")
load("Result.RData")

data <- Haberman

independent   <<- 1:3
dependent     <<- 4

# set seed
s <<- 9 


# divide the dataset into input features and output feature
X             <- data[, independent]
y             <- data[, dependent]


####################################################################################################
# AdaBoost
####################################################################################################
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

####################################################################################################
# FlexBoost
####################################################################################################
FlexBoost <- function(X, y, n_rounds, par.k, type, control = rpart.control(cp = -1, maxdepth = 1)){
  
  
  # set 3 ways of parameter k on exp loss function
  k.list      <- c(1/par.k, 1, par.k)
  
  
  # count the number of rows
  n           <- nrow(X)
  
  
  # save parameter k path, globalize path to see after function
  k.path      <- list()
  
  
  # initialize (weight, tree, alpha) on each data
  w           <- list(rep(1/n, n))
  trees       <- list()
  alphas      <- list()
  
  
  # save (weight, tree, alpha) of 3 ways
  temp.w      <- list()
  temp.trees  <- list()
  temp.alphas <- list()
  
  
  # save train accuracy of 3 ways to compare
  temp.result <- list()
  
  # save train accuracy
  acc.result  <- list()
  
  
  # build weak classifiers
  for(i in seq(n_rounds)){
    
    tree <- rpart::rpart(y ~ .,
                         data = data.frame(X), weights = w[[i]],
                         method = "class", control = control,
                         x = FALSE, y = FALSE, model = TRUE)
    
    
    pred      <- as.integer(as.character(stats::predict(tree, data.frame(X), type = "class")))
    
    
    # calculate the error of each classifiers
    e         <- sum(w[[i]] * (pred != y))
    
    
    # if error >= 0.5, flip the result
    if(e >= 0.5){e <- 1 - e}
    
    if(e < 1e-5){
      # if first base classifier predicts data perfectly, boosting process ends
      if(i == 1){
        # first base classifier's weight should be 1
        alphas[[i]]     <- 1
        trees[[i]]      <- tree
        terms           <- tree$terms
        acc.result[[i]] <- 1
        break
      }
      break
    }
    
    # count number of 3 ways
    n_count  <- 0
    
    for(i1 in k.list){
      
      # count number of 3 ways
      n_count <- n_count + 1
      
      
      # update weak classifier weight
      alpha   <- 1/(2*i1) * log((1-e)/e)
      
      
      # update and normalize weight of each data
      # save weight of 3ways
      temp.w[[n_count]]     <- w[[i]] * exp(-alpha*pred*y*i1)
      temp.w[[n_count]]     <- temp.w[[n_count]] / sum(temp.w[[n_count]])
      
      
      # Remove formulas since they waste memory
      if(i == 1)  { terms       <- tree$terms }
      
      else        { tree$terms  <- NULL }
      
      
      alphas[[i]] <- alpha
      trees[[i]]  <- tree
      
      # save alpha, tree of 3 ways
      temp.alphas[[n_count]] <- alpha
      temp.trees[[n_count]]  <-tree
      
      
      result                 <- list(terms  = terms,
                                     trees  = trees,
                                     alphas = unlist(alphas))
      
      class(result) <- "FlexBoost"
      
      y_hat                   <- stats::predict(result, X)
      
      
      # save train accuracy
      temp.result[[n_count]]  <- sum(y == y_hat) / length(y)
    }
    
    # clear waste memory
    w[[i]] <- NULL
    
    # if first classifier splits perfect, stop iteration
    if(max(unlist(temp.result)) == 1){
      acc.result[[i]] <- 1
      break
    }
    
    acc.result[[i]] <- max(unlist(temp.result)) 
    
    if (type == 1){
      ###### compare 3 train accuracy and update (weight, alphas, tree) for next iteration
      
      # (k > 1) => max
      if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] > temp.result[[3]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k = 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] > temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k < 1) => max
      else if (temp.result[[3]] > temp.result[[1]] & temp.result[[3]] > temp.result[[2]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k > 1, k = 1) => max
      else if (temp.result[[1]] > temp.result[[3]] & temp.result[[1]] == temp.result[[2]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k > 1, k < 1) => max
      else if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] == temp.result[[3]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k = 1, k < 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] == temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k > 1, k = 1, k < 1) => max
      else{
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
      }
      
      
      # initialize each 3 ways value
      temp.w       <- list()
      temp.alphas  <- list()
      temp.trees   <- list()
      temp.result  <- list()
      
    }
    
    else if(type == 2){
      ###### compare 3 train accuracy and update (weight, alphas, tree) for next iteration
      
      # (k > 1) => max
      if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] > temp.result[[3]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k = 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] > temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k < 1) => max
      else if (temp.result[[3]] > temp.result[[1]] & temp.result[[3]] > temp.result[[2]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k > 1, k = 1) => max
      else if (temp.result[[1]] > temp.result[[3]] & temp.result[[1]] == temp.result[[2]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k > 1, k < 1) => max
      else if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] == temp.result[[3]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k = 1, k < 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] == temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k > 1, k = 1, k < 1) => max
      else{
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
      }
      
      
      # initialize each 3 ways value
      temp.w       <- list()
      temp.alphas  <- list()
      temp.trees   <- list()
      temp.result  <- list()
      
    }
    
    else{
      ###### compare 3 train accuracy and update (weight, alphas, tree) for next iteration
      
      # (k > 1) => max
      if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] > temp.result[[3]]){
        
        k.path[[i]] <- 1/par.k
        w[[i+1]]    <- temp.w[[1]]
        alphas[[i]] <- temp.alphas[[1]]
        trees[[i]]  <- temp.trees[[1]]
        
      }
      
      # (k = 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] > temp.result[[3]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k < 1) => max
      else if (temp.result[[3]] > temp.result[[1]] & temp.result[[3]] > temp.result[[2]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k > 1, k = 1) => max
      else if (temp.result[[1]] > temp.result[[3]] & temp.result[[1]] == temp.result[[2]]){
        
        k.path[[i]] <- 1
        w[[i+1]]    <- temp.w[[2]]
        alphas[[i]] <- temp.alphas[[2]]
        trees[[i]]  <- temp.trees[[2]]
        
      }
      
      # (k > 1, k < 1) => max
      else if (temp.result[[1]] > temp.result[[2]] & temp.result[[1]] == temp.result[[3]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k = 1, k < 1) => max
      else if (temp.result[[2]] > temp.result[[1]] & temp.result[[2]] == temp.result[[3]]){
        
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
        
      }
      
      # (k > 1, k = 1, k < 1) => max
      else{
        k.path[[i]] <- par.k
        w[[i+1]]    <- temp.w[[3]]
        alphas[[i]] <- temp.alphas[[3]]
        trees[[i]]  <- temp.trees[[3]]
      }
      
      
      # initialize each 3 ways value
      temp.w       <- list()
      temp.alphas  <- list()
      temp.trees   <- list()
      temp.result  <- list()
      
    }
    
    
    
    
  }
  
  result        <- list(terms  = terms,
                        trees  = trees,
                        alphas = unlist(alphas),
                        acc = unlist(acc.result))
  
  class(result) <- "FlexBoost"
  
  return(result)
}


####################################################################################################
# Predict Function AdaBoost
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


####################################################################################################
# Predict Function FlexBoost
####################################################################################################
predict.FlexBoost <- function(object, X, type = c("response", "prob"), n_tree = NULL){
  
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


####################################################################################################
# K-Fold AdaBoost
####################################################################################################
kfold.ada <- function(iteration){
  # number of cross validation
  k       <- 5
  
  # fix the seed so that result could be reproducible
  set.seed(1) 
  
  # each data has its own id(1 to 5) to process k fold cross validation 
  data$id <- sample(1:k, nrow(data), replace = TRUE)
  list    <- 1:k
  
  summ <<- 0
  #function for k fold
  for(i in 1:k){
    prediction_ada <- testset_copy_ada <- data.frame()
    
    # divide the whole dataset into train and testset
    trainset     <- subset(data, id %in% list[-i])
    testset      <- subset(data, id %in% c(i))
    
    #run a adaboost model
    model_ada                <- basic.adaboost(trainset[, independent], trainset[, dependent], n_rounds = iteration)
    
    # predict
    temp_ada                 <- as.data.frame(predict(model_ada, testset))
    
    # append this iteration's prediction to the end of the prediction data frame
    prediction_ada           <- rbind(prediction_ada, temp_ada)
    
    # append this iteration's test set to the testset copy data frame
    testset_copy_ada         <- rbind(testset_copy_ada, as.data.frame(testset[, dependent]))
    
    # result
    result_ada               <- cbind(prediction_ada, testset_copy_ada[, 1])
    
    # confustion matrix and accuracy
    names(result_ada)        <- c("Actual", "Predicted")
    confusion_matrix_ada     <- table(result_ada$Actual, result_ada$Predicted)
    accuracy_ada             <- sum(diag(confusion_matrix_ada)) / sum(confusion_matrix_ada)
    summ <<- summ + accuracy_ada
  }  
  
  # store the accuracy of each model
  temp.ada <<- summ / 5
  
  return(temp.ada)
}



####################################################################################################
# K-Fold LogitBoost
####################################################################################################
kfold.logit <- function(iteration){
  # number of cross validation
  k       <- 5
  
  # fix the seed so that result could be reproducible
  set.seed(s)
  
  # each data has its own id(1 to 5) to process k fold cross validation 
  data$id <- sample(1:k, nrow(data), replace = TRUE)
  list    <- 1:k
  
  summ <<- 0
  #function for k fold
  for(i in 1:k){
    prediction_logit <- testset_copy_logit <- data.frame()
    
    # divide the whole dataset into train and testset
    trainset     <- subset(data, id %in% list[-i])
    testset      <- subset(data, id %in% c(i))
    
    #run a adaboost model
    model_logit              <- LogitBoost(trainset[, independent], trainset[, dependent], nIter = iteration)
    
    # predict
    temp_logit               <- as.data.frame(predict(model_logit, testset, type = "raw"))
    
    
    
    # if prob of expected classes are same as 0.5, randomly define them as prob 0.5 
    for (i1  in 1:nrow(testset)){
      if (temp_logit[i1, 1] > 0.5){
        prediction_logit         <- rbind(prediction_logit, -1)
      }
      else if (temp_logit[i1, 1] == 0.5){
        prediction_logit         <- rbind(prediction_logit, (sample(c(-1,1), size = 1,prob = c(0.5, 0.5)))  )
        
      }
      else{
        prediction_logit         <- rbind(prediction_logit, 1)
      }
      
    }
    
    
    
    # append this iteration's test set to the testset copy data frame
    testset_copy_logit       <- rbind(testset_copy_logit, as.data.frame(testset[, dependent]))
    
    # result
    result_logit             <- cbind(prediction_logit, testset_copy_logit[, 1])
    
    # confustion matrix and accuracy
    names(result_logit)      <- c("Actual", "Predicted")
    confusion_matrix_logit   <- table(result_logit$Actual, result_logit$Predicted)
    accuracy_logit           <- sum(diag(confusion_matrix_logit)) / sum(confusion_matrix_logit)
    
    
    summ <<- summ + accuracy_logit
  }  
  
  # store the accuracy of each model
  temp.logit  <<- summ / 5
  return (temp.logit)
}



####################################################################################################
# K-Fold GentleBoost
####################################################################################################
kfold.gentle <- function(iteration){
  # number of cross validation
  k       <- 5
  
  # fix the seed so that result could be reproducible
  set.seed(1) 
  
  # each data has its own id(1 to 5) to process k fold cross validation 
  data$id <- sample(1:k, nrow(data), replace = TRUE)
  list    <- 1:k
  
  summ <<- 0
  #function for k fold
  for(i in 1:k){
    prediction_gentle <- testset_copy_gentle <- data.frame()
    
    # divide the whole dataset into train and testset
    trainset     <- subset(data, id %in% list[-i])
    testset      <- subset(data, id %in% c(i))
    
    #run a adaboost model
    model_gentle             <- ada(trainset[, independent], trainset[, dependent], loss = "exponential", type = "gentle", iter = iteration)
    
    # predict
    temp_gentle              <- as.data.frame(predict(model_gentle, testset))
    
    # append this iteration's prediction to the end of the prediction data frame
    prediction_gentle        <- rbind(prediction_gentle, temp_gentle)
    
    # append this iteration's test set to the testset copy data frame
    testset_copy_gentle      <- rbind(testset_copy_gentle, as.data.frame(testset[, dependent]))
    
    # result
    result_gentle            <- cbind(prediction_gentle, testset_copy_gentle[, 1])
    
    # confustion matrix and accuracy
    names(result_gentle)     <- c("Actual", "Predicted")
    confusion_matrix_gentle  <- table(result_gentle$Actual, result_gentle$Predicted)
    accuracy_gentle          <- sum(diag(confusion_matrix_gentle)) / sum(confusion_matrix_gentle)
    summ <<- summ + accuracy_gentle
  }  
  
  temp.gentle <<- summ / 5
  
  return(temp.gentle)
}



####################################################################################################
# K-Fold FlexBoost
####################################################################################################
kfold.flex <- function(iteration, par.k, type){
  # number of cross validation
  k       <- 5
  
  # fix the seed so that result could be reproducible
  set.seed(s) 
  
  # each data has its own id(1 to 5) to process k fold cross validation 
  data$id <- sample(1:k, nrow(data), replace = TRUE)
  list    <- 1:k
  
  summ <<- 0
  #function for k fold
  for(i in 1:k){
    prediction_flex <- testset_copy_flex <- data.frame()
    
    # divide the whole dataset into train and testset
    trainset     <- subset(data, id %in% list[-i])
    testset      <- subset(data, id %in% c(i))
    
    #run a adaboost model
    model_flex             <- FlexBoost(trainset[, independent], trainset[, dependent], iteration, par.k, type)
    
    # predict
    temp_flex              <- as.data.frame(predict(model_flex, testset))
    
    # append this iteration's prediction to the end of the prediction data frame
    prediction_flex        <- rbind(prediction_flex, temp_flex)
    
    # append this iteration's test set to the testset copy data frame
    testset_copy_flex      <- rbind(testset_copy_flex, as.data.frame(testset[, dependent]))
    
    # result
    result_flex            <- cbind(prediction_flex, testset_copy_flex[, 1])
    
    # confustion matrix and accuracy
    names(result_flex)     <- c("Actual", "Predicted")
    confusion_matrix_flex  <- table(result_flex$Actual, result_flex$Predicted)
    accuracy_flex          <- sum(diag(confusion_matrix_flex)) / sum(confusion_matrix_flex)
    summ <<- summ + accuracy_flex
  }  
  
  temp.flex <<- summ / 5
  
  return(temp.flex)
}



####################################################################################################
# Plot Function
####################################################################################################
plot.result <- function(col){
  ff  <- result_flex[,col]
  ll  <- result_logit[,col]
  bb  <- result_ada[,col]
  gg  <- result_gentle[,col]
  plot(1:100, ff, ylim = c(min(ff,ll,bb,gg), max(ff,ll,bb,gg)), xlab = "Iteration", ylab = "Accuracy", col = "red")
  lines(1:100, ff, col = "red")
  par(new = TRUE)
  plot(1:100, ll, xlab = "", ylab = "", col = "blue", ylim = c(min(ff,ll,bb,gg), max(ff,ll,bb,gg)))
  lines(1:100, ll, col = "blue")
  par(new = TRUE)
  plot(1:100, bb, xlab = "", ylab = "", col = "black", ylim = c(min(ff,ll,bb,gg), max(ff,ll,bb,gg)))
  lines(1:100, bb, col = "black")
  par(new = TRUE)
  plot(1:100, gg, xlab = "", ylab = "", col = "green3", ylim = c(min(ff,ll,bb,gg), max(ff,ll,bb,gg)))
  lines(1:100, gg, col = "green3")
}


####################################################################################################
# Friedman Function
####################################################################################################
friedman.post.hoc <- function(formu, data, to.print.friedman = T, to.post.hoc.if.signif = T,  to.plot.parallel = T, to.plot.boxplot = T, signif.P = .05, color.blocks.in.cor.plot = T, jitter.Y.in.cor.plot =F)
{
  # formu is a formula of the shape:    Y ~ X | block
  # data is a long data.frame with three columns:    [[ Y (numeric), X (factor), block (factor) ]]
  
  # Note: This function doesn't handle NA's! In case of NA in Y in one of the blocks, then that entire block should be removed.
  
  
  # Loading needed packages
  if(!require(coin))
  {
    print("You are missing the package 'coin', we will now try to install it...")
    install.packages("coin")
    library(coin)
  }
  
  if(!require(multcomp))
  {
    print("You are missing the package 'multcomp', we will now try to install it...")
    install.packages("multcomp")
    library(multcomp)
  }
  
  if(!require(colorspace))
  {
    print("You are missing the package 'colorspace', we will now try to install it...")
    install.packages("colorspace")
    library(colorspace)
  }
  
  
  # get the names out of the formula
  formu.names <- all.vars(formu)
  Y.name <- formu.names[1]
  X.name <- formu.names[2]
  block.name <- formu.names[3]
  
  if(dim(data)[2] >3) data <- data[,c(Y.name,X.name,block.name)]   # In case we have a "data" data frame with more then the three columns we need. This code will clean it from them...
  
  # Note: the function doesn't handle NA's. In case of NA in one of the block T outcomes, that entire block should be removed.
  
  # stopping in case there is NA in the Y vector
  if(sum(is.na(data[,Y.name])) > 0) stop("Function stopped: This function doesn't handle NA's. In case of NA in Y in one of the blocks, then that entire block should be removed.")
  
  # make sure that the number of factors goes with the actual values present in the data:
  data[,X.name ] <- factor(data[,X.name ])
  data[,block.name ] <- factor(data[,block.name ])
  number.of.X.levels <- length(levels(data[,X.name ]))
  if(number.of.X.levels == 2) { warning(paste("'",X.name,"'", "has only two levels. Consider using paired wilcox.test instead of friedman test"))}
  
  # making the object that will hold the friedman test and the other.
  the.sym.test <- symmetry_test(formu, data = data,   ### all pairwise comparisons
                                teststat = "max",
                                xtrafo = function(Y.data) { trafo( Y.data, factor_trafo = function(x) { model.matrix(~ x - 1) %*% t(contrMat(table(x), "Tukey")) } ) },
                                ytrafo = function(Y.data){ trafo(Y.data, numeric_trafo = rank, block = data[,block.name] ) }
  )
  # if(to.print.friedman) { print(the.sym.test) }
  
  
  if(to.post.hoc.if.signif)
  {
    if(pvalue(the.sym.test) < signif.P)
    {
      # the post hoc test
      The.post.hoc.P.values <- pvalue(the.sym.test, method = "single-step")   # this is the post hoc of the friedman test
      
      
      # plotting
      if(to.plot.parallel & to.plot.boxplot)   par(mfrow = c(1,2)) # if we are plotting two plots, let's make sure we'll be able to see both
      
      if(to.plot.parallel)
      {
        X.names <- levels(data[, X.name])
        X.for.plot <- seq_along(X.names)
        plot.xlim <- c(.7 , length(X.for.plot)+.3)   # adding some spacing from both sides of the plot
        
        if(color.blocks.in.cor.plot)
        {
          blocks.col <- rainbow_hcl(length(levels(data[,block.name])))
        } else {
          blocks.col <- 1 # black
        }
        
        data2 <- data
        if(jitter.Y.in.cor.plot) {
          data2[,Y.name] <- jitter(data2[,Y.name])
          par.cor.plot.text <- "Parallel coordinates plot (with Jitter)"
        } else {
          par.cor.plot.text <- "Parallel coordinates plot"
        }
        
        # adding a Parallel coordinates plot
        matplot(as.matrix(reshape(data2,  idvar=X.name, timevar=block.name,
                                  direction="wide")[,-1])  ,
                type = "l",  lty = 1, axes = FALSE, ylab = Y.name,
                xlim = plot.xlim,
                col = blocks.col,
                main = par.cor.plot.text)
        axis(1, at = X.for.plot , labels = X.names) # plot X axis
        axis(2) # plot Y axis
        points(tapply(data[,Y.name], data[,X.name], median) ~ X.for.plot, col = "red",pch = 4, cex = 2, lwd = 5)
      }
      
      if(to.plot.boxplot)
      {
        # first we create a function to create a new Y, by substracting different combinations of X levels from each other.
        subtract.a.from.b <- function(a.b , the.data)
        {
          the.data[,a.b[2]] - the.data[,a.b[1]]
        }
        
        temp.wide <- reshape(data,  idvar=X.name, timevar=block.name,
                             direction="wide")    #[,-1]
        wide.data <- as.matrix(t(temp.wide[,-1]))
        colnames(wide.data) <- temp.wide[,1]
        
        Y.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 2, subtract.a.from.b, the.data =wide.data)
        names.b.minus.a.combos <- apply(with(data,combn(levels(data[,X.name]), 2)), 2, function(a.b) {paste(a.b[2],a.b[1],sep=" - ")})
        
        the.ylim <- range(Y.b.minus.a.combos)
        the.ylim[2] <- the.ylim[2] + max(sd(Y.b.minus.a.combos))   # adding some space for the labels
        is.signif.color <- ifelse(The.post.hoc.P.values < .05 , "green", "grey")
        
        boxplot(Y.b.minus.a.combos,
                names = names.b.minus.a.combos ,
                col = is.signif.color,
                main = "Boxplots (of the differences)",
                ylim = the.ylim
        )
        legend("topright", legend = paste(names.b.minus.a.combos, rep(" ; PostHoc P.value:", number.of.X.levels),round(The.post.hoc.P.values,5)) , fill =  is.signif.color )
        abline(h = 0, col = "blue")
        
      }
      
      list.to.return <- list(Friedman.Test = the.sym.test, PostHoc.Test = The.post.hoc.P.values)
      if(to.print.friedman) {print(list.to.return)}
      return(list.to.return)
      
    }   else {
      print("The results where not significant, There is no need for a post hoc test")
      return(the.sym.test)
    }
  }
  
  # Original credit (for linking online, to the package that performs the post hoc test) goes to "David Winsemius", see:
  # http://tolstoy.newcastle.edu.au/R/e8/help/09/10/1416.html
}



####################################################################################################
####################################################################################################
# Activate Functions
####################################################################################################
####################################################################################################


#kfold adaboost
kfold.ada(30)

#kfold logitboost
kfold.logit(30)

#kfold gentleboost
kfold.gentle(30)

#kfold flexboost 
kfold.flex(30,0.2,1)


#plot acc by col 1-29(except Iris_1 due to acc 1)
plot.result(3)


#Load acc step_1
res.acc.all <- read.csv(url("http://bit.ly/Flex_Boost"), header = T)

#Load acc step_2
res.ranks <- as.matrix(res.acc.all[, 3:6])

#Load acc step_3
for(i in 1:nrow(res.ranks)){res.ranks[i,] <- rank(-res.acc.all[i, 3:6], ties.method = "min")}

#Friedman test
res.fm.ph <- friedman.post.hoc(value ~ X2 | X1, melt(res.ranks))

#Friedman test p value matrix
res.p.val <- c(rep(NA, 4), res.fm.ph$PostHoc.Test[3], rep(NA, 3), 
               res.fm.ph$PostHoc.Test[c(2, 6)], rep(NA, 2),
               res.fm.ph$PostHoc.Test[c(1, 5, 4)], NA)





