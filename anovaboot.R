## One Way ANOVA for both wild bootstrap and residual bootstrap

#######################################################################################################################
# Added data option, distribution, test for quantitative. Warning on missing cases. Response has to be Vector.
# change to %*%. Checked k. 
#######################################################################################################################


## QUESTIONS:  Should wild bootstrap allow for different distributions than N(0,1)?
##             Want to try adding a data= option?
##             ... option that is associated with aov or lm?

### BL: Added data option. Added distribution term.

# distribution choose from "Normal", "Uniform", "LogNorm".

anova.boot <- function(formula, B=1000, type="residual", wild.dist="Normal", variance=1, seed=NULL, data = NULL, keep.data=FALSE){
  start <- proc.time()
  ##Check that a formula is input
  if(inherits(formula, "formula")==FALSE){
    stop("The input model must be a formula.")
  }
  
  if (variance <=0) {
    stop("Variance should be positive!")
  }
  
  
  full.model.frame <- model.frame(formula, data=data, na.action = na.pass)
  resp <- model.response(full.model.frame) ##get the response variable
  
  ### BL: Changed is.atomic -> is.vector, since is.atomic(NULL) -> TRUE.
  ### When someone fit ~Sepal.Length + Sepal.Width + Petal.Length. Response is NULL.
  if(is.matrix(resp)!=TRUE && is.vector(resp)!=TRUE){
    stop("Response must be a vector or matrix.")
  }
  else if((dim(resp)[1]==0 || dim(resp)[2]==0) && length(resp)==0){
    stop("Response must have entries.")
  }
  else if(mode(resp)!="numeric"){
    stop("Response must be of type numeric.")
  }
  else if(anyNA(resp)==TRUE){
    stop("Response must not have any missing values.")
  }
  #}
  
  
  ##This piece is not complete
  
  ##QUESTION:  Does it account for . case?
  ##COMMENT:  p number of factors (not accounting for levels)
  
  ##Check that the covariate variables are in the correct format
  if(length(rhs.vars(formula))<1){
    stop("Linear model must have at least 1 predictor variable.")
  }
  
  p <- length(rhs.vars(formula)) #number of covariates, does not include intercept -- does not account for matrices...
  
  if (p > 2) {
    stop("This function only work for 1 or 2 predictors.")
  }
  
  modelMat <- model.matrix(formula, data=data)[,]
  modelqr <- qr(modelMat)
  if (ncol(modelMat) > modelqr$rank) {
    warning("The design matrix isn't full column rank.")
  }
  
  if (modelqr$rank == 1) {
    stop("Model is same as response~1. Can't do ANOVA test")
  }
  model.pivot <- modelqr$pivot[1:modelqr$rank]
  
  ##Check that the type of bootstrap is implementable
  if(type!="residual" & type != "wild"){
    stop("Only residual bootstrap or wild boostrap is allowed for type.")
  } 
  
  n <- length(resp)

  
  if(mode(B)!="numeric"){
    stop("Number of bootstrap samples, B, must be of type numeric.")
  }
  else if(is.atomic(B)!=TRUE){
    stop("Number of bootstrap samples, B, must be a constant.")
  }
  else if( B < n){
    # Change this to a warning instead
    warning("Number of bootstrap samples is recommended to be more than the data length.")
  }
  
  
  ##Check the the seed is numeric, or set a seed if unspecified
  if(is.null(seed)==TRUE){
    seed <- sample(seq(1,100000000), size=1)
  }
  else{
    if(mode(seed)!="numeric"){
      stop("The seed must be of type numeric.")
    }
    else if(is.atomic(seed)!=TRUE){
      stop("The seed must be a constant.")
    }
  }
  set.seed(seed)
  
  ##################################################
  ## Least Squares Fit
  ##  -- update to use Type I SS...
  ##################################################
  anovaMod <- anova(lm(formula, data = data))
  fullModMat <- model.matrix(formula, data=data)[,]
  SS.list <- anovaMod$'Mean Sq'
  df.list <- anovaMod$Df
  termNames <- row.names(anovaMod) 
  df.len <- length(df.list)
  modelMat.list <- list()
  modelMat.list[[1]] <- matrix(fullModMat[, 1], ncol=1)
  count <- 1
  for (i in 1:(df.len-1)) {
    count <- count + df.list[i]
    modelMat.list[[i+1]] <- fullModMat[, model.pivot[1:count]]
  }
  
  obsDataregFit <- aov(formula, data = data)          #fit the linear model specified in formula input
  resp <- obsDataregFit$model[,1]                     #get the response variable
  #predictor <- obsDataregFit$model[,2]
  estParam <- matrix(obsDataregFit$coefficients, ncol=1)      #keep the param. estimates in a vector
  #obsDataResid <- as.vector(residuals(obsDataregFit)) #keep the original residuals
  ParamNames <- names(obsDataregFit$coefficients)     #keep the coefficient name/association
  rownames(estParam) <- ParamNames                    #name the rows for the parameters so we know what they are
  
  obsFStats <- summary(obsDataregFit)[[1]][,4]        #get the F test stats from orig fit
  obsFStats <- obsFStats[-length(obsFStats)]          #remove NA next to MSE
  respLen <- length(resp)
  #complete.model.frame <- model.frame(formula, data=data)
  
  
  ##QUESTION:  I'm a bit worried about this code translating to the 2-way ANOVA case
  ##           Does this work when a formula contains -1, .?
  ##           data= not part of this function yet, but what happens when we add this flexibility?
  ##           the modelMat variable is a list!  It contains more than just the matrix of predictor variables
  
  ### BL: Not test yet but I think the formula will work fine with -1 and ., 
  ### since we use model.frame and model.matrix to utilize formula.
  
  ##COMMENT:  I'm not convinced that the k variable works (when trying some pseudo code)
  
  ### BL: Checked k and k works for me. k = 0, means the predictor is quantitative.
  ### Also if k won't work, the tests below won't have the same results as the anova.boot function online. 
  
  ##COMMENT:  Need to check there aren't missing values in the predictor variables?  We check the response
  ##          doesn't contain missing values, but not predictors.  should see respLen=n
  
  ### BL: No, respLen won't equal to n if there are missing value in predictor variables. aov() will remove those cases.
  ### Thus, we will know how many cases removed by n-respLen
  
  if (respLen != n) {
    warning(paste("There are", n - respLen, "missing cases removed!"))
  }
  
  
  if (keep.data == TRUE) {
    bootResponseMatrix.list <- list()
  }
  
  bootFStats <- matrix(NA, nrow=B, ncol=(df.len-1))   #bootstrap F statistics
  
  bootSSEMat <- matrix(NA, nrow=B, ncol=1)
  bootSSTrMat <- matrix(NA, nrow=B, ncol=(df.len-1))
  
  for(i in 1:(df.len-1)) {
    nullMat <- modelMat.list[[i]]
    nullQuadMat <- solve(t(nullMat) %*% nullMat) %*% t(nullMat)
    y.orighat <- nullMat %*% (nullQuadMat %*% resp)
    
    modelMat <- modelMat.list[[i+1]]                    #model matrix (X)
    quadMat <- solve(t(modelMat) %*% modelMat) %*% t(modelMat) #projection matrix (X^TX)^-1 X^T
    y.hat <- modelMat %*% (quadMat %*% resp)
    obsDataResid <- resp - y.hat
    
    if (keep.data == TRUE) {
      bootResponseMatrix.list[[i]] <- matrix(NA, nrow=B, ncol=respLen)
    }
    
    df.Tr <- df.list[i]
    df.E <- respLen - (1 + sum(df.list[1:i]))
    for(j in 1:B){
      if (type == "residual") {
        ##QUESTION:  Is there a way/would it be better to remove the loop structure and do residual bootstrap with
        ##           multinomials?
        bootResid <- matrix(sample(obsDataResid, replace=TRUE), ncol=1)  #bootstrap residuals
        
        #Bootstrap under H0 to get F test statistics
        yb <- matrix(y.orighat + bootResid, nrow=respLen, ncol=1)      #bootstrap response under H0
      } else {
        if (wild.dist == "Normal") {
          bootResid <- matrix(obsDataResid * rnorm(respLen, mean=0, sd=sqrt(variance)), ncol=1)  
        } else if (wild.dist == "Uniform") {
          bootResid <- matrix(obsDataResid * runif(respLen, -sqrt(3*variance), sqrt(3*variance)), ncol=1)  
        } else if (wild.dist == "LogNorm") {
          bootResid <- matrix(obsDataResid * rlnorm(respLen, meanlog=0, sdlog=sqrt(log(1/2 + sqrt(1+4*variance)/2))) - 
                               exp(log(1/2 + sqrt(1+4*variance)/2)/2), ncol=1)  
        } else {
          stop("Distribution can only choose from Normal, Uniform, LogNorm. Default: Normal, with variance: 1")
        }
          #Bootstrap under H0 to get F test statistics
          
          ##QUESTION:  modelMat %*% beta.hat can come outside the loop to save time?
          ##           same question about crossprod()
          ### BL: Fixed crossprod(), Done modelMat %*% beta.hat.
          ##           Should we add an option to use different distributions for wild bootstrap?
          ### Added some distributions: Normal, Uniform, LogNorm.
        yb <- matrix(y.orighat + bootResid, nrow=respLen, ncol=1)      #bootstrap response under H0
      }
      if (keep.data == TRUE) {
        bootResponseMatrix.list[[i]][j, ] <- yb
      }
      y.null <- nullMat %*% (nullQuadMat %*% yb)
      
      ##QUESTION:  I'm guessing that using the crossprod() function is slower than matrix mult with %*%
      ### BL: Fixed crossprod()
      y.hat <- modelMat %*% (quadMat %*% yb)
      SSE <- t(yb - y.hat) %*% (yb - y.hat)
      SST <- t(yb - y.null) %*% (yb - y.null)
      SSTr <- SST - SSE
      MSE <- SSE / df.E
      MSTr <- SSTr / df.Tr
      Ftest <- MSTr/MSE
      bootSSTrMat[j, i] <- SSTr
      if (i == (df.len-1)){
        bootSSEMat[j, 1] <- SSE
      }
      bootFStats[j, i]<- Ftest
    }
 
    ##COMMENT:  I think the stop above is already accounted for in the front in function input checks
    ### BL: Checked already, so deleted.
    
  }
  
  pvalues <- rep(NA, df.len-1)
  for (i in 1:(df.len-1)) {
    pvalues[i] <- mean(bootFStats[,i] > obsFStats[i])
  }
  
  end <- proc.time()
  print(end-start)
  if (keep.data == TRUE) {
    structure(invisible(list(terms = termNames,df = df.list, bootFStats=bootFStats, 
                             origSSE=SS.list[length(SS.list)], origSSTr=SS.list[1:(length(SS.list)-1)],
                             bootSSE=bootSSEMat, bootSSTr=bootSSTrMat,
                             origFStats=obsFStats, "p-values" = pvalues,
                             bootResponse=bootResponseMatrix.list)))
  } else {
    structure(invisible(list(terms = termNames, df= df.list, bootFStats=bootFStats,
                             origSSE=SS.list[length(SS.list)], origSSTr=SS.list[1:(length(SS.list)-1)],
                             bootSSE=bootSSEMat, bootSSTr=bootSSTrMat,
                             origFStats=obsFStats, "p-values" = pvalues)))
  }
}
