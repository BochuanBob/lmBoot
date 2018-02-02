
PMData <- read.csv("C:/Users/lub/Documents/R_Courses/MA382/PMdata.csv")
PmMod <- lm(PMData$PM2.5~PMData$City)
anova(PmMod)

city <- PMData$City
pm <- PMData$PM2.5

myOneWayANOVA <- function(predictor, response) {
  flag <- (!is.na(response)) & (!is.na(predictor))
  predictor <- predictor[flag]
  response <- response[flag]
  n.levels <- levels(predictor)
  k <- length(n.levels)
  n <- length(response)
  
  num_list <- rep(0, k)
  xbar_list <- rep(0, k)
  se2_list <- rep(0, k)
  for (i in 1:k) {
    xk <- response[predictor == n.levels[i]]
    num_list[i] <- length(xk)
    xbar_list[i] <- mean(xk)
    se2_list[i] <- var(xk)
  }
  
  xbar <- mean(response)
  
  SSTr <- sum(num_list * (xbar_list - xbar)^2)
  SSE <- sum((num_list-1) * se2_list)
  
  MSE <- SSE / (n-k)
  MSTr <- SSTr / (k-1)
  
  Ftest <- MSTr/MSE
  
  return(Ftest)
}

myOneWayANOVA(city, pm)

myTwoWayANOVA <- function(predictorA, predictorB, response) {
  flag <- (!is.na(response)) & (!is.na(predictorA)) & (!is.na(predictorB))
  predictorA <- predictorA[flag]
  predictorB <- predictorB[flag]
  response <- response[flag]
  A.levels <- levels(predictorA)
  B.levels <- levels(predictorB)
  kA <- length(A.levels)
  kB <- length(B.levels)
  
  ybar... <- mean(response)
  ybari.. <- rep(0, kA)
  ybar.j. <- rep(0, kB)
  ybarij. <- matrix(0, nrow=kA, ncol=kB)
  
  n... <- length(response)
  ni.. <- rep(0,kA)
  n.j. <- rep(0,kB)
  nij. <- matrix(0, nrow=kA, ncol=kB)
  
  sigij <- matrix(0, nrow=kA, ncol=kB)
  for (i in 1:kA) {
    yi.. <- response[predictorA == A.levels[i]]
    ybari..[i] <- mean(yi..)
    ni..[i] <- length(yi..)
  }
  
  for (j in 1:kB) {
    y.j. <- response[predictorB == B.levels[j]]
    ybar.j.[j] <- mean(y.j.)
    n.j.[j] <- length(y.j.)
  }
  
  SSE <- 0
  for (i in 1:kA) {
    for (j in 1:kB) {
      yij. <- response[(predictorA == A.levels[i]) & (predictorB == B.levels[j])]
      ybarij.[i,j] <- mean(yij.)
      nij.[i,j] <- length(yij.)
      sigij[i,j] <- sd(yij.)
      SSE <- SSE + sum((yij. - ybarij.[i,j])^2)
    }
  }
  
  SST <- sum((response-ybar...)^2)
  
  SSA <- sum((ybari.. - ybar...)^2 * ni..)
  
  SSB <- sum((ybar.j. - ybar...)^2 * n.j.)
  
  #SSAB <- SST - SSA - SSB - SSE
  SSAB <- sum((ybarij. - matrix(ybari.., nrow=kA, ncol=kB, byrow=FALSE) - matrix(ybar.j., nrow=kA, ncol=kB, byrow=TRUE) + ybar...)^2 * nij./(sigij^2))
  
  MSA <- SSA/(kA-1)
  MSB <- SSB/(kB-1)
  MSAB <- SSAB/((kA-1)*(kB-1))
  MSE <- SSE/(n...-kA*kB)
  
  return(list(SST, SSA, SSB, SSAB, SSE, MSA/MSE, MSB/MSE, MSAB/MSE))
}


constructXMatrix <- function(predictorA, predictorB, KA, KB, n) {
  A.levels <- levels(predictorA)
  B.levels <- levels(predictorB)
  Xtotal <- matrix(0, ncol=(KA*KB), nrow=n)
  Xtotal[ , 1] <- 1
  for (i in 1:n) {
    alevel <- which(predictorA[i]==levels(predictorA))-1
    blevel <- which(predictorB[i]==levels(predictorB))-1
    if (alevel > 0) {
      Xtotal[i, 1+alevel] = 1
    }
    if (blevel > 0) {
      Xtotal[i, KA+blevel] = 1
    }
    if (alevel*blevel > 0) {
      Xtotal[i, KA+KB-1+alevel+(blevel-1)*(KA-1)] = 1
    }
  }
  XA <- Xtotal[ , 1:KA]
  XAplusB <- Xtotal[ , 1:(KA+KB-1)]
  return(list(XA,XAplusB,Xtotal))
}

myTwoWayANOVA2 <- function(predictorA, predictorB, response) {
  flag <- (!is.na(response)) & (!is.na(predictorA)) & (!is.na(predictorB))
  predictorA <- predictorA[flag]
  predictorB <- predictorB[flag]
  response <- response[flag]
  modA <- lm(response~predictorA)
  modAplusB <- lm(response~predictorA+predictorB)
  modAmultB <- lm(response~predictorA*predictorB)
  
  KA <- length(levels(predictorA))
  KB <- length(levels(predictorB))
  n <- length(response)
  
  X.list <- constructXMatrix(predictorA, predictorB, KA, KB, n)
  XA <- X.list[[1]]
  XAplusB <- X.list[[2]]
  XAmultB <- X.list[[3]]
  Y <- matrix(response, nrow=length(response), ncol=1)
  Ybar <- mean(response)
  YA.hat <- XA %*% solve(t(XA) %*% XA) %*% t(XA) %*% Y
  YAplusB.hat <- XAplusB %*% solve(t(XAplusB) %*% XAplusB) %*% t(XAplusB) %*% Y
  YAmultB.hat <- XAmultB %*% solve(t(XAmultB) %*% XAmultB) %*% t(XAmultB) %*% Y
  SSEA <- t(Y - YA.hat) %*% (Y - YA.hat)
  SST <- t(Y - Ybar) %*% (Y - Ybar)
  SSEAplusB <- t(Y - YAplusB.hat) %*% (Y - YAplusB.hat)
  SSE <- t(Y - YAmultB.hat) %*% (Y - YAmultB.hat)
  
  SSA <- SST-SSEA
  SSB <- SSEA - SSEAplusB
  SSAB <- SSEAplusB - SSE
  
  MSA <- SSA/(KA-1)
  MSB <- SSB/(KB-1)
  MSAB <- SSAB/((KA-1)*(KB-1))
  MSE <- SSE/(n-KA*KB)
  
  sumS <- c(SSA, SSB,SSAB, SSE)
  df <- c(KA-1, KB-1, (KA-1)*(KB-1), n-KA*KB)
  meanS <- c(MSA, MSB, MSAB, MSE)
  F0 <- c(MSA/MSE, MSB/MSE, MSAB/MSE)
  
  result <- list(df, sumS,meanS, F0)
  names(result) <- c("Degrees of Freedom", "Sum of Squares", "Mean Square", "F0")
  
  return(result)
}

# Test 1
speed <- rollarCoasterData$Speed
construct <- rollarCoasterData$Construction
region <- rollarCoasterData$Region
anova(lm(Speed~Construction*Region,data=rollarCoasterData))
myTwoWayANOVA(construct,region,speed)

library(car)
Anova(lm(Speed~Region*Construction,data=rollarCoasterData), type='III')

# Test 2
hours <- c(130,155,34,40,20,70,74,180,80,75,82,58,150,188,136,122,25,70,159,126,106,115,58,45,138,110,174,120,96,104,168,160,150,139,82,60)
temp <- as.factor(rep(c("15","15","70","70","125","125"),6))
material <- as.factor(rep(c("1","2","3"), each=12))
anova(lm(hours~material*temp))
myTwoWayANOVA(material, temp, hours)
