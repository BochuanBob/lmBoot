myTwoWayANOVA <- function(predictorA, predictorB, response) {
  start <- proc.time()
  flag <- (!is.na(response)) & (!is.na(predictorA)) & (!is.na(predictorB))
  predictorA <- predictorA[flag]
  predictorB <- predictorB[flag]
  response <- response[flag]
  
  KA <- length(levels(predictorA))
  KB <- length(levels(predictorB))
  n <- length(response)
  
  if (KA == 0) {
    KA <- 2
  }
  if (KB == 0) {
    KB <- 2
  }
  XAmultB <- model.matrix(response~predictorA*predictorB)[,]
  XA <- model.matrix(response~predictorA)[,]
  XAplusB <- model.matrix(response~predictorA+predictorB)[,]
  
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
  
  end <- proc.time()
  print(end-start)
  return(result)
}


# Test 1
speed <- rollarCoasterData$Speed
construct <- rollarCoasterData$Construction
region <- rollarCoasterData$Region
anova(lm(Speed~Construction*Region,data=rollarCoasterData))
myTwoWayANOVA(construct,region,speed)

a <- anova(lm(mpg~factor(gear)*factor(vs),data=mtcars))$F.stats
b <- myTwoWayANOVA(factor(mtcars$gear),factor(mtcars$vs), mtcars$mpg)

(a$`F value`[1:3]-b$F0)

library(car)
Anova(lm(Speed~Region*Construction,data=rollarCoasterData), type='III')

# Test 2
hours <- c(130,155,34,40,20,70,74,180,80,75,82,58,150,188,136,122,25,70,159,126,106,115,58,45,138,110,174,120,96,104,168,160,150,139,82,60)
temp <- as.factor(rep(c("15","15","70","70","125","125"),6))
material <- as.factor(rep(c("1","2","3"), each=12))
anova(lm(hours~material*temp))
myTwoWayANOVA(material, temp, hours)

