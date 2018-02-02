## Testing
library(bitops)
library(formula.tools)
library(RCurl)
library(stats)
data("iris")
data(mtcars)
s4 <- getURL("https://raw.githubusercontent.com/elegacy/lm.boot/master/anova.boot.R",
             ssl.verifypeer=FALSE)
eval(parse(text=s4))
remove(s4)
anova.boot.old <- anova.boot
source("C:/Users/lub/Documents/R_Courses/anovaGitHub/lmBoot/anovaboot.R")

### Test 1 Basic one (Work Well)
formula1 <- iris$Sepal.Length~iris$Species
bootres.new <- anova.boot(formula1, B =1000, seed =27, keep.data=TRUE)

# Test for FTest same as using anova()
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[1]][i,]~iris$Species))$'F value'[1] - bootres.new$bootFStats[i,1]) < 10^(-10))
}

hist(bootres.new$bootFStats[, 1])
# Make sure the result is the same using same seed.
sum(abs(bootres.new$bootFStats-bootres.old$bootFStats) <= 10^-6) == 1000
bootres.new$origFStats == bootres.old$origFStats
#bootres.new$origEstParam == bootres.old$origEstParam

### BL: I think the bootEstParam in the online version of anova.boot is incorrect since I have compare the new verson with aov()$coefficients 
#sum(bootres.new$bootEstParam == bootres.old$bootEstParam) == dim(bootres.old$bootEstParam)[1] * dim(bootres.old$bootEstParam)[2]

### Test 2 Quantitative Predictor (Work Well)
formula2 <- iris$Sepal.Length~iris$Sepal.Width
bootres.new <- anova.boot(formula2, B =100, seed =27)
bootres.old <- anova.boot.old(formula2, B =100, seed =27)

# Make sure the result is the same using same seed.
sum(abs(bootres.new$bootFStats-bootres.old$bootFStats) <= 10^-6) == 100
bootres.new$origFStats == bootres.old$origFStats
#bootres.new$origEstParam == bootres.old$origEstParam

### Test 3 data option (Work Well)
formula3 <- Sepal.Length~Species
bootres.new <- anova.boot(formula3, data=iris, B =1000, seed =27)
bootres.old <- anova.boot.old(formula1, B =1000, seed =27)
# Make sure the result is the same using same seed.
sum(abs(bootres.new$bootFStats-bootres.old$bootFStats) <= 10^-6) == 1000
bootres.new$origFStats == bootres.old$origFStats
#bootres.new$origEstParam == bootres.old$origEstParam


### Test 4 Missing values. Checked how many missing cases for predictors. (Same as cases removed, since response can't have any missing cases)
### Work well
myData <- iris
myData[1:10,-1] <- NA
bootres.new <- anova.boot(formula3, data=myData, B =1000, seed =27)


### Test 5 bootstrap response (Work Well)
bootres.new1 <- anova.boot(formula1, B =1000, seed =27, keep.data=TRUE)
bootres.new2 <- anova.boot(formula1, B =1000, seed =27, keep.data=TRUE, type = "wild")
dim(bootres.new1$bootResponse[[1]])
dim(bootres.new2$bootResponse[[1]])


### Test 6 two way ANOVA bootstrap *
formula4 <- mpg~factor(vs)*wt
bootres.new <- anova.boot(formula4, data=mtcars, B=1000, seed=27, keep.data = TRUE)
dim(bootres.new$bootFStats)
bootres.new$origSSE
bootres.new$origSSTr
length(bootres.new$bootSSE)
dim(bootres.new$bootSSTr)
# Test for FTest same as using anova()
index <- 0
for (i in 1:1000) {
  index <- index + (abs(anova(lm(bootres.new$bootResponse[[1]][i,]~factor(mtcars$vs)))$'F value'[1] - bootres.new$bootFStats[i,1]) < 10^(-10))
}
index
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[2]][i,]~factor(vs)+wt, data=mtcars))$'F value'[2] - bootres.new$bootFStats[i,2]) < 10^(-10))
}
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[3]][i,]~factor(vs)*wt, data=mtcars))$'F value'[3] - bootres.new$bootFStats[i,3]) < 10^(-10))
}

formula5 <- mpg~factor(gear)*factor(vs)
bootres.new <- anova.boot(formula5, data=mtcars, B=1000, seed=27, keep.data = TRUE)
dim(bootres.new$bootFStats)

# Test for FTest same as using anova()
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[1]][i,]~factor(mtcars$gear)))$'F value'[1] - bootres.new$bootFStats[i,1]) < 10^(-10))
}
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[2]][i,]~factor(gear)+factor(vs), data=mtcars))$'F value'[2] - bootres.new$bootFStats[i,2]) < 10^(-10))
}
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[3]][i,]~factor(gear)*factor(vs), data=mtcars))$'F value'[3] - bootres.new$bootFStats[i,3]) < 10^(-10))
}

### Test 7 two way ANOVA bootstrap +
formula6 <- mpg~factor(gear)+factor(vs)
bootres.new <- anova.boot(formula6, data=mtcars, B=1000, seed=27, keep.data = TRUE)
dim(bootres.new$bootFStats)
bootres.new$`p-values`
# Test for FTest same as using anova()
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[1]][i,]~factor(mtcars$gear)))$'F value'[1] - bootres.new$bootFStats[i,1]) < 10^(-10))
}
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[2]][i,]~factor(gear)+factor(vs), data=mtcars))$'F value'[2] - bootres.new$bootFStats[i,2]) < 10^(-10))
}



### Test 8 two way ANOVA bootstrap * (deal with exact singular cases)
formula7 <- mpg~factor(cyl)*factor(vs)
bootres.new <- anova.boot(formula7, data=mtcars, B=1000, seed=27, keep.data = TRUE)
bootres.new$`p-values`
par(mfrow=c(2,2))
hist(bootres.new$bootFStats[,1])
hist(bootres.new$bootFStats[,2])
hist(bootres.new$bootFStats[,3])
# Test for FTest same as using anova()
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[1]][i,]~factor(mtcars$cyl)))$'F value'[1] - bootres.new$bootFStats[i,1]) < 10^(-10))
}
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[2]][i,]~factor(cyl)+factor(vs), data=mtcars))$'F value'[2] - bootres.new$bootFStats[i,2]) < 10^(-10))
}
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[3]][i,]~factor(cyl)*factor(vs), data=mtcars))$'F value'[3] - bootres.new$bootFStats[i,3]) < 10^(-10))
}


### Test 9
formula8 <- Sepal.Length~Species+paste(Species, "123")
bootres.new <- anova.boot(formula8, data=iris, B=1000, seed=27, keep.data = TRUE)
for (i in 1:10) {
  print(abs(anova(lm(bootres.new$bootResponse[[1]][i,]~iris$Species))$'F value'[1] - bootres.new$bootFStats[i,1]) < 10^(-10))
}

##CHECK:  Do both bootstrap results match aov() output in estimates, F statistic, etc?
##        What if the rhs contains a quantitative variable?  Does it break?

### BL: Works and matches in both categorical and quantitative predictor cases.

##QUESTION:  Do we want to keep the estimated parameters?  the lm functions already do this?

### BL: Make decision it later.

##CHECK:  Try on more examples, examples that shouldn't work, etc.

### BL: Added some tests.

##QUESTION:  Are there other things we might want to keep from the bootstrap?  Maybe what the seed value
##           was set at?  Do we want an option to keep bootstrap samples (that the user may turn on)?
##           Do we want to keep any other pieces of the ANOVA table?

### BL: Make decision later.

# Wild bootstrap
### Should have no warnings and errors.
bootres.new.wild1 <- anova.boot(iris$Sepal.Length~iris$Species, B=1000, type = "wild", wild.dist = "Uniform", seed=27)
bootres.new.wild2 <- anova.boot(iris$Sepal.Length~iris$Species, B=1000, type = "wild", wild.dist = "Uniform", variance=10, seed=27)
bootres.new.wild3 <- anova.boot(iris$Sepal.Length~iris$Species, B=1000, type = "wild", wild.dist = "LogNorm", seed=27)

### Should get some errors.
# Variance should be positive
bootres.new.wild4 <- anova.boot(iris$Sepal.Length~iris$Species, B=1000, type = "wild", wild.dist = "Uniform", variance=-1, seed=27)
# The name of wild.dist
bootres.new.wild5 <- anova.boot(iris$Sepal.Length~iris$Species, B=1000, type = "wild", wild.dist = "wrongthing",seed=27)