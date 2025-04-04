rm(list=ls()) #clear

## (1) Neyman-based inference; binary treatment and response

## Data

x <- c(rep(0,20), rep(1,20))
y <- c(0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,
       1,1,1,1,0,0,0,1,1,0,0,1,1,0,0,1,1,0,1,1)
n <- length(x)
k <- sum(x)  # size of treatment group

mean.ctrl <- mean(y[x==0])
mean.trt <- mean(y[x==1])
hat.ace <-  mean.trt - mean.ctrl
hat.ace


## Confidence interval:
var.ctrl <- var(y[x==0])
var.trt <- var(y[x==1])

hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)

##95% asymptotic confidence interval:

hat.ace - 1.96*sqrt(hat.var.hat.ace)
hat.ace + 1.96*sqrt(hat.var.hat.ace)

#######################################################

# (2) Neyman-based inference; binary treatment and 
#     continuous response

x <- c(rep(0,20), rep(1,20))
y <- c(9.87, 12.14,9.62,8.63,11.40,7.40,6.88,10.00,7.39,11.31,
       11.56, 8.78, 3.70, 10.36, 10.68, 6.60, 11.11, 7.94, 10.97, 10.43,
       14.90, 17.06, 20.62, 18.18, 15.29, 12.38, 18.27, 16.00, 13.45, 16.29,
       14.93, 11.35, 17.02, 16.60, 13.99, 16.59, 10.76, 16.79, 13.42, 15.30)
n <- length(x)
k <- sum(x)

mean.ctrl <- mean(y[x==0])
mean.trt <- mean(y[x==1])
hat.ace <-  mean.trt - mean.ctrl
hat.ace


## Confidence interval:
var.ctrl <- var(y[x==0])
var.trt <- var(y[x==1])

hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)
hat.var.hat.ace

##95% asymptotic confidence interval:

hat.ace - 1.96*sqrt(hat.var.hat.ace)
hat.ace + 1.96*sqrt(hat.var.hat.ace)




#######################################################

# (3) Neyman-based inference; binary treatment and 
#     continuous response
#     Investigating coverage properties and additivity


n <- 40  #finite population size
k <- 20  #treatment group size


## (3a) Simulating potential outcomes under additivity
ace.true <- 3
yctrl <- rnorm(n,10,2) # simulating Y(0) [from N(10,2)]
ytrt <- yctrl + ace.true  # additivity

# summary and histogram of ICE
summary(ytrt-yctrl)
hist(ytrt-yctrl)

mean(ytrt) - mean(yctrl) # true ACE
nsims <- 5000
covered <- rep(NA,nsims) # vector for coverage indicators
hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec <- rep(NA,nsims)

for(i in 1:nsims){
	trt.grp <- sample((1:n),k)  #choose our treatment group
	x <- rep(0,n)
	x[trt.grp] <- 1  #construct treatment vector
	y <- yctrl*(x==0) + ytrt*(x==1)  #build observed data vector
### Neyman procedure:
	mean.ctrl <- mean(y[x==0])
    mean.trt <- mean(y[x==1])
    hat.ace <-  mean.trt - mean.ctrl
    hat.ace.vec[i] <- hat.ace	# store the ACE estimate
	var.ctrl <- var(y[x==0])  # variance of control grp
    var.trt <- var(y[x==1])   # varianc
    hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)
    covered[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace)) < ace.true) &
                   ((hat.ace + 1.96*sqrt(hat.var.hat.ace)) > ace.true))
    hat.var.hat.ace.vec[i] <- hat.var.hat.ace # store var(hat ACE) estimate
}
sum(covered)/nsims  ## about 95%
mean(hat.ace.vec)  # are we unbiased ?
var(hat.ace.vec)   # true variability
mean(hat.var.hat.ace.vec) # are we unbiased?


## Histogram of the estimate of the ACE
hist(hat.ace.vec,main="Distribution of hat.ace",xlab="hat.ace")
abline(v=ace.true,col="red") #True ACE
abline(v=mean(hat.ace.vec),col="blue") #mean of histogram

## Boxplot of the estimate of the ACE
boxplot(hat.ace.vec,main="Distribution of hat.ace",xlab="hat.ace")
abline(h=ace.true,col="red") #True ACE
abline(h=mean(hat.ace.vec),col="blue") #mean of histogram

## Histogram of the estimate of the variance of the ACE
hist(hat.var.hat.ace.vec,main="Distribution of hat.var.hat.ace",
xlab="hat.var.hat.ace.vec")
abline(v=var(hat.ace.vec),col="red") #True variability of hat.ace
abline(v=mean(hat.var.hat.ace.vec),col="blue") #mean of histogram

## Boxplot of the estimate of the variance of the ACE
boxplot(hat.var.hat.ace.vec,main="Distribution of hat.var.hat.ace",
xlab="hat.var.hat.ace.vec")
abline(h=var(hat.ace.vec),col="red") #True variability of hat.ace
abline(h=mean(hat.var.hat.ace.vec),col="blue") #mean of histogram

## Scatterplot of the estimate of ACE vs. estimate of the variance of the ACE
plot(hat.ace.vec,hat.var.hat.ace.vec, main="hat(ACE) vs. hat(var(hat(ACE)))",
ylab="hat.var.hat.ace", xlab="hat.ace" )
abline(h=mean(hat.var.hat.ace.vec),col="blue") #mean of estimates of variance
abline(h=var(hat.ace.vec),col="red") #True variability of hat.ace
abline(v=mean(hat.ace.vec), col="blue") #Mean of hat.ACE
abline(v=ace.true,col="red") #True ACE


## (3b) Simulating potential outcomes under non-additivity

ace.true <- 3
yctrl <- rnorm(n,10,2) # simulating from N(10,2)
heterog <- rnorm(n,0,3) # building disturbance vector
heterog <- heterog - mean(heterog) # ensuring sums to 0
ytrt <- yctrl + ace.true + heterog # building ytrt

# summary and histogram of ICE
summary(ytrt-yctrl)
hist(ytrt-yctrl)

mean(ytrt) - mean(yctrl) # True ACE


### simulation over different random assignments
nsims <- 5000
covered <- rep(NA,nsims) # vector for coverage indicators
hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec <- rep(NA,nsims)



for(i in 1:nsims){
	trt.grp <- sample((1:n),k)  #choose our treatment group
	x <- rep(0,n)
	x[trt.grp] <- 1  #construct treatment vector
	y <- yctrl*(x==0) + ytrt*(x==1)  #build observed data vector
### Neyman procedure:
	mean.ctrl <- mean(y[x==0])
    mean.trt <- mean(y[x==1])
    hat.ace <-  mean.trt - mean.ctrl
    hat.ace.vec[i] <- hat.ace	# store the ACE estimate
	var.ctrl <- var(y[x==0])  # variance of control grp
    var.trt <- var(y[x==1])   # varianc
    hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)
    covered[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace)) < ace.true) &
                   ((hat.ace + 1.96*sqrt(hat.var.hat.ace)) > ace.true))
    hat.var.hat.ace.vec[i] <- hat.var.hat.ace # store var(hat ACE) estimate
}
sum(covered)/nsims  ## 
mean(hat.ace.vec)  # are we unbiased ?
var(hat.ace.vec)   # true variability
mean(hat.var.hat.ace.vec) # are we unbiased?

## Histogram of the estimate of the ACE
hist(hat.ace.vec,main="Distribution of hat.ace",xlab="hat.ace")
abline(v=ace.true,col="red") #True ACE
abline(v=mean(hat.ace.vec),col="blue") #mean of histogram

## Boxplot of the estimate of the ACE
boxplot(hat.ace.vec,main="Distribution of hat.ace",xlab="hat.ace")
abline(h=ace.true,col="red") #True ACE
abline(h=mean(hat.ace.vec),col="blue") #mean of histogram

## Histogram of the estimate of the variance of the ACE
hist(hat.var.hat.ace.vec,main="Distribution of hat.var.hat.ace",
xlab="hat.var.hat.ace.vec")
abline(v=var(hat.ace.vec),col="red") #True variability of hat.ace
# Recall true variability of hat.ace (red line) given by:
var(hat.ace.vec)
abline(v=mean(hat.var.hat.ace.vec),col="blue") #mean of histogram

## Boxplot of the estimate of the variance of the ACE
boxplot(hat.var.hat.ace.vec,main="Distribution of hat.var.hat.ace",
xlab="hat.var.hat.ace.vec")
abline(h=var(hat.ace.vec),col="red") #True variability of hat.ace
abline(h=mean(hat.var.hat.ace.vec),col="blue") #mean of histogram

## Scatterplot of the estimate of ACE vs. estimate of the variance of the ACE
plot(hat.ace.vec,hat.var.hat.ace.vec, main="hat(ACE) vs. hat(var(hat(ACE)))",
ylab="hat.var.hat.ace", xlab="hat.ace" )
abline(h=mean(hat.var.hat.ace.vec),col="blue") #mean of estimates of variance
abline(h=var(hat.ace.vec),col="red") #True variability of hat.ace
abline(v=mean(hat.ace.vec), col="blue") #Mean of hat.ACE
abline(v=ace.true,col="red") #True ACE


## (3c) Simulating potential outcomes: Extreme scenario

ace.true <- 3
yctrl <- rnorm(n,0,2) # simulating from N(0,2)
yctrl <- yctrl - mean(yctrl)
ytrt <- ace.true - yctrl  # building ytrt
mean(ytrt) - mean(yctrl) # True ACE

summary(ytrt-yctrl) #ICE distribution non-degenerate
hist(ytrt-yctrl)
summary((ytrt+yctrl)/2) #Average of yi(1) and yi(0) is constant!
						# (Weird simulated example.)

nsims <- 5000
covered <- rep(NA,nsims) # vector for coverage indicators
hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec <- rep(NA,nsims)

for(i in 1:nsims){
	trt.grp <- sample((1:n),k)  #choose our treatment group
	x <- rep(0,n)
	x[trt.grp] <- 1  #construct treatment vector
	y <- yctrl*(x==0) + ytrt*(x==1)  #build observed data vector
### Neyman procedure:
	mean.ctrl <- mean(y[x==0])
    mean.trt <- mean(y[x==1])
    hat.ace <-  mean.trt - mean.ctrl
    hat.ace.vec[i] <- hat.ace	# store the ACE estimate
	var.ctrl <- var(y[x==0])  # variance of control grp
    var.trt <- var(y[x==1])   # varianc
    hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)
    covered[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace)) < ace.true) &
                   ((hat.ace + 1.96*sqrt(hat.var.hat.ace)) > ace.true))
    hat.var.hat.ace.vec[i] <- hat.var.hat.ace # store var(hat ACE) estimate
}
sum(covered)/nsims  ##  100% coverage!
mean(hat.ace.vec)  # are we unbiased ?
summary(hat.ace.vec)
var(hat.ace.vec)   # true variability = 0 (!)
mean(hat.var.hat.ace.vec) # are we unbiased?

## Histogram of the estimate of the ACE (degenerate!)
hist(hat.ace.vec,main="Distribution of hat.ace",xlab="hat.ace")
abline(v=ace.true,col="red") #True ACE
abline(v=mean(hat.ace.vec),col="blue") #mean of histogram

## Boxplot of the estimate of the ACE (degenerate!)
boxplot(hat.ace.vec,main="Distribution of hat.ace",xlab="hat.ace")
abline(h=ace.true,col="red") #True ACE
abline(h=mean(hat.ace.vec),col="blue") #mean of histogram

## Histogram of the estimate of the variance of the ACE
hist(hat.var.hat.ace.vec,main="Distribution of hat.var.hat.ace",
xlab="hat.var.hat.ace.vec",xlim=c(0,max(hat.var.hat.ace.vec)))
abline(v=var(hat.ace.vec),col="red") #True variability of hat.ace
abline(v=mean(hat.var.hat.ace.vec),col="blue") #mean of histogram

## Boxplot of the estimate of the variance of the ACE
boxplot(hat.var.hat.ace.vec,main="Distribution of hat.var.hat.ace",
xlab="hat.var.hat.ace.vec")
abline(h=var(hat.ace.vec),col="red") #True variability of hat.ace
abline(h=mean(hat.var.hat.ace.vec),col="blue") #mean of histogram

## Scatterplot of the estimate of ACE vs. estimate of the variance of the ACE
plot(hat.ace.vec,hat.var.hat.ace.vec, main="hat(ACE) vs. hat(var(hat(ACE)))",
ylab="hat.var.hat.ace", xlab="hat.ace" )
abline(h=mean(hat.var.hat.ace.vec),col="blue") #mean of estimates of variance
abline(h=var(hat.ace.vec),col="red") #True variability of hat.ace (not on plot!)
abline(v=mean(hat.ace.vec), col="blue") #Mean of hat.ACE
abline(v=ace.true,col="red") #True ACE


