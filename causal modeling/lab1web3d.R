###############
# Class lab 1 #
###############

rm(list=ls())  # Clear objects

set.seed(as.numeric(Sys.Date()))  # Set the random seed

# 0. SET UP
# load our set of subjects
cntrfact <- read.table("https://www.stat.washington.edu/tsr/s566/labs/cntrfact1.dat")

# see the types of first 10 patients
cntrfact[1:10,]

# make vectors for the counterfactual outcomes
ydrug <- cntrfact[,2]
yplac <- cntrfact[,1]
## the researcher typically will not know these

#find out how many subjects there are:
n <- length(yplac)

# 1. RANDOMLY ASSIGN TREATMENT
# let's assign some treatments

pdrug <- 0.5 #probability to receive drug
tmp <- rbinom(n,1,pdrug) 
trt <- sapply(tmp+1,switch,"Plac","Drug"); rm(tmp)

#2. GET THE OBSERVED OUTCOMES
#let's get the observed outcomes

yobs <- ydrug*(trt == "Drug") + yplac*(trt == "Plac")
# try to work out what this did 
# Hint: "TRUE" corresponds to 1; "False" to 0

## Let's see what the researcher will get to see
## after an experiment is performed:
trt[1:10]
yobs[1:10]

hist(yobs)
table(yobs)

table(trt)

table(trt,yobs)

#3. FIND PROP. SUCCESS IN THE DRUG GROUP
# let's compute the proportion of successes in the drug group
mean(yobs[trt=="Drug"])

# compare this to the prop of successes if everyone got the drug
mean(ydrug)  ## this is not something the researcher can see

#4. FIND PROP. SUCCESS IN THE CONTROL GROUP
# and in the control group
mean(yobs[trt=="Plac"])

# compare this to the prop of successes if everyone got placebo
mean(yplac)  ## again, not something the researcher can see

#5. COMPUTE OUR ESTIMATE OF THE AVERAGE CAUSAL EFFECT:
# hence our estimate of the average causal effect is:
mean(yobs[trt=="Drug"])  - mean(yobs[trt=="Plac"])

# compare this to the true mean average causal effect:
mean(ydrug) - mean(yplac)
# any difference between these two is due to sampling variability
# i.e. it would disappear as n goes to infinity

#6. LOOK AT THE TYPES OF PATIENT IN OUR STUDY:
# Note that we could have done this before we assigned treatment
# again, this is not something the researcher can do

helped <- sum(ydrug*(1-yplac))
hurt <- sum((1-ydrug)*yplac)
alwaysrec <- sum(ydrug*yplac)
doomed <- sum((1-ydrug)*(1-yplac))

# how about a nice table:

noquote( cbind(c("helped","hurt","alwaysrec","doomed"),c(helped,hurt,alwaysrec,doomed),c(helped,hurt,alwaysrec,doomed)/n) )

#alternative display
table(ydrug,yplac)

table(ydrug,yplac)/n  #compare to table on slide 19

####### 3d plot

require(rgl)  #3d plotting package

source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")

## Set up interactive 3d plot
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1)) 
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)

## Add a blue plane for data from treatment arm:

prop.rec.trt <- mean(yobs[trt=="Drug"])
vec4.trt <- c(0,1-prop.rec.trt,0,prop.rec.trt)
do.polytope.rgl(vec4.trt,red="blue", blue="red")

## Note: don't worry about the 4-vector that we are inputting here;
##       this will be explained later in the course.
##       (The function is originally designed to handle observational studies.)

## Add a purple plane for data from control arm:

prop.rec.ctrl <- mean(yobs[trt=="Plac"])
vec4.ctrl <- c(1-prop.rec.ctrl,0,prop.rec.ctrl,0)
do.polytope.rgl(vec4.ctrl,red="purple", blue="red")

## Add a green line indicating the average causal effect 
## and showing the set of distributions compatible with the
## data from both arms. 

do.itt.line.rgl(c(vec4.ctrl,vec4.trt))

## Note: Here we are ignoring sampling variability by treating the
## observed proportion recovering as if they were the true population values.


#### (Extra)
####
#### Non-interactive version 
#### (Useful for including plots in documents)

par(mfrow=c(1,2))
simp <- do.simplex(phi=30,theta=120,r=1000,main="Experimental study") # Set up plot
do.polytope(vec4.trt,simp,red="blue", blue="red")
do.polytope(vec4.ctrl,simp,red="purple", blue="red")
do.itt.line(c(vec4.ctrl,vec4.trt),simp)

simp <- do.simplex(phi=-90,theta=90,r=500,main="Exp Study from below") # Set up plot
do.polytope(vec4.trt,simp,red="blue", blue="red")
do.polytope(vec4.ctrl,simp,red="purple", blue="red")
do.itt.line(c(vec4.ctrl,vec4.trt),simp)

#####################################################
## How to close RGL so that you can the quit XQuartz
##
## rgl.quit()
##
#####################################################

#################################################
## Further exercises:

# (I)
# Now repeat from the step 1.
# When you repeat the subsequent steps, what is the same?
# What is (slightly) different?

# (II)
# Now repeat from step 1. but this time change pdrug
# to a value (strictly) between 0 and 1:
pdrug <-  #your favourite probability here#

# What has changed?
# Is our estimate of the average causal effect very different?


# (III)
# Now repeat from step 0. with the data set cntrfact2.dat:
cntrfact <- read.table("https://www.stat.washington.edu/tsr/s566/labs/cntrfact2.dat")

#repeat all steps
# pay particular attention to:
# 5. the true and estimated ACE 
# 6. the matrix of patient types


# (IV)
# Now repeat from step 0. with the data set cntrfact3.dat:
cntrfact <- read.table("https://www.stat.washington.edu/tsr/s566/labs/cntrfact3.dat")

#repeat all steps
# pay particular attention to:
# 5. the true and estimated ACE 
# 6. the matrix of patient types
# contrast these to those obtained from (III)

# (V)
# Looking at your answers to (II) and (III) is it correct to say
# "If the average causal effect is zero, then the drug 
#   is doing nothing"? 
#  Why or why not?
# Hint: think about the perspective of an individual (or his lawyer)
#       vs. the perspective of a policy analyst.
