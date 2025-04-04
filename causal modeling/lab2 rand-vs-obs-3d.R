### Comparing the set of distributions compatible with data from
### a randomized study vs. an observational study.

rm(list=ls()) # delete anything currently in memory

require(rgl)  #load 3d plotting package

source("https://www.stat.washington.edu/tsr/s566/labs/y0y1polytopenew-rgl-col.R")

## Set up interactive 3d plot
plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1)) 
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)

### Data from 2x2 table, lecture 1, slide 28:

###### Placebo	Drug
######   X=0		X=1
# Y=0    7/20	4/20
# Y=1	 3/20	6/20

### Viewed as a Randomized experiment (was also included last week)

prop.rec.ctrl <- (3/20)/((7/20)+(3/20))   #P(Y=1 | X=0)
prop.rec.trt <-  (6/20)/((4/20)+(6/20))   #P(Y=1 | X=1)

## Blue plane for treatment arm
vec4.trt <- c(0,1-prop.rec.trt,0,prop.rec.trt)
do.polytope.rgl(vec4.trt,red="blue", blue="red")

## Purple plane for control arm
vec4.ctrl <- c(1-prop.rec.ctrl,0,prop.rec.ctrl,0)
do.polytope.rgl(vec4.ctrl,red="purple", blue="red")

## Add a green line indicating the average causal effect 
## and showing the set of distributions compatible with the
## data from both arms. 

do.itt.line.rgl(c(vec4.ctrl,vec4.trt))


################################

## Reset plot

plot3d(c(1,0,0,1),c(0,1,0,0),c(0,0,1,0), ylab="%HU",xlab="%HE",zlab="%AR",par3d(FOV=1)) 
lines3d( c(1,0,0,1),c(0,1,0,0),c(0,0,1,0),col="red",lty=3)

## shape (aka polytope) describing the set of distributions 
## compatible with the same data viewed as an observational study


# p(y0, x0), p(y0, x1), p(y1, x0), p(y1, x1)

vec4.obs <- c(7/20,4/20,3/20,6/20)
do.polytope.rgl(vec4.obs,red="violet", blue="red")

## Add the set of distributions implied by a randomized experiment 
## (This line is inside the polytope. (?Why?))
do.itt.line.rgl(c(vec4.ctrl,vec4.trt))

##

