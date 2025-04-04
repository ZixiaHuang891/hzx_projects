# Clear everything
rm(list=ls())

### Load required packages
##  Separate commands because some are large and the
##  download command may be slow or fail.
##  (However, if you re-run the commands it will not reload
##   packages that are already there.)

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")}
    
if (!requireNamespace("Rgraphviz", quietly = TRUE))
    BiocManager::install("Rgraphviz",update=FALSE)
library("Rgraphviz")

if (!requireNamespace("RBGL", quietly = TRUE)){
    BiocManager::install("RBGL",update=FALSE)}
library("RBGL")

if (!requireNamespace("abind", quietly = TRUE)){
    install.packages("abind",update=FALSE)}
library("abind")

if (!requireNamespace("corpcor", quietly = TRUE)){
    install.packages("corpcor",update=FALSE)}
library("corpcor")

if (!requireNamespace("sfsmisc", quietly = TRUE)){
    install.packages("sfsmisc",update=FALSE)}
library("sfsmisc")

if (!requireNamespace("robustbase", quietly = TRUE)){
    install.packages("robustbase",update=FALSE)}    
library("robustbase")

if (!requireNamespace("pcalg", quietly = TRUE)){
    install.packages("pcalg",update=FALSE)}
library("pcalg")

if (!requireNamespace("graph", quietly = TRUE)){
    install.packages("graph",update=FALSE)}
library("graph")

## Boolean checking whether Rgraphviz loaded
## (So people can continue in case not!)

plotcpdag <- "Rgraphviz" %in% print(.packages(lib.loc = .libPaths()[1])) 
  
##################################################
## PC Algorithm example, using discrete data
##################################################
## Load data (provided with package)
data(gmD)
V <- colnames(gmD$x) # gmD$x gives a set of 5 columns

# Look at the data
gmD$x[1:10,]
str(gmD$x)

## The number of unique values of each variable:
sapply(gmD$x, function(v) nlevels(as.factor(v)))

# X1 X2 X3 X4 X5 
#  3  2  3  4  2 

## define sufficient statistics
suffStat <- list(dm = gmD$x, nlev = c(3,2,3,4,2), adaptDF = FALSE)

## estimate CPDAG
pc.D <- pc(suffStat,
           ## independence test: G^2 statistic
           indepTest = disCItest, alpha = 0.01, labels = V, verbose = TRUE)
if (require(Rgraphviz)) {
  ## show estimated CPDAG
  par(mfrow = c(1,2))
  plot(pc.D, main = "Estimated CPDAG",labels=c("X1","X2","X3","X4","X5"))
#######
# Plot true graph for comparison	 
	 true.graph <- gmD$g
	# change labels
	tmp <- c("X1","X2","X3","X4","X5")
	names(tmp) <- nodes(gmD$g)
	nAttrs <- list()
	nAttrs$label <- tmp
	rm(tmp)
	# plot with node labels
	plot(true.graph, nodeAttrs=nAttrs, attrs=list(node=list(shape="ellipse",fontsize=8)), main = "True DAG")
}

# Note that the X4 -> X5 edge has the wrong orientation in the PC Output.
# This is likely due to a faithfulness violation.



