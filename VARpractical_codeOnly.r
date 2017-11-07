#####################################################################
#### 1.1 Libraries
#####################################################################

# load libraries
library(Biobase)
library(lattice)
library(longitudinal)
library(rags2ridges)
library(ragt2ridges)
library(SparseTSCGM)


#####################################################################
#### 1.2 Data
#####################################################################

# load data
data(hpvP53)

# reformat data
Y <- longitudinal2array(t(exprs(hpvP53rna)))

# zero center data, variate- and cell line-wise
Y <- centerVAR1data(Y)


#####################################################################
#### 1.3 Data exploration
#####################################################################

# plot time-courses of all variates
plotVAR1data(Y)

# plot time-courses of a single variate
plotVAR1data(Y[5, , , drop=FALSE])

# K-means clustering of the variates
cellLine <- 1
kClust   <- kmeans(Y[, , cellLine], centers=7, nstart=100)$cluster

# heatmap of reshuffled data
edgeHeat(Y[unlist(lapply(1:max(kClust), 
                         function(id, clusters){ which(clusters==id) }, 
                         kClust)), , cellLine])

# plot time-courses of smallest cluster
plotVAR1data(Y[kClust==which.min(table(kClust)), , , drop=FALSE])


#####################################################################
#### 1.4 Penalty selection
#####################################################################

# find optimal penalty parameters
optLambdas <- optPenaltyVAR1(Y, lambdaMin =c(0.01, 0.00001), 
                                lambdaMax =c(1000, 1), 
                                lambdaInit=c(100 , 0.1))

# print optimal penalty parameters
print(optLambdas)

# define grid for both penalty parameters
lambdaAgrid <- seq(400,    450, length.out=20)
lambdaPgrid <- seq(0.001, 0.01, length.out = 20)

# evaluate LOOCV log-likelihood over the grid
LOOCVres    <- loglikLOOCVcontourVAR1(lambdaAgrid, lambdaPgrid, 
                                      Y, verbose=FALSE)

# plot LOOCV log-likelihood somewhat nicer
contour(lambdaAgrid, lambdaPgrid, LOOCVres$llLOOCV, xlab="lambdaA", 
        ylab="lambdaP", main="cross-validated log-likelihood", nlevels=25)

# add optimal penalty parameters as red dot
points(optLambdas[1], optLambdas[2], pch=20, cex=2, col="red")


#####################################################################
#### 1.5 Penalized estimation
#####################################################################

# fit the model
VAR1hat        <- ridgeVAR1(Y=Y, 
                            lambdaA=optLambdas[1], 
                            lambdaP=optLambdas[2])

# extract parameter estimates
Ahat           <- VAR1hat$A
Phat           <- VAR1hat$P

# add row and column names
rownames(Ahat) <- colnames(Ahat) <- rownames(Phat) <-
                  colnames(Phat) <- rownames(hpvP53rna)

# heatmap of estimate of A
edgeHeat(Ahat, main="ridge estimate of A")

# heatmap of partial correlation matrix estimate
edgeHeat(pcor(Phat), main="ridge partial correlation estimate", diag=FALSE)


#####################################################################
#### 1.6 Stability analysis
#####################################################################

# leave each time point out once and re-estimate
Aperturb <- numeric()
for (u in 1:((dim(Y)[2]-1)*dim(Y)[3])){
  # leave one sample out
	unbalanced <- cbind(rep(2:dim(Y)[2], dim(Y)[3]), 
                      sort(rep(1:dim(Y)[3], dim(Y)[2]-1)))[u, , drop=FALSE]

	# choose penalty parameters and estimate model
	newLambda <- optPenaltyVAR1(Y, 
                              lambdaMin =c(100, 0.001),
                              lambdaMax =c(1000, 1), 
                              lambdaInit=optLambdas, 
                              optimizer ="optim", 
                              unbalanced=unbalanced)

	# refit
	Aperturb  <- cbind(Aperturb, 
                     as.numeric(ridgeVAR1(Y, 
                                          lambdaA   =newLambda[1], 
                                          lambdaP   =newLambda[2], 
                                          unbalanced=unbalanced)$A))
}

# select hundred equidistant autoregression parameters
ids <- match(ceiling(seq(1, nrow(Aperturb), length.out=100)), 
             rank(Ahat))

# plot boxplots of the left-out estimates
boxplot(as.numeric(Aperturb[ids,]) ~ rep(1:length(ids), ncol(Aperturb)), 
        col ="blue", border="lightblue", pch=20, cex=0.5, 
        xlab="original estimate ranking, proportionally", 
        ylab="entry of perturbed A estimate")


#####################################################################
#### 1.7 Support determination
#####################################################################

# support determination of A
zerosA <- sparsifyVAR1(A=Ahat, SigmaE=symm(solve(Phat)), 
                       threshold="top", top=25, 
                       statistics=FALSE, verbose=FALSE)$zeros

# support determination of precision matrix
zerosP <- sparsify(Phat, threshold="top", top=10, 
                   output="light", verbose=FALSE)$zeros


#####################################################################
#### 1.8 Re-estimation with inferred support
#####################################################################

# format precision support
supportP <- support4ridgeP(zeros=zerosP, nNodes=nrow(Y))

# optimal penalty parameter determination
optLambdas <- optPenaltyVAR1(Y, lambdaMin=c(10^(-5), 10^(-5)), 
                             lambdaMax=c(10, 0.1), lambdaInit=c(5, 0.01), 
                             zerosA=zerosA, zerosP=zerosP, 
                             cliquesP=supportP$cliques, 
                             separatorsP=supportP$separators,                              
                             zerosAfit="sparse")
print(optLambdas)

# define grid for both penalty parameters
lambdaAgrid <- seq(0.1, 5, length.out=20)
lambdaPgrid <- seq(0.001, 0.02, length.out = 20)

# evaluate LOOCV log-likelihood over the grid
LOOCVres    <- loglikLOOCVcontourVAR1(lambdaAgrid, lambdaPgrid, 
                                      Y, zerosA=zerosA, zerosP=zerosP, 
                                      cliquesP=supportP$cliques, 
                                      separatorsP=supportP$separators,
                                      zerosAfit="sparse")

# plot LOOCV log-likelihood somewhat nicer
contour(lambdaAgrid, lambdaPgrid, LOOCVres$llLOOCV, 
        xlab="lambdaA", ylab="lambdaP", main="cross-validated 
        log-likelihood", nlevels=25)

# add optimal penalty parameters as red dot 
points(optLambdas[1], optLambdas[2], pch=20, cex=2, col="red")

# re-fit the model
VAR1hat        <- ridgeVAR1(Y=Y, lambdaA=optLambdas[1], 
                            lambdaP=optLambdas[2], zerosA=zerosA,
                            cliquesP=supportP$cliques, 
                            separatorsP=supportP$separators,
                            zerosP=zerosP, zerosAfit="sparse")

# extract parameter estimates
Ahat           <- VAR1hat$A
Phat           <- VAR1hat$P

# add row and column names
rownames(Ahat) <- colnames(Ahat) <- rownames(Phat) <- 
                  colnames(Phat) <- rownames(hpvP53rna)

# heatmap of support of A
edgeHeat(adjacentMat(Ahat), legend=FALSE, main="inferred support of A")

# heatmap of support of precision matrix
edgeHeat(adjacentMat(Phat), legend=FALSE, main="inferred support of
                            precision matrix", diag=FALSE)

# heatmap of estimate of A
edgeHeat(Ahat, main="ridge re-estimate of A with inferred support")

# heatmap of support of the partial correlation matrix
edgeHeat(pcor(Phat), main="ridge re-estimate of the partial correlation matrix
         with inferred support", diag=FALSE)

# time-series chain graph
graphVAR1(Ahat, Phat, nNames=rownames(Ahat), type="TSCG", 
          vertex.label.cex=0.5, vertex.label.font=1, vertex.size=4, 
          vertex.label.color.T0="darkblue", 
          vertex.label.color.T1="darkblue", 
          vertex.frame.color="steelblue", vertex.color.T0="lightblue", 
          vertex.color.T1="lightblue", edge.width=1.5, main = "")

# determine adjacency matrix
adjMatCIG <- CIGofVAR1(Ahat, Phat)

# plot global conditional independence graph
graphVAR1(Ahat, Phat, nNames=rownames(Ahat), type="globalPC")


#####################################################################
#### 1.9 Fit assessment
#####################################################################

# calculate fits and fit-vs-observation correlations
Yhat <- array(dim=dim(Y))
for (i in 1:dim(Y)[3]){ Yhat[, -1, i] <- Ahat %*% Y[, -dim(Y)[2], i] }

# specify molecular entity of interest
entityID <- 1

# format data for plotting
label                <- paste("CellLine", sort(rep(c(1:dim(Y)[3]), 
                                                   dim(Y)[2]-1)))
time                 <- rep(2:dim(Y)[2], dim(Y)[3])
entityData           <- data.frame(as.numeric(Y[entityID, -1, ]), 
                                   time, label)
colnames(entityData) <- c("Data", "TimePoints", "CellLine")

# plot data and fit for each cell line separately 
xyplot(Data ~ TimePoints | CellLine, data=entityData, 
       layout=c(dim(Y)[3], 1), col="blue", pch=20, aspect=1)
for (i in 1:dim(Y)[3]){
    trellis.focus("panel", i, 1)
    llines(Yhat[entityID, -1, i] ~ c(2:dim(Y)[2]), col="red", lwd=2)
}
trellis.unfocus()

# calculate fits and fit-vs-observation correlations
for (i in 1:dim(Y)[3]){ Yhat[, -1, i] <- Ahat %*% Y[, -dim(Y)[2], i] }
corFit <- numeric()
for (j in 1:dim(Y)[1]){
    slHelper <- numeric()
    for (i in 1:dim(Y)[3]){ 
      slHelper <- c(slHelper, cor(Yhat[j, -1, i], Y[j, -1, i], m="s")) 
    }
    corFit <- rbind(corFit, slHelper)
}

# histogram of the fit-vs-observation correlations
hist(corFit, xlab="Correlation", ylab="Frequency", main="Histogram of the
     correlation fit vs. observation", n=20, col="blue", border="lightblue", 
     xlim = c(-1, 1))

# "previous observation"-forecasts
for (i in 1:dim(Y)[3]){ Yhat[, -1, i] <- Y[, -dim(Y)[2], i] }

# "other cell lines"-forecasts
for (i1 in 1:dim(Y)[3]){
  Yslh <- 0 * Y[, , 1]
  for (i2 in 1:dim(Y)[3]){
    if (i2 != i1){ Yslh <- Yslh + Y[, , i2] }
  }
	Yhat[, , i1] <- Yslh / (dim(Y)[3]-1)
}


#####################################################################
#### 1.10 Down-stream analysis
#####################################################################

# calculate node-wise network stats 
nodeStats <- nodeStatsVAR1(Ahat, Phat, as.table=TRUE)

# show node degree table
print(nodeStats[, 1:7])

# specify time lag
lag <- 1

# evaluate mutual informations with specified lag
MIs <- mutualInfoVAR1(Ahat, solve(Phat), lag)

# specify time lag
lag <- 1

# evaluate impulse response with specified lag
IRs <- impulseResponseVAR1(Ahat, lag)

#####################################################################
#### 1.11 Penalized estimation -- lasso
#####################################################################

# reformat data
Y <- as.longitudinal(array2longitudinal(Y[1:10,,]))

# find BIC-optimal penalty parameters and fit VAR(1) model
VAR1lasso <- sparse.tscgm(Y, optimality = "bic")

# extract parameters
AhatL           <- VAR1lasso$gamma
PhatL           <- VAR1lasso$theta
rownames(AhatL) <- colnames(AhatL) <- rownames(PhatL) <- 
                   colnames(PhatL) <- rownames(hpvP53rna)[1:10]

#####################################################################
#### 1.12 Simulate data
#####################################################################

# set dimensions (p=covariates, n=individuals, T=time points)
p <- 3; n <- 4; T <- 10

# set model parameters
A      <- -createA(p, "clique", nCliques=1, nonzeroA=0.1)
SigmaE <- diag(p)/4

# generate data
Y <- dataVAR1(n, T, A, SigmaE)

