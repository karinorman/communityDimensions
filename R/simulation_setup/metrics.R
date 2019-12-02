
#### kl divergence ####
klDivergence <- function(cov1, cov2){
  ## natural log
  wantT <- solve(cov2) %*% cov1
  
  log(det(cov2)/det(cov1))/2 - nrow(cov1)/2 + sum(diag(wantT))/2

}
## this isn't going to work for low rank matrices, not invertible

helperCreateMatrixNice <- function(sparsity, signal, d) {
  precision.mat <- diag(d)
  precision.mat[upper.tri(precision.mat)] <- signal * (runif(d * (d - 1) / 2) < sparsity)
  
  precision.mat <- precision.mat + t(precision.mat)
  
  ## Check Pos-Def
  eig <- eigen(precision.mat)
  
  if (sum(eig$values > 0) < d) {
    precision.mat <- nearPD(precision.mat) ## adjust if not pos-def
    # eig <- eigen(precision.mat$mat) # can print to see what's happening
    return(precision.mat$mat)
  } else {
    # eig <- eigen(precision.mat)
    return(precision.mat)
  }
}

klDivergence(helperCreateMatrixNice(0.5, 0.5, 5), helperCreateMatrixNice(0.5, 0.5, 5))

## common sense checks 
test = helperCreateMatrixNice(0.5, 0.5, 5)
test2 = test
test2[1,1]=5
test3 = test2
test3[1,1]=3
klDivergence(test, test)
klDivergence(test, test2)
klDivergence(test, test3)


#### frobenius norm ####

## probably not technically the frobenius norm, but inspired by it

frobeniusNorm <- function(estimatedMat, truthMat){
  sqrt(sum(abs(truthMat - estimatedMat)))
}

frobeniusNorm(matrix(rnorm(25), 5, 5),  matrix(rnorm(25), 5, 5))

test = matrix(rnorm(25), 5, 5)
test2 = test
test2[1,1]=-100
frobeniusNorm(test, test)
frobeniusNorm(test, test2)
