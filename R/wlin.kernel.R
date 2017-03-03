wlin.kernel <-
function(X, W.beta) {
MAF=colMeans(X)/2  ## colMeans is in package MSKAT
a1=W.beta[1]
a2=W.beta[2]
w=dbeta(MAF,a1,a2)
W=diag(w)
return( X%*%W%*%W%*%t(X) )  
}
