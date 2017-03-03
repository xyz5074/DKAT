\name{wlin.kernel}
\alias{wlin.kernel}
\title{ Weighted Linear Kernel}
\description{
Calculates the weighted linear kernel matrix for genotypes
}
\usage{ wlin.kernel(X, W.beta) }
\arguments{
  \item{X}{Genotype matrix, each row is a sample and each column is a genetic variant}
  \item{W.beta}{ two-dimensional weights as in the beta density function}
}
\details{
Let W=diag(\eqn{w_1,\ldots,w_p}) be the diagonal matrix containing the weights of the p genetic variants, where \eqn{\sqrt{w_j}=beta(MAF_j,a_1,a_2)}, \eqn{MAF_j} is the minor allele frequency of variant j, and (\eqn{a_1,a_2}) are the weights. Then the weighted linear kernel matrix is calculated as 
\eqn{K=XWWX^T}.
}
\value{
A n by n kernel matrix, where n is the number of subjects.
}
\references{
Wu, M. C. et al. (2011). Rare--variant association testing for sequencing data with the sequence kernel associaiton test. The American Journal of Human Genetics, 89, 82--93.
}



