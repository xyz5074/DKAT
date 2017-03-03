\name{DKAT}
\alias{DKAT}
\title{ Dual Kernel-based Association Test}
\description{
Test the association between genotypes and phenotypes using two kernel matrices.
}
\usage{ DKAT(K, L) }

\arguments{
  \item{K}{ A n by n kernel matrix calculated from genotypes}
  \item{L}{ A n by n kernel matrix calculated from phenotypes}
}
\details{The DKAT test statistic is the normalized frobenius inner product between two centered kernel matrices, and permutations are used to evaluate the significance based on this DKAT statistic. Instead of explicitly 
drawing permutations and calculating the permuted DKAT statistics, we approxiamte the permutation null distribution to a Pearson type III distribution. Details of this approximation can be found in Zhan et al., (2017).
}
\value{
A p-value
}
\references{
Zhan, X. et al. (2017). Powerful genetic association analysis for common or rare varaints with high--dimensional structured tratis. Genetics, submitted.
}

\examples{
X = matrix(rbinom(200,2,0.3), nrow=200,ncol=1) 	## 200*1 SNP vector
Y = matrix(rnorm(1000,0,1),nrow=200,ncol=5)	## 200*5 normal outcomes
W.beta=c(1,1) ## Beta-weights: flat for SNPs and beta(1,25) for RV
KX=wlin.kernel(X,W.beta)
KY=pheno.kernel(Y)
pv=DKAT(KX,KY)

}
