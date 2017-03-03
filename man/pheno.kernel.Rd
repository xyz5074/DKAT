\name{pheno.kernel}
\alias{pheno.kernel}
\title{ Phenotype Kernel}
\description{
Calculates the kernel matrix for multivariate (potentially high-dimensional and structured) phenotypes 
}
\usage{pheno.kernel(Y, rho = 0.1)}

\arguments{
  \item{Y}{Phenotype matrix, each row is a sample and each column is a phenotype}
  \item{rho}{Graphical lasso regularization parameter used in estimating the precision matrix of phenotypes}
}
\details{ Let \eqn{\Theta} be the graphical lasso estimator of the precision matrix of phenotypes. Then the phenotype kernel matrix is calculated as  \eqn{K=Y \Theta Y^T}.

}
\value{ 
	A n by n kernel matrix, where n is the number of subjects.
}
\references{
Friedman, J. et al. (2008). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 9, 432--441. \cr
Zhan, X. et al. (2017). Powerful genetic association analysis for common or rare varaints with high--dimensional structured tratis. Genetics, submitted.
}

