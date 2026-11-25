#
#' Computes Quasi-ANCOVA F-test statistic and its bootstrap p-value.
#'
#' \code{quasidriver} function computes the Quasi-ANCOVA analysis of one-way treatment design with
#' covariates, some of which may have been affected by the treatments and will
#' be transformed to the space orthogonal to the oneway designi space.  Besides 
#' the test of effects tables of estimates of contrasts of treatment effects and regression coefficients
#' based on quasi-ANCOVA fit.
#' 
#' \code{quasibs3} also computes the  Quasi-ANCOVA analysis.   The difference is that
#' \code{quasidriver} formats the output of \code{quasibs3}.
#' 
#' @param mat matrix with columns: treatment indicator (integers), responses, covariates.
#' @param vecss vector of treatment sample sizes.
#' @param work with jth component 1 if jth covariate is to be transformed, else 0.
#' @param nbs is the number of bootstraps, default is 3000.
#' @param tabs logical, if TRUE prints the table of effects with standard errors.
#' @param intseed integer, Initial seed for R's random number generator.
#' @return Fstat Quasi-ANCOVA test statistic.
#' @return Bootpval the bootstrap p-value of the Quasi-ANCOVA test of treatment effect.
#' @return Regest the vector of estimated regression coefficients.
#' @return Tableregrcoef table of regression coefficients and their bootstrapped standard errors.
#' @return Tableeffects table of estimated effects with bootstrapped standard errors.
#' @return Tablemeans table of estimated treatment means with bootstrapped standard errors.
#' @return covmat bootstrapped estimate of the covariance matrix of the regressioin coefficients.
#
#

quasidriver <- function(mat,work,nbs=10000,intseed=832659){
	print("Argument #1 is a matrix",quote=FALSE)
	print("      First column is the treatment indicator, components must be positive integers",quote=FALSE)
	print("      Second column contains the responses",quote=FALSE)
	print("      The remaining columns contain the covariates",quote=FALSE)
	print("",quote=FALSE)
	print("Argument #2 is the work vector with length the number of covariates.",quote=FALSE)
	print("      If the jth component is 1 the jth covariate is transformed to the null space",quote=FALSE)
	print("      of the design (treatment) space (covariate may have been affected by the treatments)",quote=FALSE)
	print("      If the jth component is 0 the jth covariate is not transformed.",quote=FALSE)
	print("Argument #3 is the number of bootstrap iterations; default is 10,000",quote=FALSE)
	print("",quote=FALSE)
	print("Argument #4 is the seed that initializes R's random generator, it must be a positive integer.",quote=FALSE)

	set.seed(intseed)
	vecss <- table(mat[,1])
	qac <- quasibs3(mat,vecss,work,nbs,tabs=TRUE,intseed)
	print("",quote=FALSE)
	print("",quote=FALSE)
	print(" Quasi-ANCOVA Test Statistic and p-value for testing no treatment effects",quote=FALSE)
	stuff1 <- c(qac$Fstat,qac$Bootpval); print(stuff1)
	print("",quote=FALSE)
	print(" Table of Estimated Effects with Standard Errors",quote=FALSE)
	print(qac$Tableeffects)
	print("",quote=FALSE)
	print(" Table of Estimated Regression Coefficients with Standard Errors",quote=FALSE)
	print(qac$Tableregrcoef)
	print("",quote=FALSE)
	print(" Table of Estimated Means Based on QuasiANCOVA Fit",quote=FALSE)
	print(qac$Tablemeans)
	print("",quote=FALSE)
	print(" Estimated Covariance Matrix Based on QuasiANCOVA Fit",quote=FALSE)
	print(round(qac$covmat,digits=4))
        return()
}


