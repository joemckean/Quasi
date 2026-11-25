#
#' Computes Quasi-ANCOVA F-test statistic and its bootstrap p-value.
#'
#' \code{quasibs3} function computes the Quasi-ANCOVA analysis of one-way treatment design with
#' covariates, some of which may have been affected by the treatments and will
#' be transformed to the space orthogonal to the oneway design.  Besides 
#' the test of effects tables of estimates of contrasts of treatment effects and regression coefficients
#' based on quasi-ANCOVA fit.
#'
#' \code{quasidriver} also computes the  Quasi-ANCOVA analysis.   The difference is that
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
quasibs3 <-
function(mat,vecss,work,nbs=3000,tabs=TRUE,intseed){
#
#       mat contains the columns: group indicator (ind),
#       response (y), and the last columns are the covariates.
#       vecss is the vector of group sample sizes.
#       Work is a vector with length number of covariates.
#       If work[i] = 1, the ith covariate will be
#       transformed to the space that is orthogonal to effect space,
#       otherwise work[i] is 0.
#
#
#   Fit to obtain the quasi fit on the given data 
	set.seed(intseed)
	fit0 <- quasif3(mat,vecss,work)
	betaF <- fit0$fitF$coef
	sigF <- var(fit0$fitF$resid)
#   The quasi F-test statistic of no treatment effect
	qF <- fit0$qF
	qpv <- fit0$qpv
#   Number of groups
	ngp <- length(vecss)
#   Total sample size
	n <- sum(vecss)
#   mat2 is the original mat sorted by the treatment indicator.
	mat2 <- fit0$mat2
#   Row number for the resampling.
	int <- 1:n
#   Treatment indicator
	ind <- mat2[,1]
	kp2 <- length(mat2[1,])

	indbs <- 0
	pF <- length(betaF)
	covmat <- matrix(rep(0,pF^2),ncol=pF)
	for(i in 1:nbs){
		collint <- c()
                for(j in 1:ngp){
			inttemp <- sample(int,vecss[j],replace=TRUE)
			collint <- c(collint,inttemp)
		}
#   mats contains the resampled (bootstrapped) mat2
		mats <- mat2[collint,2:kp2]
		mats <- cbind(ind,mats)
#   quasi fit of the bootstraped data.
                fits <- quasif3(mats,vecss,work)
		sigs <- var(fits$fitF$resid)
		betas <- fits$fitF$coef
#   Update of bootstrap estimate of the var-cov matrix of the quasi betahat.
		covmat <- covmat + (betas-betaF)%*%t(betas-betaF)/sigs
#   The bootstrapped quasi F-test statistic of no treatment effect.
		qfs <- fits$qF
		if(qfs > qF){indbs <- indbs + 1}
	}

#   Bootstrap estimate of the var-cov matrix of the quasi betahat.
	covmat <- (sigF/nbs)*covmat
#   The bootstrapped p-value for the quasi F-test of no treatment effect.
	pvbs <- indbs/nbs
	if(!tabs){
	         list(Fstat=qF,Bootpval=pvbs,Regest=betaF,covmat=covmat,qpv=qpv)
	} else {
		t123 <- effectstab(ngp,betaF,covmat)
		list(Fstat=qF,Bootpval=pvbs,Regest=betaF,Tableregrcoef=t123$Tableregrcoef,Tableeffects=t123$Tableeffects,Tablemeans=t123$Tablemeans,covmat=covmat)

	}
}
