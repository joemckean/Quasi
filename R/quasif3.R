#
#' Computes Quasi-ANCOVA Test Statistic.
#'
#' @param mat matrix with columns: treatment indicator (integers), responses, covariates.
#' @param vecss vector of treatment sample sizes.
#' @param work with jth component 1 if jth covariate is to be transformed, else 0.
#' @return fitf the LS fit of the Quasi-ANCOVA model.
#' @return des the design matrix of the Quasi-ANCOVA model.
#' @return qF the Quasi-ANCOVA F-statistic.
#' @return mat2 the input matrix mat reordered by treatment.
#' @return design matrix for treatments: references the first treatment.
#' @return work the input work vector.
#
#
quasif3 <-
function(mat,vecss,work){
#
#   n is the total sample size
	n <- sum(vecss)
#   ngp is the number of treatments
        ngp <- length(vecss)
	ord <- order(mat[,1])
	mat2 <- mat[ord,]
	ind <- mat2[,1]
	kp2 <- length(mat[1,])
#   k is the number of covariates
	k <- kp2 - 2
#   y is the response vector
	y <- mat2[,2]
	xd <- model.matrix(~factor(ind))
#   xd is the design matrix for the effects
        numdf <- length(xd[1,]) - 1
	xd <- xd[,-1]
	xtr <- c()
	xc <- c()
#   transformations of the indicted, by work, covariates
	for(j in 1:k){
		if(work[j]==1){
			fit0 <- lm(mat2[,j+2]~xd)
			xtr <- cbind(xtr,fit0$resid)
		} else {
			xc <- cbind(xc,mat2[,j+2])
		}
	}
	des <- cbind(xd,xc,xtr)
	dendf <- n - length(des[1,]) - 1
#   ANCOVA fit.
	fitF <- lm(y~des)
#   Reduced design matrix for H_0: no treatment effect.
	desred <- cbind(xc,xtr)
	fitr <- lm(y~desred)
	qF <- ((sum(fitr$resid^2) - sum(fitF$resid^2))/numdf)/(sum(fitF$resid^2)/dendf)
#       qpv <- 1 - pf(qF,numdf,dendf)
	list(fitF=fitF,des=des,qF=qF,mat2=mat2,xd=xd,work=work)
}
