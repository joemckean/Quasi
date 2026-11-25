#
#' Computes Tables of the effects computed by the Quasi-ANCOVA analysis.
#'	
#' @param ngp number of treatment groups.
#' @param beta vector of estimated regression coefficients from quasif.R.
#' @param cmat is the bootstrap estimated covariance matrix of beta.
#' @return Table of regression estimates with standard errors.
#' @return Table of estimated effects with standard errors.
#' @return Table of estimated treatment means with standard errors.
effectstab <- function(ngp,beta,cmat){

	p <- length(beta)

	sebeta <- sqrt(diag(cmat))
	trat <- beta/sebeta
	tab1 <- round(cbind(beta,sebeta,trat),digits=4)
	tab1 <- cbind(1:p,tab1)
	colnames(tab1) <- c("Param","Regr Cf","SE","t-ratio")

#       write.table(tab1)

	mueff <- beta[1:ngp]
	ceff <- cmat[1:ngp,1:ngp]
	vecs <- rep(0,ngp)
	hmat <- c(); hind <- c()
	for(j in 1:(ngp-1)){
		if(j==1){
			for(i in 1:(ngp-1)){
				vec <- vecs; vec[i+1] = 1
				hmat <- rbind(hmat,vec)
				hind <- rbind(hind,c(i+1,j))
			} 
		} else {
			for(i in (j + 1):ngp){
				veca <- vecs; veca[j] <- 1
				vecb <- vecs; vecb[i] <- 1
				hmat <- rbind(hmat,vecb - veca)
				hind <- rbind(hind,c(i,j))
			}
		}
	}


        effs <- hmat%*%mueff; tab2 <- cbind(hind,effs)
        seeffs <- sqrt(diag(hmat%*%ceff%*%t(hmat)))
        teffs <- effs/seeffs
	tab2 <- cbind(tab2,round(seeffs,digits=4),round(teffs,digits=4))
	colnames(tab2) <- c("_j","_i","mu_j-mu_i","SE","t-ratio")
#        write.table(tab2)

        v1 <- c(1,rep(0,ngp-1))
	hmat1 <- v1
	for(j in 2:ngp){
		vec <- v1; vec[j] <- 1
		hmat1 <- rbind(hmat1,vec)
	}
	mus <- hmat1%*%mueff
	semus <- sqrt(diag(hmat1%*%ceff%*%t(hmat1)))
	tmus <- mus/semus
	tab3 <- cbind(1:ngp,round(mus,digits=4),round(semus,digits=4),round(tmus,digits=4))
	colnames(tab3) <- c("_i","mu_i","SE","t-ratio")



	list(Tableregrcoef=tab1,Tableeffects=tab2,Tablemeans=tab3)
}
