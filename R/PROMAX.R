

# Promax rotation -- Marcus, 1993

PROMAX <- function (loadings, ppower=4, display=TRUE) {

# the loadings are first rotated to an orthogonal varimax solution, then the
# orthogonality of the factors is relaxed to better fit simple structure 
# (Mulaik 1972 p 300; Rummel 1970 p 419)


if (is.list(loadings) == 'TRUE') {loadings <- loadings$loadings}

if (ncol(loadings) == 1 )  {	
	promaxOutput <- list( structure=loadings )
	return(invisible(promaxOutput))
	if (display == TRUE) cat("\n\nWARNING: There was only one factor. Rotation was not performed.\n")
}

if (ncol(loadings) > 1) {

	loadings <- as.matrix(VARIMAX(loadings,display=FALSE))

	B <- loadings

	Bstar <- B ^ppower 

	Tr <- solve(t(B) %*% B)%*% t(B) %*% Bstar

	Tr <- Tr %*% sqrt(solve(diag(diag(t(Tr)%*%Tr))))  # Normalizes columns of Tr

	Tpp <- solve(Tr)        # Tp' from definition

	Tpp <- solve(sqrt(diag(diag(Tpp%*%t(Tpp)))))%*%Tpp # Normalizes rows of Tp'

	Sr <- B %*% Tr           # Oblique reference structure

	Phip <- Tpp %*% t(Tpp)     # Correlation between primary factors  # Table 6.VI

	Sp <- B %*% t(Tpp)    # Primary Factor structure matrix Table 6.V

	Pp <- B %*% solve(Tpp) # Primary Factor pattern matrix Table 6.IV

	Pr <- B %*% solve(t(Tr))  # Reference pattern matrix

	Phir <- t(Tr) %*% Tr    # Correlations between reference axes

	colnames(Sp)   <-  c(paste("  Factor ", 1:ncol(Sp),   sep="") )
	colnames(Pp)   <-  c(paste("  Factor ", 1:ncol(Pp),   sep="") )
	colnames(Phip) <-  c(paste("  Factor ", 1:ncol(Phip), sep="") )
	rownames(Phip) <-  c(paste("  Factor ", 1:ncol(Phip), sep="") )


	if (display == TRUE) {
		# cat("\n\nUnrotated Loadings:\n\n")
		# print(round(B,2))
		cat("\n\nPromax Rotation Structure Matrix:\n\n")
		print(round(Sp,2))
		cat("\n\nPromax Rotation Pattern Matrix:\n\n")
		print(round(Pp,2));cat("\n")
		cat("\n\nPromax Rotation Factor Correlations:\n\n")
		print(round(Phip,2));cat("\n\n")
	}

	promaxOutput <- list( structure=Sp, pattern=Pp, correls=Phip  )

	return(invisible(promaxOutput))
}

}
