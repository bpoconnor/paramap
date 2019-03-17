
# parallel analysis eigenvalues (no real data required)

PARALLEL <- function ( Nvars=50, Ncases=300, Ndatasets=100, extract='PCA', percentile=95,
                       corkind='pearson', display=TRUE){

evals <- matrix(0, nrow = Nvars, ncol = Ndatasets)
pb <- utils::txtProgressBar(min = 0, max = Ndatasets, style = 3) 

for (nds in 1:Ndatasets) { 

	utils::setTxtProgressBar(pb,nds)

#	cat("Progress - Random Dataset: ", nds, "of ", Ndatasets, "\r"); utils::flush.console()

	randat <- matrix(rnorm(Ncases*Nvars),nrow=Ncases,ncol=Nvars)

	# random data correlation matrix
	if (corkind=='pearson')      Rrand <- cor(randat, method="pearson") 
	if (corkind=='kendall')      Rrand <- cor(randat, method="kendall") 
	if (corkind=='spearman')     Rrand <- cor(randat, method="spearman") 
		
	# random data eigenvalues
	if (extract=='PCA')  evals[,nds ] <- eigen(Rrand) $values 

	if (extract=='PAF') {
		smc = 1 - (1 / diag(solve(Rrand)))
		diag(Rrand) <- smc
		evals[,nds] <- eigen(Rrand) $values 
	}

	if (extract=='image') {
		d <-  diag( 1 / diag(solve(Rrand)) )
		gvv <- Rrand + d %*% solve(Rrand) %*% d - 2 * d
		s <- sqrt(d)                  #  Velicer 1974 p 565 formula (7)
		r2 <- solve(s) %*%  gvv  %*% solve(s)  # Velicer 1974 p 565 formula (5)
		evals[,nds] <- eigen(r2) $values
		} 	
} 


# mean & percentile eigenvalues for each position
means <- apply(evals, 1, mean) 
# sorting the eigenvalues for each root
for (luper in 1:Nvars) { evals[luper,] <- sort(evals[luper,]) }
percentiles <- as.matrix(evals[,round((percentile*Ndatasets)/100)])


results <- cbind(1:Nvars,means,percentiles)
rownames(results) <- 1:Nvars
colnames(results) <- c("    Root", "    Mean", "  Percentile")

if (display == TRUE) {
	cat("\n\n\nPARALLEL ANALYSIS\n")

	# specification notices
	if (corkind=='pearson')    ctype <- 'pearson' 
	if (corkind=='kendall')    ctype <- 'kendall' 
	if (corkind=='spearman')   ctype <- 'spearman' 
	cat("\nType of correlations specified for the random data eigenvalues: ", ctype)

	if (extract=='PCA')    cat("\n\nExtraction Method: Principal Components\n") 
	if (extract=='PAF')    cat("\n\nExtraction Method: Common Factor Analysis\n")
	if (extract=='image')  cat("\n\nExtraction Method: Image Factor Extraction\n")

	cat("\nVariables  = ", Nvars) 
	cat("\nCases      = ", Ncases) 
	cat("\nNdatasets  = ", Ndatasets) 
	cat("\nPercentile = ", percentile, "\n\n") 

	print(round(results,3))
}

parOutput <- list( eigenvalues=results  )

return(invisible(parOutput))

}
