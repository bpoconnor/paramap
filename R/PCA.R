
# Principal Components Analysis

PCA <- function (data, corkind='pearson', Nfactors=NULL, rotate='promax', ppower=3, verbose=TRUE) {

cnoms <- colnames(data) # get colnames

# determine whether data is a correlation matrix
if (nrow(data) == ncol(data)) {
	if (all(diag(data==1))) {datakind = 'correlations'}
} else{ datakind = 'notcorrels'}


if (datakind == 'correlations')  {
	cormat <- data 
	ctype <- 'from user'
}

if (datakind == 'notcorrels') {
	Ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}
	if (corkind=='pearson')     {cormat <- cor(data, method="pearson");  ctype <- 'Pearson'}
	if (corkind=='kendall')     {cormat <- cor(data, method="kendall");  ctype <- 'Kendall'}
	if (corkind=='spearman')    {cormat <- cor(data, method="spearman"); ctype <- 'Spearman'} 
	if (corkind=='polychoric')  {cormat <- POLYCHORIC_R(data);           ctype <- 'Polychoric'}
}

if (is.null(Nfactors)) {		
	nfactsMAP <- MAP(cormat, verbose=FALSE)
	Nfactors <- nfactsMAP$nfMAP
	NfactorsWasNull <- TRUE
} else {NfactorsWasNull <- FALSE}

eigval <- diag(eigen(cormat) $values)
eigvect <- eigen(cormat) $vectors
if (Nfactors == 1) {loadings <- eigvect[,1:Nfactors] * sqrt(eigval[1:Nfactors,1:Nfactors])
}else {loadings <- eigvect[,1:Nfactors] %*% sqrt(eigval[1:Nfactors,1:Nfactors])}
loadings <- as.matrix(loadings)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste("  Factor ", 1:Nfactors, sep=""))


evalpca  <-  cbind(diag(eigval))
rownames(evalpca) <- 1:nrow(evalpca)
colnames(evalpca) <- "Eigenvalues"

if (rotate=='none')  { pcaOutput <- list(eigenvalues=evalpca, loadingsNOROT=loadings) }

if (rotate=='promax' | rotate=='varimax') {
	
	if (Nfactors==1)  pcaOutput <- list(eigenvalues=evalpca, loadingsNOROT=loadings, 
	                                    loadingsROT=loadings, structure=loadings, pattern=loadings)  

	if (Nfactors > 1) {
		if (rotate=='varimax') { 
			loadingsROT <- paramap::VARIMAX(loadings,verbose=FALSE)
			pcaOutput <- list(eigenvalues=evalpca, loadingsNOROT=loadings, loadingsROT=loadingsROT) 
			}  
		if (rotate=='promax')  { 
			loadingsROT <- paramap::PROMAX(loadings,verbose=FALSE)
			pcaOutput <- list(eigenvalues=evalpca, structure=loadingsROT$structure, 
			                  pattern=loadingsROT$pattern, correls=loadingsROT$correls) 
		}
	}
}

if (verbose == TRUE) {
	cat("\n\nPrincipal Components Analysis\n\n")
	cat("\nSpecified kind of correlations for this analysis: ", ctype, "\n\n")
	if (NfactorsWasNull == TRUE) {
		cat('\nNfactors was not specified and so the MAP test was conducted to determine')
		cat('\nthe number of factors to extract: Nfactors =', Nfactors,'\n\n\n')		
	} else if (NfactorsWasNull == FALSE) {
		cat('\nThe specified number of factors to extract =', Nfactors,'\n\n\n')
	}
	print(round(evalpca,2))
	cat("\n\nUnrotated PCA Loadings:\n\n")
	print(round(loadings[,1:Nfactors],2));cat("\n")
		if (Nfactors==1) { cat("\n\nNo rotation because there is only one component\n\n") }
		if (Nfactors > 1) {
			if (rotate=='none')    {cat("\n\nRotation Procedure:  No Rotation")}
			if (rotate=='varimax') {cat("\n\nVarimax Rotated Loadings:\n\n"); print(round(loadingsROT,2));cat("\n\n") }
			if (rotate=='promax')  { 
				cat("\n\nPromax Rotation Structure Matrix:\n\n");    print(round(loadingsROT$structure,2));cat("\n")
				cat("\n\nPromax Rotation Pattern Matrix:\n\n");      print(round(loadingsROT$pattern,2));cat("\n")
				cat("\n\nPromax Rotation Factor Correlations:\n\n"); print(round(loadingsROT$correls,2));cat("\n\n")
			}
		}
}

return(invisible(pcaOutput))

}
