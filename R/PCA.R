
# Principal Components Analysis

PCA <- function (data, corkind='pearson', nfactors=2, rotate='promax', ppower=3, display=TRUE) {

cnoms <- colnames(data) # get colnames

# determine whether data is a correlation matrix
if ( nrow(data) == ncol(data) ) {
	if ( all(diag(data==1)) ) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}


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

eigval <- diag(eigen(cormat) $values)
eigvect <- eigen(cormat) $vectors
if (nfactors == 1) {loadings <- eigvect[,1:nfactors] * sqrt(eigval[1:nfactors,1:nfactors])
}else {loadings <- eigvect[,1:nfactors] %*% sqrt(eigval[1:nfactors,1:nfactors])}
loadings <- as.matrix(loadings)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste("  Factor ", 1:nfactors, sep="") )



evalpca  <-  cbind(diag(eigval))
rownames(evalpca) <- 1:nrow(evalpca)
colnames(evalpca) <- "Eigenvalues"

if (rotate=='none')  { pcaOutput <- list( eigenvalues=evalpca, loadingsNOROT=loadings ) }

if (rotate=='promax' | rotate=='varimax') {
	
	if (nfactors==1)  pcaOutput <- list( eigenvalues=evalpca, loadingsNOROT=loadings, loadingsROT=loadings, structure=loadings, pattern=loadings )  

	if (nfactors > 1) {
		if (rotate=='varimax') { 
			loadingsROT <- paramap::VARIMAX(loadings,display=FALSE)
			pcaOutput <- list( eigenvalues=evalpca, loadingsNOROT=loadings, loadingsROT=loadingsROT ) 
			}  
		if (rotate=='promax')  { 
			loadingsROT <- paramap::PROMAX(loadings,display=FALSE)
			pcaOutput <- list( eigenvalues=evalpca, structure=loadingsROT$structure, pattern=loadingsROT$pattern, correls=loadingsROT$correls ) 
			}
}}

if (display == TRUE) {
	cat("\n\nPrincipal Components Analysis\n\n")
	cat("\nSpecified kind of correlations for this analysis: ", ctype, "\n\n\n")
	print(round(evalpca,2))
	cat("\n\nUnrotated PCA Loadings:\n\n")
	print(round(loadings[,1:nfactors],2));cat("\n")
		if (nfactors==1) { cat("\n\nNo rotation because there is only one component\n\n") }
		if (nfactors > 1) {
			if (rotate=='none')    {cat("\n\nRotation Procedure:  No Rotation")}
			if (rotate=='varimax') {cat("\n\nVarimax Rotated Loadings:\n\n"); print(round(loadingsROT,2));cat("\n\n") }
			if (rotate=='promax')  { 
				cat("\n\nPromax Rotation Structure Matrix:\n\n");    print(round(loadingsROT$structure,2));cat("\n")
				cat("\n\nPromax Rotation Pattern Matrix:\n\n");      print(round(loadingsROT$pattern,2));cat("\n")
				cat("\n\nPromax Rotation Factor Correlations:\n\n"); print(round(loadingsROT$correls,2));cat("\n\n")
}}}

return(invisible(pcaOutput))

}
