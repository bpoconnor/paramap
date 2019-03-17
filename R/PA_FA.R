

# CFA / PAF  (Bernstein p 189; smc = from Bernstein p 104)

PA_FA <- function (data, corkind='pearson', nfactors=2, iterpaf=100, rotate='promax', ppower=3, display=TRUE ) {

cnoms <- colnames(data) # get colnames

# determine whether data is a correlation matrix
if ( nrow(data) == ncol(data) ) {
	if ( all(diag(data==1)) ) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}


if (datakind == 'correlations')  cormat <- data 

if (datakind == 'notcorrels') {
	Ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
		}
	if (corkind=='pearson')     {cormat <- cor(data, method="pearson");  ctype <- 'Pearson'}
	if (corkind=='kendall')     {cormat <- cor(data, method="kendall");  ctype <- 'Kendall'}
	if (corkind=='spearman')    {cormat <- cor(data, method="spearman"); ctype <- 'Spearman'} 
	if (corkind=='polychoric')  {cormat <- POLYCHORIC_R(data);            ctype <- 'Polychoric'}
	}

converge  <- .001
rpaf <- as.matrix(cormat)
smc <- 1 - (1 / diag(solve(rpaf)) )
for (iter in 1:(iterpaf + 1)) {

	diag(rpaf) <- smc # putting smcs on the main diagonal of r

	eigval <-  diag((eigen(rpaf) $values))
	# substituting zero for negative eigenvalues
	for (luper in 1:nrow(eigval)) { if ( eigval[luper,luper] < 0 ) { eigval[luper,luper] <- 0 }}

eigvect <- eigen(rpaf) $vectors

if (nfactors == 1) {
loadings <- eigvect[,1:nfactors] * sqrt(eigval[1:nfactors,1:nfactors])
communal <- loadings^2
}else {
loadings <- eigvect[,1:nfactors] %*% sqrt(eigval[1:nfactors,1:nfactors])
communal <- rowSums(loadings^2) }

if ( max(max(abs(communal-smc))) < converge) { break }
if ( max(max(abs(communal-smc))) >= converge  & iter < iterpaf) { smc <- communal }

}
loadings <- as.matrix(loadings)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste("  Factor ", 1:nfactors, sep="") )

evalpaf  <-  cbind(diag(eigval))
rownames(evalpaf) <- 1:nrow(evalpaf)
colnames(evalpaf) <- "Eigenvalues"

if (rotate=='none')  { pafOutput <- list( eigenvalues=evalpaf, loadingsNOROT=loadings ) }

if (rotate=='promax' | rotate=='varimax') {
	
	if (nfactors==1) { 
		pafOutput <- list( eigenvalues=evalpaf, loadingsNOROT=loadings, loadingsROT=loadings, structure=loadings, pattern=loadings ) 
		} 

	if (nfactors > 1) {
		if (rotate=='varimax') { 
			loadingsROT <- paramap::VARIMAX(loadings,display=FALSE)
			pafOutput <- list( eigenvalues=evalpaf, loadingsNOROT=loadings, loadingsROT=loadingsROT )  
			} 
	  # if (rotate=='varimax') { 
			# loadingsROT <- stats::VARIMAX(loadings)
		    # pafOutput <- list( eigenvalues=evalpaf, loadingsNOROT=loadings, loadingsROT=loadingsROT$loadings )  
			# } 
		if (rotate=='promax')  { 
			loadingsROT <- paramap::PROMAX(loadings,display=FALSE)
			pafOutput <- list( eigenvalues=evalpaf, structure=loadingsROT$structure, pattern=loadingsROT$pattern, correls=loadingsROT$correls ) 
			}
}}

if (display == TRUE) {
	cat("\n\nPrincipal Axis Factor Analysis\n\n")
	cat("\nSpecified kind of correlations for this analysis: ", ctype, "\n\n")
	if ( max(max(abs(communal-smc))) < converge) {
		cat("\nPAF converged in iterations = ", iter, "\n\n")	
		print(round(evalpaf,2))
		# cat("\n\nCommunalities: \n")
		# print(round(communal,2))
		cat("\nUnrotated PAF Loadings:\n\n")
		print(round(cbind(loadings, communal),2))
		} else { cat("\nPAF did not converge in the following number of iterations:  ", (iter-1))
		}
	if (nfactors==1) { cat("\n\nNo rotation because there is only one factor\n\n") }
	if (nfactors > 1) {
		if (rotate=='none')    {cat("\n\nRotation Procedure:  No Rotation")}
		if (rotate=='varimax') {cat("\n\nVarimax Rotated Loadings:\n\n"); print(round(loadingsROT,2));cat("\n\n")}
		if (rotate=='promax')  { 
		cat("\n\nPromax Rotation Structure Matrix:\n\n");    print(round(loadingsROT$structure,2));cat("\n")
		cat("\n\nPromax Rotation Pattern Matrix:\n\n");      print(round(loadingsROT$pattern,2));cat("\n")
		cat("\n\nPromax Rotation Factor Correlations:\n\n"); print(round(loadingsROT$correls,2));cat("\n\n")
}}}

return(invisible(pafOutput))

}



