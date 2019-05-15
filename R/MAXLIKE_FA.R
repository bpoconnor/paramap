

# Maximum likelihood factor analysis algorithm -- Marcus, 1993

MAXLIKE_FA <- function (data, corkind='pearson', Nfactors=NULL, tolerml=.001, iterml=100, rotate='promax', ppower=3, display=TRUE ) {

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
	if (corkind=='polychoric')  {cormat <- POLYCHORIC_R(data);            ctype <- 'Polychoric'}
}

if (is.null(Nfactors)) {		
	nfactsMAP <- MAP(cormat, display='no')
	Nfactors <- nfactsMAP$nfMAP
	NfactorsWasNull <- TRUE
}

Rho <- cormat
k <- Nfactors
p <- nrow(Rho)

# preliminary singular value decomposition of Rho
L <- diag(svd(Rho) $d)
A <- svd(Rho) $u

if (Nfactors == 1) { A1 <- A[,1:k] * sqrt(L[1:k,1:k])
}else { A1 <- A[,1:k] %*% sqrt(L[1:k,1:k]) }   # Prin. Comp. loadings 

Uni <- diag(diag(Rho-A1%*%t(A1)))     # Uniqueness matrix
Rh1 <- sqrt(solve(Uni))%*%(Rho-Uni)%*%sqrt(solve(Uni)) # Matrix to iterate

# First estimate of Maximum Likelihood Loadings
L <- diag(svd(Rh1) $d)
A <- svd(Rh1) $u

A1 <- sqrt(Uni)%*%A[,1:k]%*%sqrt(L[1:k,1:k])

check <- tolerml
for (i in 1:iterml) {
   Uni <- diag(diag(Rho-A1%*%t(A1)))
   Rh1 <- sqrt(solve(Uni))%*%(Rho-Uni)%*%sqrt(solve(Uni))
   L <- diag(svd(Rh1) $d)
   A <- svd(Rh1) $u
   A2 <- sqrt(Uni)%*%A[,1:k]%*%sqrt(L[1:k,1:k])

   if (max(max(abs(A1-A2))) < check) {break}
   A1 <- A2 }

FacVar <- diag(t(A1)%*%A1)

Com <- as.matrix(diag(A1%*%t(A1)))  # communalities
rownames(Com) <- cnoms
colnames(Com) <- "Communalities"

Uniq <- matrix(1,p,1) - Com

Resid <- Rho - A1%*%t(A1)

loadings <- as.matrix(A1)
rownames(loadings) <- cnoms
colnames(loadings) <-  c(paste("  Factor ", 1:Nfactors, sep="") )

evalmax <- as.matrix(diag(L))
rownames(evalmax) <- 1:nrow(evalmax)
colnames(evalmax) <- "Eigenvalues"

if (rotate=='none')  { maxlikeOutput <- list( eigenvalues=evalmax, loadingsNOROT=loadings ) }

if (rotate=='promax' | rotate=='varimax') {

	if (Nfactors==1) { maxlikeOutput <- list( eigenvalues=evalmax, loadingsNOROT=loadings, loadingsROT=loadings, structure=loadings, pattern=loadings ) } 

	if (Nfactors > 1) {
		if (rotate=='varimax') { 
			loadingsROT <- paramap::VARIMAX(loadings,display=FALSE)
			maxlikeOutput <- list( eigenvalues=evalmax, loadingsNOROT=loadings, loadingsROT=loadingsROT ) 
			} 
		if (rotate=='promax')  { 
			loadingsROT <- paramap::PROMAX(loadings,display=FALSE)
			maxlikeOutput <- list( eigenvalues=evalmax, structure=loadingsROT$structure, pattern=loadingsROT$pattern, correls=loadingsROT$correls ) 
			}
}}

if (display == TRUE) {
	cat("\n\nMaximum likelihood factor analysis:\n\n")
	cat("\nSpecified kind of correlations for this analysis: ", ctype, "\n")
	if (NfactorsWasNull <- TRUE) {
		cat('\nNfactors was not specified and so the MAP test was conducted to determine')
		cat('\nthe number of factors to extract: Nfactors =', Nfactors,'\n\n')		
	}
	cat("\n\nNumber of iterations = ", i, "\n\n")	
	print(round(evalmax,2));cat("\n\n")
	print(round(Com,2))
	cat("\n\nUnrotated Maximum Likelihood Loadings:\n\n")
	print(round(loadings[,1:Nfactors],2));cat("\n\n")
	if (Nfactors==1) { cat("\n\nNo rotation because there is only one factor\n\n") }
	if (Nfactors > 1) {
		if (rotate=='none')    {cat("\n\nRotation Procedure:  No Rotation") }
		if (rotate=='varimax') {cat("\n\nVarimax Rotated Loadings:\n\n"); print(round(loadingsROT,2));cat("\n\n") }
		if (rotate=='promax')  { 
		cat("\n\nPromax Rotation Structure Matrix:\n\n");    print(round(loadingsROT$structure,2));cat("\n")
		cat("\n\nPromax Rotation Pattern Matrix:\n\n");      print(round(loadingsROT$pattern,2));cat("\n")
		cat("\n\nPromax Rotation Factor Correlations:\n\n"); print(round(loadingsROT$correls,2));cat("\n\n")
}}}

return(invisible(maxlikeOutput))

}
