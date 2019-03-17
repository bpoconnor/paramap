

NEVALSGT1 <- function (data, corkind='pearson', display=FALSE) {

# Number of eigenvalues > 1

nvars  <- ncol(data)

# determine whether data is a correlation matrix
if ( nrow(data) == ncol(data) ) {
	if ( max(diag(data)) == 1 & min(diag(data)) == 1 ) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}

if (datakind == 'correlations')  rdata <- data 

if (datakind == 'notcorrels') {
	ncases <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}
	if (corkind=='pearson')     rdata <- cor(data, method="pearson") 
	if (corkind=='kendall')     rdata <- cor(data, method="kendall") 
	if (corkind=='spearman')    rdata <- cor(data, method="spearman") 
	if (corkind=='polychoric')  rdata <- POLYCHORIC_R(data) 
}


eigvals <- cbind(eigen(rdata) $values)

nfnevalsgt1 <- 0

for (nev in 1:nrow(eigvals)) {if (eigvals[nev,] > 1) nfnevalsgt1 <- nfnevalsgt1 + 1}


if (display == TRUE) { 

if (datakind == 'correlations') cat("\n\n The entered data is a correlation matrix.") 

if (datakind == 'notcorrels') {
	cat("\nNumber of cases in the data file =       ", ncases)
	cat("\nNumber of variables in the data file =   ", nvars)

	# specification notices
	if (corkind=='pearson')    {cat("\nCorrelations to be Analyzed: Pearson")}
	if (corkind=='kendall')    {cat("\nCorrelations to be Analyzed: Kendall")}
	if (corkind=='spearman')   {cat("\nCorrelations to be Analyzed: Spearman")}
	if (corkind=='polychoric') {cat("\nCorrelations to be Analyzed: Polychoric")}
}

cat("\n\n\n")
colnames(eigvals) <- 'Eigenvalues'
rownames(eigvals) <- 1:length(eigvals)
print(round(eigvals,5))

cat('\n\nThe number of eigenvalues greater than one = ', nfnevalsgt1, '\n\n')
}

return(invisible(nfnevalsgt1))

cat("\n\n")

}
