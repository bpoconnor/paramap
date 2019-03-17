
# three methods of assessing the factorability of a correlation matrix or raw data set

FACTORABILITY <- function (data, corkind='pearson', Ncases=NULL, display=TRUE) {

cnoms <- colnames(data) # get colnames

# determine whether data is a correlation matrix
if ( nrow(data) == ncol(data) ) {
	if ( all(diag(data==1)) ) {datakind = 'correlations'}} else{ datakind = 'notcorrels'}

if (datakind == 'correlations') {
	cormat <- data 
	if (is.null(Ncases)) {
		Ncases = 200
		cat('\n\n"data" is a correlation matrix but Ncases was not specified, so Ncases was set = 200')		
	}
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


#cat('\n\nThree methods of assessing the factorability of a correlation matrix or raw data set:')


# the determinant of the correlation matrix should be > 0.00001 for factorability
detcor <- det(cormat)
# cat('\n\n\nThe determinant of the correlation matrix should be > 0.00001 for factorability.')
# if (detcor > 0.00001) cat('\n\nThe determinant is',round(detcor,7),'which is > 0.00001, indicating factorability.')
# if (detcor <= 0.00001) cat('\n\nThe determinant is',round(detcor,7),'which is NOT > 0.00001, indicating NON factorability.')


# Bartlett's identity matrix test
# cortest.bartlett(cormat)  # from the psych package, same results
# cat("\n\n\nBartlett's test of whether a correlation matrix is significantly different")
# cat('\nfrom an identity matrix (wherein all of the off-diagonal elements are zero).') 
if (datakind == 'notcorrels') Ncases <- nrow(data)
Nvars  <- ncol(data)
chi2 <- (Ncases - 1 - (2 * Nvars + 5) / 6) * -log(detcor) 
df <- Nvars * (Nvars - 1) / 2
pvalue <- pchisq(chi2, df, lower.tail=F)
# cat('\n\nchisq =',chi2,'    df=',df,'     p =',round(pvalue,8))
# cat('\n\nA significant difference is required for factorability.')



# the Kaiser-Meyer-Olkin measure of sampling adequacy (MSA) -- Kaiser & Rice (1974) 
#cat('\n\n\n\nThe Kaiser-Meyer-Olkin measure of sampling adequacy (MSA):')

# # using the KMO function from the psych package
# KMOobj <- KMO(cormat)
# cat('\n\nOverall measure of sampling adequacy (MSA) =',round(KMOobj$MSA,2),'\n\n')
# itemMSA <- matrix(KMOobj$MSAi,length(KMOobj$MSAi),1)
# rownames(itemMSA) <- cnoms
# colnames(itemMSA) <- '  Item MSA'
# print(round(itemMSA,2))

Rinv <- solve(cormat)
Rpart <- cov2cor(Rinv) 
cormat_sq <- cormat^2
Rpart_sq  <- Rpart^2

# overall KMO
KMOnum <- sum(cormat_sq) - sum(diag(cormat_sq))
KMOdenom <- KMOnum + (sum(Rpart_sq) - sum(diag(Rpart_sq))) 
KMO <- KMOnum / KMOdenom
#cat('\n\nOverall measure of sampling adequacy (MSA) =',round(KMO,2),'\n\n')

# variable KMOs
diag(cormat_sq) <- 0
diag(Rpart_sq)  <- 0
KMOvars <- colSums(cormat_sq)/(colSums(cormat_sq) + colSums(Rpart_sq))
KMOvars <- matrix(KMOvars,length(KMOvars),1)
rownames(KMOvars) <- cnoms
colnames(KMOvars) <- '  Variable MSA'


if (display) {
cat('\n\nThree methods of assessing the factorability of a correlation matrix or raw data set:')

cat('\n\n\nThe determinant of the correlation matrix should be > 0.00001 for factorability.')
if (detcor > 0.00001) cat('\n\nThe determinant is',round(detcor,7),'which is > 0.00001, indicating factorability.')
if (detcor <= 0.00001) cat('\n\nThe determinant is',round(detcor,7),'which is NOT > 0.00001, indicating NON factorability.')

cat("\n\n\nBartlett's test of whether a correlation matrix is significantly different")
cat('\nfrom an identity matrix (wherein all of the off-diagonal elements are zero):')
cat('\n\nchisq =',chi2,'    df=',df,'     p =',round(pvalue,8))
cat('\n\nA significant difference is required for factorability.')

cat('\n\n\n\nThe Kaiser-Meyer-Olkin measure of sampling adequacy (MSA):')

cat('\n\nOverall measure of sampling adequacy (MSA) =',round(KMO,2),'\n\n')

print(round(KMOvars,2))

cat("\n\nKaiser & Rice's (1974) interpretation guidelines for MSA values:")
cat('\n
   KMO >= .9 is marvelous
   KMO in the .80s is mertitorious
   KMO in the .70s is middling
   KMO in the .60s is medicore
   KMO in the .50s is miserable
   KMO < .5 is unacceptable')
cat('\n\nConsider excluding items with KMO values < .5 and then re-run the FACTORABILITY analyses.\n')

cat('\n\nThe overall KMO coefficient indicates the proportion of')
cat('\nvariance in the variables that might be caused by underlying')
cat('\nfactors. If the variables share common factors, then the')
cat('\noverall KMO coefficient should be close to 1.0. The overall')
cat('\nKMO indicates the extent to which there is at least one')
cat('\nlatent factor underlying the variables. The overall KMO')
cat('\nindex is considered particularly meaningful when the cases')
cat('\nto variables ratio is less than 1:5. The KMO coefficient for')
cat('\na variable is a kind of summary index of how much a')
cat('\nvariable overlaps with the other variables.\n\n\n')
}

factOutput <- list( chisq=chi2, df=df, pvalue=pvalue, Rimage=Rpart, KMO=KMO, KMOvars=KMOvars )

return(invisible(factOutput))

}
