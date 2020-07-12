
ROOTFIT <- function (data, corkind='pearson', Ncases=NULL, extract='PAF', verbose = 'TRUE') {

Nvars <- ncol(data)

# determine whether data is a correlation matrix
if (nrow(data) == Nvars) {
	if (all(diag(data==1))) {datakind = 'correlations'}
} else{datakind = 'notcorrels'} 
		
if (datakind == 'correlations') { 
	cormat <- data 	
	ctype <- 'from user'
	if (verbose=='TRUE') cat("\n\nThe entered data is a correlation matrix.")	
	if(is.null(Ncases)) {
		cat('\n\nA correlation matrix was entered as input without specifying the number of cases.\n\n')
		cat('An error message will be generated. Enter a value for Ncases.\n\n\n')
	}
	n.obs <- Ncases
}

if (datakind == 'notcorrels') {
	Ncases <- n.obs <- nrow(data)
	if (anyNA(data) == TRUE) {
		data <- na.omit(data)
		cat('\n\nCases with missing values were found and removed from the data matrix.\n\n')
	}
	if (corkind=='pearson')     cormat <- cor(data, method="pearson");  ctype <- 'Pearson' 
	if (corkind=='kendall')     cormat <- cor(data, method="kendall");  ctype <- 'Kendall' 
	if (corkind=='spearman')    cormat <- cor(data, method="spearman"); ctype <- 'Spearman' 
	if (corkind=='polychoric')  cormat <- POLYCHORIC_R(data);           ctype <- 'Polychoric' 
}

 
# factor extraction
if (extract=='PAF')    fm <- 'pa' 
if (extract=='PCA')    fm <- 'pc' 
if (extract=='ML')     fm <- 'mle' 


if (fm=='mle' | fm=='pa') {

	fits <- matrix(-9999,20,12)

	fits[,1] <- 1:20

	evals <- eigen(cormat)$values

	if (min(evals) < 0) {
	     cat('\n\nThe correlation matrix is not positive definite, which means that')
	     cat('\nthere are eigenvalues less than zero. Expect errors.\n')
	}

	fits[,2] <- evals[1:20]


#	for (root in 1:(Nvars)) {
	for (root in 1:20) {

		dof <- 0.5 * ((Nvars - root)^2 - Nvars - root)
		if (dof < 1) {
		fits <- fits[1:(root-1),]
		if (verbose == 'TRUE') {
			cat('\n\nThe degrees of freedom for the model with',
			    root,'factors is < 1 and so the procedure was stopped.\n\n\n')}
		break
		}


		#faOUT <- factanal(covmat=cormat, n.obs=Ncases, factors=root, rotation='none', maxit=1000)
       
		faOUT <- psych::fa(cormat, n.obs=n.obs, nfactors=root, rotate="none", fm=fm, max.iter=1000, SMC=TRUE)


		chisq <- faOUT$STATISTIC
		df <- faOUT$dof
		rmsea <- sqrt(max(((chisq-df)/(n.obs-1)),0)/df)
		pvalue <- faOUT$PVAL

		TLI <- faOUT$TLI
		BIC <- faOUT$BIC

		chisqNULL <- faOUT$null.chisq
		dfNULL <- faOUT$null.dof
		#faOUT$null.model

		CFI <- ((chisqNULL - dfNULL) - (chisq - df)) / (chisqNULL - dfNULL)

		# MacDonald & Marsh (1990) MFI = an absolute fit index that does not depend  
		# on comparison with another model  (T&F, 2001, p 700)
		# MFI <- exp (-.5 * ( (chisq - df) / Ncases))


		fitRevelle <- psych::factor.fit(cormat,faOUT$loadings)

		residuals <- psych::factor.residuals(cormat,faOUT$loadings) 
		residuals <- as.matrix(residuals[upper.tri(residuals)])
		mnsqdresid <- mean(residuals^2) # mean of the off-diagonal squared residuals (as in Waller's MicroFact)
		rmsr <- sqrt(mean(residuals^2)) # rmr is perhaps the more common term for this stat
		# no srmsr computation because it requires the SDs for the variables in the matrix

		# GFI from Waller's MicroFact: 1 - (mean-squared residual / mean-squared correlation)
		correls <- residuals <- as.matrix(cormat[upper.tri(cor(data))])
		mnsqdcorrel <- mean(correls^2) 
		GFI <- 1 - (mnsqdresid / mnsqdcorrel)


		# # AIC Akaike Information Criteria (T&F, 2001, p 700)
		# # not on a 0-1 scale; & the value/formula varies across software
		# AIC <- chisq - 2 * df

		# # CAIC Consistent Akaike Information Criteria (T&F, 2001, p 700)
		# # not on a 0-1 scale; & the value/formula varies across software
		# CAIC <- chisq - (log(Ncases) + 1) * df


		# fits[root,3:16] <- cbind(chisq, df, pvalue, rmsea, fitRevelle, mnsqdresid, rmsr, GFI, MFI, CFI, TLI, BIC, AIC, CAIC)
		# colnames(fits) <- c("Root","  Eigenvalue","  Chi-Squared","   df","       p","   RMSEA","   fitRevelle",
	    # "  mnsqdresid", "     rmsr","      GFI","     MFI","     CFI","     TLI","     BIC","     AIC","     CAIC")
		# rownames(fits) <- matrix((""),nrow(fits),1)


		# fits[root,3:14] <- cbind(chisq, df, pvalue, rmsea, fitRevelle, mnsqdresid, rmsr, GFI, MFI, CFI, TLI, BIC)
		# colnames(fits) <- c("Root","  Eigenvalue","  Chi-Squared","   df","       p","   RMSEA","   fitRevelle",
		# "  mnsqdresid", "     rmsr","      GFI","     MFI","     CFI","     TLI","     BIC")
		# rownames(fits) <- matrix((""),nrow(fits),1)

		fits[root,3:12] <- cbind(chisq, df, pvalue, rmsea, TLI, CFI, rmsr, GFI, BIC, fitRevelle)
		colnames(fits) <- c("Root","  Eigenvalue","  Chi-Squared","   df","       p","   RMSEA",
				       "     TLI","     CFI","     RMSR","      GFI","     BIC", "   Revelle.fit")
		rownames(fits) <- matrix((""),nrow(fits),1)
	}	
}


if (fm=='pc') {

	fits <- matrix(-9999,Nvars,5)

	fits[,1] <- 1:Nvars

	evals <- eigen(cormat)$values

	if (min(evals) < 0) {
	     cat('\n\nThe correlation matrix is not positive definite, which means that')
	     cat('\nthere are eigenvalues less than zero. Expect errors.\n')
	}
	
	fits[,2] <- evals

	for (root in 1:(Nvars)) {
       
		pcaout  <- psych::principal(cormat, n.obs=n.obs, nfactors=root, rotate="none", fm=fm)

		fitRevelle <- pcaout$fit
		fitRevelle.off <- pcaout$fit.off

		residuals <- pcaout$residual
		residuals <- as.matrix(residuals[upper.tri(residuals)])
		mnsqdresid <- mean(residuals^2) # mean of the off-diagonal squared residuals (as in Waller's MicroFact)
		rmsr <- sqrt(mean(residuals^2)) # rmr is perhaps the more common term for this stat
		# no srmsr computation because it requires the SDs for the variables in the matrix

		# GFI from Waller's MicroFact: 1 - mean-squared residual / mean-squared correlation
		correls <- residuals <- as.matrix(cormat[upper.tri(cor(data))])
		mnsqdcorrel <- mean(correls^2) 
		GFI <- 1 - (mnsqdresid / mnsqdcorrel)


		# fits[root,3:7] <- cbind(fitRevelle, fitRevelle.off, mnsqdresid, rmsr, GFI)
		# colnames(fits) <- c("Root","  Eigenvalue","   Revelle.fit","   fitRevelle.off",
		# "  mnsqdresid", "     RMSR","      GFI")
		# rownames(fits) <- matrix((""),nrow(fits),1)


		fits[root,3:5] <- cbind(rmsr, GFI, fitRevelle)
	}
	colnames(fits) <- c("Root","  Eigenvalue","     RMSR","      GFI","   Revelle.fit")
	rownames(fits) <- matrix((""),nrow(fits),1)
}


if (verbose == 'TRUE') {
	
	cat('\n\nFit coefficients for N-factor solutions\n\n')

	cat("\nType of correlations used in the analyses: ", ctype)

	if (extract=='PCA')   { cat("\n\nExtraction Method: Principal Components\n") }
	if (extract=='PAF')   { cat("\n\nExtraction Method: Common Factor Analysis\n")}
	if (extract=='ML')    { cat("\n\nExtraction Method: Maximum Likelihood Estimation\n") } 

	cat("\nThe number of cases      = ", Ncases, "\n")
	cat("\nThe number of variables  = ", Nvars,  "\n\n\n") 

	print(round(fits,3))
}


if (verbose == 'TRUE' & (fm=='mle' | fm=='pa')) {

	cat("\n\n\nEigenvalue\n\n")
	
	cat("An eigenvalue is the variance of the factor. There are as many eigenvalues for a\n")
	cat("correlation or covariance matrix as there are variables in the matrix. The sum\n")
	cat("of the eigenvalues is equal to the number of variables. An eigenvalue of one\n")
	cat("means that a factor explains as much variance as one variable.\n\n\n")
	
	cat("TLI -- Tucker Lewis Index (incremental fit)\n\n")
	
	cat("The Tucker-Lewis index, TLI, is also sometimes called the non-normed fit index,\n")
	cat("NNFI, or the Bentler-Bonett non-normed fit index, or RHO2. The TLI penalizes for\n")
	cat("model complexity.\n\n")
	
	cat("Schermelleh-Engel (2003): 'The (TLI or) NNFI ranges in general from zero to one,\n")
	cat("but as this index is not normed, values can sometimes leave this range, with\n")
	cat("higher (TLI or) NNFI values indicating better fit. A rule of thumb for this\n")
	cat("index is that .97 is indicative of good fit relative to the independence model,\n")
	cat("whereas values greater than .95 may be interpreted as an acceptable fit. An\n")
	cat("advantage of the (TLI or) NNFI is that it is one of the fit indices less\n")
	cat("affected by sample size (Bentler, 1990; Bollen, 1990; Hu & Bentler, 1995,\n")
	cat("1998).'\n\n")
	
	cat("Kenny (2015): 'The TLI (and the CFI) depends on the average size of the\n")
	cat("correlations in the data. If the average correlation between variables is not\n")
	cat("high, then the TLI will not be very high.'\n\n\n")
	
	cat("CFI - Comparative Fit Index (incremental fit)\n\n")
	
	cat("Schermelleh-Engel (2003): 'The CFI ranges from zero to one with higher values\n")
	cat("indicating better fit. A rule of thumb for this index is that .97 is indicative\n")
	cat("of good fit relative to the independence model, while values greater than .95\n")
	cat("may be interpreted as an acceptable fit. Again a value of .97 seems to be more\n")
	cat("reasonable as an indication of a good model fit than the often stated cutoff\n")
	cat("value of .95. Comparable to the NNFI, the CFI is one of the fit indices less\n")
	cat("affected by sample size.'\n\n")
	
	cat("Hooper (2008): 'A cut-off criterion of CFI _ 0.90 was initially advanced\n")
	cat("however, recent studies have shown that a value greater than 0.90 is needed in\n")
	cat("order to ensure that misspecified models are not accepted (Hu and Bentler,\n")
	cat("1999). From this, a value of CFI _ 0.95 is presently recognised as indicative of\n")
	cat("good fit (Hu and Bentler, 1999). Today this index is included in all SEM\n")
	cat("programs and is one of the most popularly reported fit indices due to being one\n")
	cat("of the measures least effected by sample size (Fan et al, 1999).'\n\n")
	
	cat("Kenny (2015): 'Because the TLI and CFI are highly correlated only one of the two\n")
	cat("should be reported. The CFI is reported more often than the TLI, but I think the\n")
	cat("CFI's penalty for complexity of just 1 is too low and so I prefer the TLI even\n")
	cat("though the CFI is reported much more frequently than the TLI.'\n\n\n")
	
	cat("RMSEA - Root Mean Square Error of Approximation (absolute fit)\n\n")
	
	cat("Schermelleh-Engel (2003): 'The Root Mean Square Error of Approximation (RMSEA;\n")
	cat("Steiger, 1990) is a measure of approximate fit in the population and is\n")
	cat("therefore concerned with the discrepancy due to approximation. _ Steiger (1990)\n")
	cat("as well as Browne and Cudeck (1993) define a 'close fit' as a RMSEA value less\n")
	cat("than or equal to .05. According to Browne and Cudeck (1993), RMSEA values _ .05\n")
	cat("can be considered as a good fit, values between .05 and .08 as an adequate fit,\n")
	cat("and values between .08 and .10 as a mediocre fit, whereas values > .10 are not\n")
	cat("acceptable. Although there is general agreement that the value of RMSEA for a\n")
	cat("good model should be less than .05, Hu and Bentler (1999) suggested an RMSEA of\n")
	cat("less than .06 as a cutoff criterion.'\n\n")
	
	cat("Kenny (2015): 'The measure is positively biased (i.e., tends to be too large)\n")
	cat("and the amount of the bias depends on smallness of sample size and df, primarily\n")
	cat("the latter. The RMSEA is currently the most popular measure of model fit. _\n")
	cat("MacCallum, Browne and Sugawara (1996) have used 0.01, 0.05, and 0.08 to indicate\n")
	cat("excellent, good, and mediocre fit respectively. However, others have suggested\n")
	cat("0.10 as the cutoff for poor fitting models. These are definitions for the\n")
	cat("population. That is, a given model may have a population value of 0.05 (which\n")
	cat("would not be known), but in the sample it might be greater than 0.10.'\n\n")
	
	cat("Hooper (2008): 'In recent years it has become regarded as 'one of the most\n")
	cat("informative fit indices' (Diamantopoulos and Siguaw, 2000: 85) due to its\n")
	cat("sensitivity to the number of estimated parameters in the model. In other words,\n")
	cat("the RMSEA favours parsimony in that it will choose the model with the lesser\n")
	cat("number of parameters.'\n\n\n")
	
	cat("RMSR -- Root Mean Square Residual (absolute fit)\n\n")
	
	cat("RMSR (or perhaps more commonly, RMR) is an index of the overall badness-of-fit.\n")
	cat("It is the square root of the mean of the squared residuals (the residuals being\n")
	cat("the simple differences between original correlations and the correlations\n")
	cat("implied by the N-factor model). RMSR is 0 when there is perfect model fit. A\n")
	cat("value less than .08 is generally considered a good fit. A standardized version\n")
	cat("of the RMSR is often recommended over the RMSR in structural equation modeling\n")
	cat("analyses. This is because the values in covariance matrices are scale-dependent.\n")
	cat("However, the RMSR coefficient that is provided in this package is based on\n")
	cat("correlation coefficients (not covariances) and therefore does not have this\n")
	cat("problem. The RMSR is provided rather than the standardized version because the\n")
	cat("variable standard deviations are not required for the RMSR.\n\n\n")
	
	cat("GFI (absolute fit)\n\n")
	
	cat("Waller, MicroFACT 2.0 Manual\n")
	cat("(http://www.psych.umn.edu/faculty/waller/downloads.htm): 'A GFI index is equal\n")
	cat("to 1.0 - mean-squared residual/mean-squared correlation'\n\n\n")
	
	cat("BIC -- Bayesian Information Criterion (degree of parsimony fit index) Kenny\n")
	cat("(2015): 'The BIC is a comparative measure of fit and so it is meaningful only\n")
	cat("when two different models are estimated. Lower values indicate a better fit and\n")
	cat("so the model with the lowest BIC is the best fitting model. Whereas the AIC has\n")
	cat("a penalty of 2 for every parameter estimated, the BIC increases the penalty as\n")
	cat("sample size increases. The BIC places a high value on parsimony (perhaps too\n")
	cat("high).'\n\n\n")
	
	cat("Revelle.fit (absolute fit)\n\n")
	
	cat("This is the factor.fit coefficient from Revelle's psych package. From the\n")
	cat("package documentation: 'How well does the factor model reproduce the correlation\n")
	cat("matrix? This fit is a plausible estimate of the amount of reduction in a\n")
	cat("correlation matrix given a factor model. Note that it is sensitive to the size\n")
	cat("of the original correlations. That is, if the residuals are small but the\n")
	cat("original correlations are small, that is a bad fit. The basic factor or\n")
	cat("principal components model is that a correlation or covariance matrix may be\n")
	cat("reproduced by the product of a factor loading matrix times its transpose. One\n")
	cat("simple index of fit is the 1 - sum squared residuals/sum squared original\n")
	cat("correlations. The sums are taken for the off diagonal elements.'\n\n\n")
	
	cat("Hooper, D., Coughlan, J., & Mullen, M. (2008). Structural Equation Modelling:\n")
	cat("Guidelines for Determining Model Fit. Electronic Journal of Business Research\n")
	cat("Methods, 6(1), 53-60.\n\n")
	
	cat("Kenny, D. A. (2015). Measuring model fit. http://davidaKenny.net/cm/fit.htm\n\n")
	
	cat("Schermelleh-Engel, K., & Moosbrugger, H. (2003). Evaluating the Fit of\n")
	cat("Structural Equation Models: Tests of Significance and Descriptive\n")
	cat("Goodness-of-Fit Measures. Methods of Psychological Research Online, Vol.8(2),\n")
	cat("pp. 23-74. http://www.mpr-online.de\n\n\n\n")
}


if (verbose == 'TRUE' & fm=='pc') {
	
	cat("\n\n\nEigenvalue\n\n")
	
	cat("An eigenvalue is the variance of the factor. There are as many eigenvalues for a\n")
	cat("correlation or covariance matrix as there are variables in the matrix. The sum\n")
	cat("of the eigenvalues is equal to the number of variables. An eigenvalue of one\n")
	cat("means that a factor explains as much variance as one variable.\n\n\n")
	
	cat("RMSR -- Root Mean Square Residual (absolute fit)\n\n")
	
	cat("RMSR (or perhaps more commonly, RMR) is an index of the overall badness-of-fit.\n")
	cat("It is the square root of the mean of the squared residuals (the residuals being\n")
	cat("the simple differences between original correlations and the correlations\n")
	cat("implied by the N-factor model). RMSR is 0 when there is perfect model fit. A\n")
	cat("value less than .08 is generally considered a good fit. A standardized version\n")
	cat("of the RMSR is often recommended over the RMSR in structural equation modeling\n")
	cat("analyses. This is because the values in covariance matrices are scale-dependent.\n")
	cat("However, the RMSR coefficient that is provided in this package is based on\n")
	cat("correlation coefficients (not covariances) and therefore does not have this\n")
	cat("problem. The RMSR is provided rather than the standardized version because the\n")
	cat("variable standard deviations are not required for the RMSR.\n\n\n")
	
	cat("GFI (absolute fit)\n\n")
	
	cat("Waller, MicroFACT 2.0 Manual\n")
	cat("(http://www.psych.umn.edu/faculty/waller/downloads.htm): 'A GFI index is equal\n")
	cat("to 1.0 - mean-squared residual/mean-squared correlation'\n\n\n")
	
	cat("BIC -- Bayesian Information Criterion (degree of parsimony fit index) Kenny\n")
	cat("(2015): 'The BIC is a comparative measure of fit and so it is meaningful only\n")
	cat("when two different models are estimated. Lower values indicate a better fit and\n")
	cat("so the model with the lowest BIC is the best fitting model. Whereas the AIC has\n")
	cat("a penalty of 2 for every parameter estimated, the BIC increases the penalty as\n")
	cat("sample size increases. The BIC places a high value on parsimony (perhaps too\n")
	cat("high).'\n\n\n")
	
	cat("Revelle.fit (absolute fit)\n\n")
	
	cat("This is the factor.fit coefficient from Revelle's psych package. From the\n")
	cat("package documentation: 'How well does the factor model reproduce the correlation\n")
	cat("matrix? This fit is a plausible estimate of the amount of reduction in a\n")
	cat("correlation matrix given a factor model. Note that it is sensitive to the size\n")
	cat("of the original correlations. That is, if the residuals are small but the\n")
	cat("original correlations are small, that is a bad fit. The basic factor or\n")
	cat("principal components model is that a correlation or covariance matrix may be\n")
	cat("reproduced by the product of a factor loading matrix times its transpose. One\n")
	cat("simple index of fit is the 1 - sum squared residuals/sum squared original\n")
	cat("correlations. The sums are taken for the off diagonal elements.\n\n\n")
}

return(invisible(fits))

}




