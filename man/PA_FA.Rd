\name{PA_FA}
\alias{PA_FA}
\title{principal axis (common) factor analysis}
\description{principal axis (common) factor analysis with squared multiple 
	correlations as the initial communality estimates}
\usage{PA_FA(data, corkind, nfactors, iterpaf, rotate, ppower, display)}
\arguments{
  \item{data}{
  \code{}An all-numeric dataframe where the rows are cases & the columns are the 
  variables, or a correlation matrix with ones on the diagonal.The function internally 
  determines whether the data are a correlation matrix.}
  \item{corkind}{
  \code{}The kind of correlation matrix to be used if data is not a correlation matrix. 
  The options are 'pearson', 'kendall', 'spearman', and 'polychoric'. Required only if 
  the entered data is not a correlation matrix.}
  \item{nfactors}{
  \code{}The number of factors to extract.}
  \item{iterpaf}{
  \code{}The maximum number of iterations.}
  \item{rotate}{
  \code{}The factor rotation method. The options are: 'promax', 'varimax', and 'none'.}
  \item{ppower}{ 
  \code{}The power value to be used in a promax rotation (required only if
         rotate = 'promax'). Suggested value: 3}
   \item{display}{
  \code{}Display the results? The options are TRUE or FALSE.}
}
\value{ A list with the following elements:
 \item{eigenvalues}{The eigenvalues}
 \item{loadingsNOROT}{The unrotated factor loadings}
 \item{loadingsROT}{The rotated factor loadings (for varimax rotation)}
 \item{structure}{The structure matrix (for promax rotation)}
 \item{pattern}{The pattern matrix (for promax rotation)}
 \item{correls}{The correlations between the factors (for promax rotation)}
}
\author{Brian P. O'Connor }
\examples{
PA_FA(data_RSE, corkind="pearson", nfactors = 2, iterpaf = 50, 
      rotate='promax', ppower=3, display="yes")
}
