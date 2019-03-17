\name{PROMAX}
\alias{PROMAX}
\title{promax rotation}
\description{promax rotation}
\usage{PROMAX(loadings, ppower, display)}
\arguments{
  \item{loadings}{
  \code{}A loading matrix.}
  \item{ppower}{
  \code{}The exponent for the promax target matrix. 'ppower' must be 1 or greater. '4' is a conventional value.}
  \item{display}{
  \code{}Display the results? The options are TRUE or FALSE.}
}
\value{ A list with the following elements:
 \item{structure}{The structure matrix (for promax rotation)}
 \item{pattern}{The pattern matrix (for promax rotation)}
 \item{correls}{The correlations between the factors (for promax rotation)}
}
\author{Brian P. O'Connor}
\examples{
\dontrun{
loadings <- PCA(data_NEOPIR, corkind='pearson', nfactors = 5, rotate='none', display=TRUE)

PROMAX(loadings, ppower = 3, display=TRUE)
}
}
