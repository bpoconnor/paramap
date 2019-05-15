\name{VARIMAX}
\alias{VARIMAX}
\title{varimax rotation}
\description{varimax rotation}
\usage{VARIMAX(loadings, display)}
\arguments{
  \item{loadings}{
  \code{}A loading matrix.}
  \item{display}{
  \code{}Display the results? The options are TRUE or FALSE.}
}
\value{The varimax-rotated loadings}
\author{Brian P. O'Connor }
\examples{
\dontrun{
loadings <- PCA(data_NEOPIR, corkind='pearson', Nfactors = 5, rotate='none', display=TRUE)

VARIMAX(loadings, display=TRUE)
}
}
