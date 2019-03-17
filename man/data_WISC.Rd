\name{data_WISC}
\alias{data_WISC}
\docType{data}
\title{data_WISC}
\description{A data frame with scores for 175 cases on 10 WISC-R subscales, used 
	by Tabacknick & Fidell (2013, p. 737) in their section on confirmatory factor analysis.}
\usage{data(data_WISC)}
\references{Tabachnik, B. G., & Fidell, L. S. (2014). \emph{Using multivariate statistics.} New York, NY: Pearson.}
\examples{
head(data_WISC) 

\dontrun{
# MAP test on the data_WISC data
MAP(data_WISC, corkind='pearson', display=TRUE)

# parallel analysis of the data_WISC data
RAWPAR(data_WISC, extract='PCA', Ndatasets=100, percentile=95,
       corkind='pearson', display=TRUE)
}
}
