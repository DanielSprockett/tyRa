% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeling.R
\name{fit_sncm}
\alias{fit_sncm}
\title{Fits \code{sncm} to an OTU table}
\usage{
fit_sncm(spp = spp, pool = NULL, taxon = NULL)
}
\arguments{
\item{spp}{A community table for communities of interest with local communities/samples as rows and taxa as columns.
All samples must be rarefied to the same depth.}

\item{pool}{(optional) A community table for defining source community.}

\item{taxon}{(optional) A table listing the taxonomic calls for each otu, with OTU ids as row names
and taxonomic classifications as columns.}
}
\value{
This function returns list of two elements.
The first element, spp.out$fitstats, contains various fitting stats.
The second element contains the predicted occurrence frequencies for each OTU/ASV, as well as their \code{fit_class}
}
\description{
Fits the neutral model from
\href{http://onlinelibrary.wiley.com/doi/10.1111/j.1462-2920.2005.00956.x/abstract}{Sloan \emph{et al.} (2006)}
to an OTU table and returns several fitting statistics as well as predicted occurrence frequencies for each OTU
based on their abundance in the metacommunity. The author of this function is Adam Burns (\email{aburns2@uoregon.edu}),
and was originally published in \href{https://www.nature.com/articles/ismej2015142}{Burns \emph{et al.} (2016)}.
}
\seealso{
\code{\link{plot_sncm_fit}}
}
