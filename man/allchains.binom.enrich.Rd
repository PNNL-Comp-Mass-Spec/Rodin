% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_binom_enrichment_functions.R
\name{allchains.binom.enrich}
\alias{allchains.binom.enrich}
\title{allchains.binom.enrich()}
\usage{
allchains.binom.enrich(X, Y, pval = 0.05, adjpval = 1)
}
\arguments{
\item{X}{any <Query.allchains> object}

\item{Y}{any <Universe.allchains> object}

\item{pval}{is a regular pvalue to use as a cutoff for the enrichment}

\item{adjpval}{is a BH adjusted pvalue to use as a cutoff for the enrichment}
}
\value{
This function will return a data frame with the results of the tests
}
\description{
Provide a rapid way to perform binom.tests on <.allchains> objects created using the lipid.miner function.
}
\details{
The difference between the allchains.binom and the allchains.binom.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a Pvalue or a BH adjusted pvalue.

the universe file db.allchains is provided if you are lacking an universe file
}
\examples{
allchains.binom.enrich(queryExample.allchains, universeExample.allchains)

}
\author{
Geremy Clair
}
