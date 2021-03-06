% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_binom_enrichment_functions.R
\name{intact.binom.enrich}
\alias{intact.binom.enrich}
\title{intact.binom.enrich()}
\usage{
intact.binom.enrich(X, Y, pval = 0.05, adjpval = 1)
}
\arguments{
\item{X}{any column of a <Query.intact> object}

\item{Y}{any column of a <Universe.intact> object (has to be the same column as X)}

\item{pval}{is a regular pvalue to use as a cutoff for the enrichment}

\item{adjpval}{is a BH adjusted pvalue to use as a cutoff for the enrichment}
}
\value{
This function will return a data frame with the results of the tests
}
\description{
Provide a rapid way to perform binom.tests on columns of the <.intact> objects created using the lipid.miner function.
}
\details{
The difference between the intact.binom and the intact.binom.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a Pvalue or a BH adjusted pvalue.

The universe file db.intact is provided if you are lacking an universe file.
}
\examples{
intact.binom.enrich(queryExample.intact$`Main class`, universeExample.intact$`Main class`)

}
\author{
Geremy Clair
}
