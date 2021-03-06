% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_subset_lipidclass.R
\name{subsetmainclass}
\alias{subsetmainclass}
\title{subsetmainclass()}
\usage{
subsetmainclass(X, mainclass = "MainClassname")
}
\arguments{
\item{X}{object of class *.intact, *.chain, or *.allchains}

\item{mainclass}{character string specifying the main class of lipid to subset}
}
\value{
An object of the same class as X.
}
\description{
Subset of the lipid chain object based on the main class
}
\details{
Subset an object based on its lipid main class

(text prior the parenthesis in the lipid name).

The main classes supported are in the main.colors object.
}
\examples{
subsetmainclass(queryExample.intact,mainclass="TG")

}
\author{
Geremy Clair
}
