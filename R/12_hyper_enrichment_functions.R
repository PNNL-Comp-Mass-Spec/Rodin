#' intact.hyper.enrich()
#'
#' Provide a rapid way to perform hyper.tests on columns of the <.intact> objects created using the lipid.miner function.
#'
#' @param X any column of a <Query.intact> object
#' @param Y any column of a <Universe.intact> object (has to be the same column as X)
#' @param p is a regular pue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the intact.hyper and the intact.hyper.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value.
#' @details The universe file db.intact is provided if you are lacking an universe file.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.hyper.enrich(queryExample.intact$`Main class`, universeExample.intact$`Main class`)
#'
#' @author Geremy Clair
#' @export
#'
intact.hyper.enrich<- function(X,Y,p=0.05, q=1.0){
  #X is the query list
  #Y is the Universe list
  if (!exists("p")){p=0.05}
  if (!exists("q")){q=1.0}
  if (p<0){p=0.05}
  if (a<0){q=1.0}
  lipidEnrich <- intact.hyper(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}
#' allchains.hyper()
#'
#' Provide a rapid way to perform hyper.tests on <.allchains> objects created using the lipid.miner function.
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#' @param p is a regular pue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the allchains.hyper and the allchains.hyper.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value.
#' @details the universe file db.allchains is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.hyper.enrich(queryExample.allchains, universeExample.allchains)
#'
#' @author Geremy Clair
#' @export
#'
allchains.hyper.enrich<- function(X,Y,p=0.05, q=1.0){
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- allchains.hyper(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}

#' chain.hyper()
#'
#' Provide a rapid way to perform hyper.tests on <.chain> objects created using the lipid.miner function.
#'
#' @param X any <Query.chain> object
#' @param Y any <Universe.chain> object
#' @param p is a regular pue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the chain.hyper and the chain.hyper.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.hyper.enrich(queryExample.chain, universeExample.chain)
#'
#' @author Geremy Clair
#' @export
#'
chain.hyper.enrich<- function(X,Y,p=0.05, q=1.0){
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- chain.hyper(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<a,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}
