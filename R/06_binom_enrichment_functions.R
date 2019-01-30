#' intact.binom.enrich()
#'
#' Provide a rapid way to perform binom.tests on columns of the <.intact> objects created using the lipid.miner function.
#'
#' @param X any column of a <Query.intact> object
#' @param Y any column of a <Universe.intact> object (has to be the same column as X)
#' @param p is a regular pvalue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the intact.binom and the intact.binom.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or q values.
#' @details The universe file db.intact is provided if you are lacking an universe file.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.binom.enrich(queryExample.intact$`Main class`, universeExample.intact$`Main class`)
#'
#' @author Geremy Clair
#' @export
#'
intact.binom.enrich<- function(X,Y,p=0.05, q=1.0){
  #X is the query list
  #Y is the Universe list
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- intact.binom(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}

#' allchains.binom.enrich()
#'
#' Provide a rapid way to perform binom.tests on <.allchains> objects created using the lipid.miner function.
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#' @param p is a regular pvalue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the allchains.binom and the allchains.binom.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or q values.
#' @details the universe file db.allchains is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.binom.enrich(queryExample.allchains, universeExample.allchains)
#'
#' @author Geremy Clair
#' @export
#'
allchains.binom.enrich<- function(X,Y,p=0.05, q=1.0){
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- allchains.binom(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}

#' chain.binom()
#'
#' Provide a rapid way to perform binom.tests on <.chain> objects created using the lipid.miner function.
#'
#' @param X any <Query.chain> object
#' @param Y any <Universe.chain> object
#' @param p is a regular pvalue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the chain.binom and the chain.binom.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or q values.
#' @details the universe file db.chain is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.binom.enrich(queryExample.chain, universeExample.chain)
#'
#' @author Geremy Clair
#' @export
#'
chain.binom.enrich<- function(X,Y,p=0.05, q=1){
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- chain.binom(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}
