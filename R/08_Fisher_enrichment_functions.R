#' intact.fisher.enrich()
#'
#' Provide a rapid way to perform fisher.tests on columns of the <.intact> objects created using the lipid.miner function.
#'
#' @param X any column of a <Query.intact> object
#' @param Y any column of a <Universe.intact> object (has to be the same column as X)
#' @param p is a regular pvalue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the intact.fisher and the intact.fisher.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value.
#' @details The universe file db.intact is provided if you are lacking an universe file.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.fisher.enrich(queryExample.intact$`Main class`, universeExample.intact$`Main class`)
#'
#' @author Geremy Clair
#' @export
#'
intact.fisher.enrich<- function(X,Y,p=0.05, q=1.0){
  #X is the query list
  #Y is the Universe list
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- intact.fisher(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}

#' allchains.fisher()
#'
#' Provide a rapid way to perform fisher.tests on <.allchains> objects created using the lipid.miner function.
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#' @param p is a regular pvalue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the allchains.fisher and the allchains.fisher.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value
#' @details the universe file db.allchains is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.fisher.enrich(queryExample.allchains, universeExample.allchains)
#'
#' @author Geremy Clair
#' @export
#'
allchains.fisher.enrich<-function(X,Y,p=0.05, q=1.0){
  #X is the query list
  #Y is the Universe list
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- allchains.fisher(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}

#' chain.fisher()
#'
#' Provide a rapid way to perform fisher.tests on <.chain> objects created using the lipid.miner function.
#'
#' @param X any <Query.chain> object
#' @param Y any <Universe.chain> object
#' @param p is a regular pvalue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the chain.fisher and the chain.fisher.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value.
#' @details the universe file db.chain is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.fisher.enrich(queryExample.chain, universeExample.chain)
#'
#' @author Geremy Clair
#' @export
#'
chain.fisher.enrich<- function(X,Y,p=0.05, q=1.0){
  #X is the query list
  #Y is the Universe list
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- chain.fisher(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}
