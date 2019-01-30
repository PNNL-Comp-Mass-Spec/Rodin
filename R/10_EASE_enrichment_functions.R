#' intact.EASE.enrich()
#'
#' Provide a rapid way to perform EASE.tests on columns of the <.intact> objects created using the lipid.miner function.
#'
#' @param X any column of a <Query.intact> object
#' @param Y any column of a <Universe.intact> object (has to be the same column as X)
#' @param p is a regular pue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the intact.EASE and the intact.EASE.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value.
#' @details The universe file db.intact is provided if you are lacking an universe file.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.EASE.enrich(queryExample.intact$`Main class`, universeExample.intact$`Main class`)
#'
#' @author Geremy Clair
#' @export
#'
intact.EASE.enrich<- function(X,Y,p=0.05, q=1.0){
  #X is the query list
  #Y is the Universe list
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- intact.EASE(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}

#' allchains.EASE()
#'
#' Provide a rapid way to perform EASE.tests on <.allchains> objects created using the lipid.miner function.
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#' @param p is a regular pue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue
#'
#' @details The difference between the allchains.EASE and the allchains.EASE.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value.
#' @details the universe file db.allchains is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.EASE.enrich(queryExample.allchains, universeExample.allchains)
#'
#' @author Geremy Clair
#' @export
#'
allchains.EASE.enrich<- function(X,Y,p=0.05, q=1.0){
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- allchains.EASE(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}

#' chain.EASE()
#'
#' Provide a rapid way to perform EASE.tests on <.chain> objects created using the lipid.miner function.
#'
#' @param X any <Query.chain> object
#' @param Y any <Universe.chain> object
#' @param p is a regular pue to use as a cutoff for the enrichment
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue

#'
#' @details The difference between the chain.EASE and the chain.EASE.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a p or a q value.
#' @details the universe file db.chain is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.EASE.enrich(queryExample.chain, universeExample.chain)
#'
#' @author Geremy Clair
#' @export
#'
chain.EASE.enrich<- function(X,Y,p=0.05, q=1.0){
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  lipidEnrich <- chain.EASE(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$`p-value`<p,]
  lipidEnrich <- lipidEnrich[lipidEnrich$`FDR.q-value`<q,]
  lipidEnrich <- lipidEnrich[lipidEnrich$Fold.change>1,]
  lipidEnrich
}
