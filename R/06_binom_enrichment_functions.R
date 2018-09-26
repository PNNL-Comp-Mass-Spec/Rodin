#' intact.binom.enrich()
#'
#' Provide a rapid way to perform binom.tests on columns of the <.intact> objects created using the lipid.miner function.
#'
#' @param X any column of a <Query.intact> object
#' @param Y any column of a <Universe.intact> object (has to be the same column as X)
#' @param pval is a regular pvalue to use as a cutoff for the enrichment
#' @param adjpval is a BH adjusted pvalue to use as a cutoff for the enrichment

#'
#' @details The difference between the intact.binom and the intact.binom.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a Pvalue or a BH adjusted pvalue.
#' @details The universe file db.intact is provided if you are lacking an universe file.
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.binom.enrich(queryExample.intact$`Main class`, universeExample.intact$`Main class`)
#'
#' @author Geremy Clair
#' @export
#'
intact.binom.enrich<- function(X,Y,pval=0.05, adjpval=1.0){
  #X is the query list
  #Y is the Universe list
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1.0}
  if (pval<0){pval=0.05}
  if (adjpval<0){adjpval=1.0}
  lipidEnrich <- intact.binom(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$Pvalue<pval,]
  lipidEnrich <- lipidEnrich[lipidEnrich$BHadjustPvalue<adjpval,]
  lipidEnrich <- lipidEnrich[lipidEnrich$fold.change>1,]
  lipidEnrich
}

#' allchains.binom.enrich()
#'
#' Provide a rapid way to perform binom.tests on <.allchains> objects created using the lipid.miner function.
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#' @param pval is a regular pvalue to use as a cutoff for the enrichment
#' @param adjpval is a BH adjusted pvalue to use as a cutoff for the enrichment

#'
#' @details The difference between the allchains.binom and the allchains.binom.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a Pvalue or a BH adjusted pvalue.
#' @details the universe file db.allchains is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.binom.enrich(queryExample.allchains, universeExample.allchains)
#'
#' @author Geremy Clair
#' @export
#'
allchains.binom.enrich<- function(X,Y,pval=0.05, adjpval=1.0){
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1.0}
  if (pval<0){pval=0.05}
  if (adjpval<0){adjpval=1.0}
  lipidEnrich <- allchains.binom(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$Pvalue<pval,]
  lipidEnrich <- lipidEnrich[lipidEnrich$BHadjustPvalue<adjpval,]
  lipidEnrich <- lipidEnrich[lipidEnrich$fold.change>1,]
  lipidEnrich
}

#' chain.binom()
#'
#' Provide a rapid way to perform binom.tests on <.allchains> objects created using the lipid.miner function.
#'
#' @param X any <Query.allchains> object
#' @param Y any <Universe.allchains> object
#' @param pval is a regular pvalue to use as a cutoff for the enrichment
#' @param adjpval is a BH adjusted pvalue to use as a cutoff for the enrichment

#'
#' @details The difference between the chain.binom and the chain.binom.enrich function is that the <.enrich> function allows to extract only the element that are enriched according to a Pvalue or a BH adjusted pvalue.
#' @details the universe file db.allchains is provided if you are lacking an universe file
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.binom.enrich(queryExample.chain, universeExample.chain)
#'
#' @author Geremy Clair
#' @export
#'
chain.binom.enrich<- function(X,Y,pval=0.05, adjpval=1.0){
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1.0}
  if (pval<0){pval=0.05}
  if (adjpval<0){adjpval=1.0}
  lipidEnrich <- chain.binom(X,Y)
  lipidEnrich <- lipidEnrich[lipidEnrich$Pvalue<pval,]
  lipidEnrich <- lipidEnrich[lipidEnrich$BHadjustPvalue<adjpval,]
  lipidEnrich <- lipidEnrich[lipidEnrich$fold.change>1,]
  lipidEnrich
}
