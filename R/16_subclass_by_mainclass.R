#' subclass.mainclass()
#'
#' Function to evaluate the enrichment of subclasses within each main class for the following tests:"Fisher", "Binom", "Hyper", or "EASE". For the KS tests, use the function subclass.mainclass.KS()
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the results for the tests
#'
#' @examples subclass.mainclass(queryExample.intact,universeExample.intact,test= "Fisher", enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

subclass.mainclass<- function(X,Y, test= "Fisher", enrich=FALSE, p = 0.05, q = 1.0){
  if (!exists("test")){test="Fisher"}
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Main class`)

  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")

  for (i in 1:length(unique.query)){
    query.sub <- X[X$`Main class`==unique.query[i],]
    universe.sub <- Y[Y$`Main class`==unique.query[i],]
    sub.name<- unique.query[i]
    if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-intact.fisher(query.sub$`Sub class`, universe.sub$`Sub class`)}
      if (test=="Binom"){result.sub<-intact.binom(query.sub$`Sub class`, universe.sub$`Sub class`)}
      if (test=="Hyper"){result.sub<-intact.hyper(query.sub$`Sub class`, universe.sub$`Sub class`)}
      if (test=="EASE"){result.sub<-intact.EASE(query.sub$`Sub class`, universe.sub$`Sub class`)}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(query.sub$`Sub class`,universe.sub$`Sub class`, p=p, q=q)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(query.sub$`Sub class`, universe.sub$`Sub class`, p=p, q=q)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(query.sub$`Sub class`, universe.sub$`Sub class`, p=p, q=q)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(query.sub$`Sub class`, universe.sub$`Sub class`, p=p, q=q)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste("subclass",result.sub$Classifier,"within the main class", sub.name)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}

#' subclass.mainclass.KS()
#'
#' Function to evaluate the enrichment of subclasses within each main class using a ks.test. For other statistical tests, use subclass.mainclass()
#'
#' @param X has to be a <.intact> object that was created using the lipid.miner function
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a p or q value
#' @param p a maximum pvalue cutoff (set at 0.05 by default)
#' @param q a maximum FDR q-value cutoff (set at 1 by default)
#'
#' @return This function will return a data.frame containing the results for the tests
#'
#' @examples subclass.mainclass.KS(RTexample.intact,cleaned.RTexample, enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

subclass.mainclass.KS<- function(X,rankingTable, order="ascending",enrich=FALSE, p = 0.05, q = 1.0){
  if(!ncol(X)==18 & !sum(colnames(X)[1:3]==c("Lipids","Category","Main class"))==3){stop("X should be an <.intact> object")}
  if(!is.data.frame(rankingTable)){stop("<rankingTable> should be a data frame containing the lipid list as first column and the ranking values as second column, we recommend that you use the function clean.RankingTable() to clean the table")}
  if(ncol(rankingTable)!=2){stop("<rankingTable> should contain two columns")}
  rankingTable[rankingTable==""]<-NA
  if(sum(is.na(rankingTable))>0){stop("<rankingTable> should not contain missing values")}

  if(missing(order)){order=="ascending"}
  if(!order %in% c("ascending", "descending")){stop("<order> should either be 'ascending' or 'descending'")}

  if (missing(enrich)){enrich=FALSE}

  if (missing(p)){p<-0.05}
  if (missing(q)){q<-1.0}
  if (p<0|!is.numeric(p)){stop("<p> should be numeric and positive")}
  if (q<0|!is.numeric(q)){stop("<q> should be numeric and positive")}

  unique.query <- unique(X$`Main class`)
  final.table <- data.frame(matrix(NA, nrow=0, ncol=4))
  colnames(final.table)<- c("Classifier","Count.in.list", "p-value", "FDR.q-value")

  for (i in 1:length(unique.query)){
    query.sub <- subsetmainclass(X, mainclass = unique.query[i])
    sub.name<- unique.query[i]
    rankingTable.sub<-rankingTable[rankingTable[,1] %in% query.sub$Lipid,]
    if(enrich==FALSE){
      result.sub<-intact.KS(query.sub, colname="Sub.class", rankingTable.sub, order=order)
    }else{
      result.sub<-intact.KS.enrich(query.sub, colname="Sub.class", rankingTable, order=order, p=p, q=q)}
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"within the main class",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  return(final.table)
}

