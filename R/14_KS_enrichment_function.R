#' intact.KS.enrich()
#'
#' Performs weighted Kolmogorov-Smirnov test on a column of an <.intact> object based on Ranking values contained in a rankingTable. After the test is performed,  the results are filtered based on the ks.test() p or by the FDR q as specified by the user.
#'
#' @param X a <Rodin.intact> object
#' @param colname a colname from the <Rodin.intact> object
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param p a maximum pvalue cutoff (set at 0.05 by default)
#' @param q a maximum FDR q-value cutoff (set at 1 by default)
#'
#' @details Peforms a ks test enrichment test on a column of a <Query.intact> based on ranking values contained in a ranking table
#' @details colname has to be in :"Lipid", "Category","Main class","Sub class","Total Number of Carbon","Double Bonds","Chain1","Chain2","Chain3","Chain4","NCarbonChain1","NCarbonChain2","NCarbonChain3","NCarbonChain4","DBChain1","DBChain2","DBChain3","DBChain4"
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples intact.KS.enrich(RTexample.intact, colname="Main class", cleaned.RTexample, p=0.05)
#'
#' @author Geremy Clair
#' @export

intact.KS.enrich<- function(X,colname="Main class", rankingTable, order="ascending", p=0.05, q=1.0){
  if(!ncol(X)==18 & !sum(colnames(X)[1:3]==c("Lipids","Category","Main class"))==3){stop("X should be an <.intact> object")}
  if(!is.data.frame(rankingTable)){stop("rankingTable should be a data frame of type <Ranking Table> containing the lipid list as first column and the ranking values as second column, we recommend that you use the function clean.RankingTable() to clean the table")}
  if(ncol(rankingTable)!=2){stop("rankingTable should contain two columns")}
  rankingTable[rankingTable==""]<-NA
  if(sum(is.na(rankingTable))>0){stop("rankingTable should not contain missing values")}

  if(!colname %in% colnames(X)){
    print("colname is not present in the object X, the columns available are the folling:")
    print(colnames(X))
    stop()
  }

  if(missing(order)){order=="ascending"}
  if(!order %in% c("ascending", "descending")){stop("order should either be 'ascending' or 'descending'")}
  if (missing(p)){p<-0.05}
  if (missing(q)){q<-1.0}
  if (p<0|!is.numeric(p)){stop("p should be numeric and positive")}
  if (q<0|!is.numeric(q)){stop("q should be numeric and positive")}

  final_table<- intact.KS(X,colname=colname,rankingTable,order=order)
  final_table<- final_table[final_table$`p-value`<p & final_table$`FDR.q-value`<q,]
  return(final_table)

}

#' chain.KS.enrich()
#'
#' Performs weighted Kolmogorov-Smirnov test on a column of an <.chain> object based on Ranking values contained in a rankingTable. After the test is performed,  the results are filtered based on the ks.test() p or by the FDR q as specified by the user.
#'
#' @param X a <Rodin.chain> object
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param p a maximum pvalue cutoff (set at 0.05 by default)
#' @param q a maximum FDR q-value cutoff (set at 1 by default)
#'
#' @details Peforms a ks test enrichment test on a column of a <Query.intact> based on ranking values contained in a ranking table
#' @details colname has to be in :"Lipid", "Category","Main class","Sub class","Total Number of Carbon","Double Bonds","Chain1","Chain2","Chain3","Chain4","NCarbonChain1","NCarbonChain2","NCarbonChain3","NCarbonChain4","DBChain1","DBChain2","DBChain3","DBChain4"
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples chain.KS.enrich(RTexample.chain, cleaned.RTexample, p=0.05,q=1)
#'
#' @author Geremy Clair
#' @export

chain.KS.enrich<- function(X, rankingTable, order="ascending", p=0.05, q=1){
  if(!ncol(X)==9 & !sum(colnames(X)[1:3]==c("Lipids","SCFA","MCFA"))==3){stop("X should be an <.chain> object")}

  if(!is.data.frame(rankingTable)){stop("rankingTable should be a data frame of type <Ranking Table> containing the lipid list as first column and the ranking values as second column, we recommend that you use the function clean.RankingTable() to clean the table")}
  if(ncol(rankingTable)!=2){stop("rankingTable should contain two columns")}
  rankingTable[rankingTable==""]<-NA
  if(sum(is.na(rankingTable))>0){stop("rankingTable should not contain missing values")}

  if(missing(order)){order=="ascending"}
  if(!order %in% c("ascending", "descending")){stop("order should either be 'ascending' or 'descending'")}
  if (missing(p)){p<-0.05}
  if (missing(q)){q<-1}
  if (p<0|!is.numeric(p)){stop("p should be numeric and positive")}
  if (q<0|!is.numeric(q)){stop("q should be numeric and positive")}

  final_table<- chain.KS(X,rankingTable,order=order)
  final_table<- final_table[final_table$`p-value`<p & final_table$`FDR.q-value`<q,]
  return(final_table)

}

#' chain.KS.enrich()
#'
#' Performs weighted Kolmogorov-Smirnov test on a column of an <.chain> object based on Ranking values contained in a rankingTable. After the test is performed,  the results are filtered based on the ks.test() p or by the FDR q as specified by the user.
#'
#' @param X a <Rodin.chain> object
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param p a maximum pvalue cutoff (set at 0.05 by default)
#' @param q a maximum FDR q-value cutoff (set at 1 by default)
#'
#' @details Peforms a ks test enrichment test on a column of a <Query.intact> based on ranking values contained in a ranking table
#' @details colname has to be in :"Lipid", "Category","Main class","Sub class","Total Number of Carbon","Double Bonds","Chain1","Chain2","Chain3","Chain4","NCarbonChain1","NCarbonChain2","NCarbonChain3","NCarbonChain4","DBChain1","DBChain2","DBChain3","DBChain4"
#'
#' @return This function will return a data frame with the results of the tests
#'
#' @examples allchains.KS.enrich(RTexample.allchains, cleaned.RTexample, p=0.05,q=1)
#'
#' @author Geremy Clair
#' @export

allchains.KS.enrich<- function(X, rankingTable, order="ascending", p=0.05, q=1){
  if(!ncol(X)==2 & !sum(colnames(X)[1:2]==c("Lipids","Chain"))==2){stop("X should be an <.allchains> object")}

  if(!is.data.frame(rankingTable)){stop("rankingTable should be a data frame of type <Ranking Table> containing the lipid list as first column and the ranking values as second column, we recommend that you use the function clean.RankingTable() to clean the table")}
  if(ncol(rankingTable)!=2){stop("rankingTable should contain two columns")}
  rankingTable[rankingTable==""]<-NA
  if(sum(is.na(rankingTable))>0){stop("rankingTable should not contain missing values")}

  if(missing(order)){order=="ascending"}
  if(!order %in% c("ascending", "descending")){stop("order should either be 'ascending' or 'descending'")}
  if (missing(p)){p<-0.05}
  if (missing(q)){q<-1}
  if (p<0|!is.numeric(p)){stop("p should be numeric and positive")}
  if (q<0|!is.numeric(q)){stop("q should be numeric and positive")}

  final_table<- allchains.KS(X,rankingTable,order=order)
  final_table<- final_table[final_table$`p-value`<p & final_table$`FDR.q-value`<q,]
  return(final_table)

}
