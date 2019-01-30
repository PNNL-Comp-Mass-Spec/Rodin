#' Total.Carbon.cat()
#'
#' Total number of chain carbon enriched within each category
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function total.carbon.cat
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE" for "KS test" please the function total.carbon.cat.KS()
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() p-value, q-value and a fold change
#'
#' @examples total.carbon.cat(queryExample.intact,universeExample.intact,test="Fisher",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.carbon.cat<- function(X,Y, test= "Fisher", enrich=FALSE, p = 0.05, q = 1.0){
  if (!exists("test")){test="Fisher"}
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Category`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")

  for (i in 1:length(unique.query)){
    query.sub <- subsetcat(X, cat = unique.query[i])
    universe.sub <- subsetcat(Y, cat = unique.query[i])
    sub.name<- unique.query[i]
    if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-intact.fisher(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`))}
      if (test=="Binom"){result.sub<-intact.binom(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`))}
      if (test=="Hyper"){result.sub<-intact.hyper(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`))}
      if (test=="EASE"){result.sub<-intact.EASE(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`))}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`), p=p, q=q)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`), p=p, q=q)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`), p=p, q=q)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`), p=p, q=q)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain carbon of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}

#' Total.DB.cat()
#'
#' Total number of chain unsaturation enriched within each category
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE" for "KS test" please the function total.DB.cat.KS()
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.DB.cat(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.DB.cat<- function(X,Y, test= "Fisher", enrich=FALSE, p = 0.05, q = 1.0){
  if (!exists("test")){test="Fisher"}
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Category`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")

  for (i in 1:length(unique.query)){
    query.sub <- subsetcat(X, cat = unique.query[i])
    universe.sub <- subsetcat(Y, cat = unique.query[i])
    sub.name<- unique.query[i]
    if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-intact.fisher(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`))}
      if (test=="Binom"){result.sub<-intact.binom(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`))}
      if (test=="Hyper"){result.sub<-intact.hyper(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`))}
      if (test=="EASE"){result.sub<-intact.EASE(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`))}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`), p=p, q=q)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`), p=p, q=q)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`), p=p, q=q)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`), p=p, q=q)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of unsaturation of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}

#' Total.Carbon.main()
#'
#' Total number of chain carbon enriched within each main lipid class
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE" for "KS test" please the function total.carbon.main.KS
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.carbon.main(queryExample.intact,universeExample.intact,test="Binom",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.carbon.main<- function(X,Y, test= "Fisher", enrich=FALSE, p = 0.05, q = 1.0){
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
  query.sub <- subsetmainclass(X, mainclass = unique.query[i])
  universe.sub <- subsetmainclass(Y, mainclass = unique.query[i])
  sub.name<- unique.query[i]
  if(enrich==FALSE){
  if (test=="Fisher"){result.sub<-intact.fisher(as.character(query.sub$Total.Number.of.Carbon),as.character(universe.sub$Total.Number.of.Carbon))}
  if (test=="Binom"){result.sub<-intact.binom(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon))}
  if (test=="Hyper"){result.sub<-intact.hyper(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon))}
  if (test=="EASE"){result.sub<-intact.EASE(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon))}
  }
  if(enrich==TRUE){
    if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), p=p, q=q)}
    if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), p=p, q=q)}
    if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), p=p, q=q)}
    if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), p=p, q=q)}
  }
  if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain carbon of",result.sub$Classifier)}
  final.table<- rbind(final.table,result.sub)
  }
if (nrow(final.table)>1){final.table <- final.table[-1,]}
return(final.table)
}

#' Total.DB.main()
#'
#' Total number of chain unsaturation enriched within each main lipid class
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE" for "KS test" please the function total.DB.main.KS()
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.DB.main(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.DB.main<- function(X,Y, test= "Fisher", enrich=FALSE, p = 0.05, q = 1.0){
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
    query.sub <- subsetmainclass(X, mainclass = unique.query[i])
    universe.sub <- subsetmainclass(Y, mainclass = unique.query[i])
    sub.name<- unique.query[i]
    if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-intact.fisher(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds))}
      if (test=="Binom"){result.sub<-intact.binom(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds))}
      if (test=="Hyper"){result.sub<-intact.hyper(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds))}
      if (test=="EASE"){result.sub<-intact.EASE(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds))}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), p=p, q=q)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), p=p, q=q)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), p=p, q=q)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), p=p, q=q)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain unsaturation of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}

#' Total.Carbon.sub()
#'
#' Total number of chain carbon enriched within each lipid subclass
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE" for "KS test" please the function total.carbon.sub.KS()
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.carbon.sub(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.carbon.sub<- function(X,Y, test= "Fisher", enrich=FALSE, p = 0.05, q = 1.0){
  if (!exists("test")){test="Fisher"}
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Sub class`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")

  for (i in 1:length(unique.query)){
    query.sub <- subsetsubclass(X, subclass = unique.query[i])
    universe.sub <- subsetsubclass(Y, subclass = unique.query[i])
    sub.name<- unique.query[i]
    if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-intact.fisher(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon))}
      if (test=="Binom"){result.sub<-intact.binom(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon))}
      if (test=="Hyper"){result.sub<-intact.hyper(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon))}
      if (test=="EASE"){result.sub<-intact.EASE(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon))}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), p=p, q=q)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), p=p, q=q)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), p=p, q=q)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), p=p, q=q)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain carbon of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}

#' Total.DB.sub()
#'
#' Total number of chain unsaturation enriched within each lipid subclass
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE" for "KS test" please the function total.DB.sub.KS()
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.DB.sub(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.DB.sub<- function(X,Y, test= "Fisher", enrich=FALSE, p = 0.05, q = 1.0){
  if (!exists("test")){test="Fisher"}
  if (missing(p)){p=0.05}
  if (missing(q)){q=1}
  if (p<0){p=0.05}
  if (q<0){q=1}
  if (!exists("enrich")){enrich=F}
  unique.query <- unique(X$`Sub class`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","p-value", "FDR.q-value","Fold.change")

  for (i in 1:length(unique.query)){
    query.sub <- subsetsubclass(X, subclass = unique.query[i])
    universe.sub <- subsetsubclass(Y, subclass = unique.query[i])
    sub.name<- unique.query[i]
    if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-intact.fisher(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds))}
      if (test=="Binom"){result.sub<-intact.binom(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds))}
      if (test=="Hyper"){result.sub<-intact.hyper(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds))}
      if (test=="EASE"){result.sub<-intact.EASE(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds))}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), p=p, q=q)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), p=p, q=q)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), p=p, q=q)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), p=p, q=q)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain unsaturation of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}

#' Total.Carbon.cat.KS()
#'
#' Total number of chain carbon enriched within each category only for the ks.test, for the other enrhchment tests, use the function Total.Carbon.cat()
#'
#' @param X has to be a <.intact> object that was created using the lipid.miner function
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the results for the tests
#'
#' @examples total.carbon.cat.KS(RTexample.intact,cleaned.RTexample,enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

total.carbon.cat.KS<- function(X,rankingTable, order="ascending", enrich=FALSE, p = 0.05, q = 1.0){
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

  unique.query <- unique(X$`Category`)
  final.table <- data.frame(matrix(NA, nrow=0, ncol=4))
  colnames(final.table)<- c("Classifier","Count.in.list", "p-value", "FDR.q-value")

  for (i in 1:length(unique.query)){
    query.sub <- subsetcat(X, cat = unique.query[i])
    sub.name<- unique.query[i]
    rankingTable.sub<-rankingTable[rankingTable[,1] %in% query.sub$Lipid,]
    if(enrich==FALSE){
      result.sub<-intact.KS(query.sub, colname="Total Number of Carbon", rankingTable.sub, order=order)
    }else{
      result.sub<-intact.KS.enrich(query.sub, colname="Total Number of Carbon", rankingTable, order=order, p=p, q=q)}
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain carbon of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  return(final.table)
}

#' Total.Carbon.main.KS()
#'
#' Total number of chain carbon enriched within each main class only for the ks.test, for the other enrhchment tests, use the function Total.Carbon.main()
#'
#' @param X has to be a <.intact> object that was created using the lipid.miner function
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a p or q value
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the results for the tests
#' @examples total.carbon.main.KS(RTexample.intact,cleaned.RTexample,enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

total.carbon.main.KS<- function(X,rankingTable, order="ascending", enrich=FALSE, p = 0.05, q = 1.0){
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
      result.sub<-intact.KS(query.sub, colname="Total.Number.of.Carbon", rankingTable.sub, order=order)
    }else{
      result.sub<-intact.KS.enrich(query.sub, colname="Total.Number.of.Carbon", rankingTable, order=order, p=p, q=q)}
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain carbon of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  return(final.table)
}

#' Total.Carbon.sub.KS()
#'
#' Total number of chain carbon enriched within each main class only for the ks.test, for the other enrhchment tests, use the function Total.Carbon.sub()
#'
#' @param X has to be a <.intact> object that was created using the lipid.miner function
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a p or q value
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the results for the tests
#' @examples total.carbon.sub.KS(RTexample.intact,cleaned.RTexample,enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

total.carbon.sub.KS<- function(X,rankingTable, order="ascending", enrich=FALSE, p = 0.05, q = 1.0){
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

  unique.query <- unique(X$`Sub class`)
  final.table <- data.frame(matrix(NA, nrow=0, ncol=4))
  colnames(final.table)<- c("Classifier","Count.in.list", "p-value", "FDR.q-value")

  for (i in 1:length(unique.query)){
    query.sub <- subsetsubclass (X, subclass = unique.query[i])
    sub.name<- unique.query[i]
    rankingTable.sub<-rankingTable[rankingTable[,1] %in% query.sub$Lipid,]
    if(enrich==FALSE){
      result.sub<-intact.KS(query.sub, colname="Total.Number.of.Carbon", rankingTable.sub, order=order)
    }else{
      result.sub<-intact.KS.enrich(query.sub, colname="Total.Number.of.Carbon", rankingTable, order=order, p=p, q=q)}
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain carbon of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  return(final.table)
}

#' Total.DB.cat.KS()
#'
#' Total number of unsaturations enriched within each category only for the ks.test, for the other enrhchment tests, use the function Total.DB.cat()
#'
#' @param X has to be a <.intact> object that was created using the lipid.miner function
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the results for the tests
#'
#' @examples total.DB.cat.KS(RTexample.intact,cleaned.RTexample,enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

total.DB.cat.KS<- function(X,rankingTable, order="ascending", enrich=FALSE, p = 0.05, q = 1.0){
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

  unique.query <- unique(X$`Category`)
  final.table <- data.frame(matrix(NA, nrow=0, ncol=4))
  colnames(final.table)<- c("Classifier","Count.in.list", "p-value", "FDR.q-value")

  for (i in 1:length(unique.query)){
    query.sub <- subsetcat(X, cat = unique.query[i])
    sub.name<- unique.query[i]
    rankingTable.sub<-rankingTable[rankingTable[,1] %in% query.sub$Lipid,]
    if(enrich==FALSE){
      result.sub<-intact.KS(query.sub, colname="Double Bonds", rankingTable.sub, order=order)
    }else{
      result.sub<-intact.KS.enrich(query.sub, colname="Double Bonds", rankingTable, order=order, p=p, q=q)}
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of unsaturation of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  return(final.table)
}

#' Total.DB.main.KS()
#'
#' Total number unsaturations enriched within each main class only for the ks.test, for the other enrhchment tests, use the function Total.DB.main()
#'
#' @param X has to be a <.intact> object that was created using the lipid.miner function
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a p or q value
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the results for the tests
#' @examples total.DB.main.KS(RTexample.intact,cleaned.RTexample,enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

total.DB.main.KS<- function(X,rankingTable, order="ascending", enrich=FALSE, p = 0.05, q = 1.0){
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
      result.sub<-intact.KS(query.sub, colname="Double.Bonds", rankingTable.sub, order=order)
    }else{
      result.sub<-intact.KS.enrich(query.sub, colname="Double.Bonds", rankingTable, order=order, p=p, q=q)}
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of unsaturation of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  return(final.table)
}

#' Total.Carbon.sub.KS()
#'
#' Total number of unsaturation enriched within each main class only for the ks.test, for the other enrhchment tests, use the function Total.DB.sub()
#'
#' @param X has to be a <.intact> object that was created using the lipid.miner function
#' @param rankingTable a RankingTable containing the lipid list as first column and the ranking values as second column
#' @param order indicates if the order of the weighting should be done in ascending order (e.g. pvalues) or descending order (e.g. fold change)
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a p or q value
#' @param p regular pvalue cutoff (0.05 by default)
#' @param q is a FDR qvalue adjusted using the bioconductor package qvalue (1 by default)
#'
#' @return This function will return a data.frame containing the results for the tests
#' @examples total.DB.sub.KS(RTexample.intact,cleaned.RTexample,enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

total.DB.sub.KS<- function(X,rankingTable, order="ascending", enrich=FALSE, p = 0.05, q = 1.0){
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

  unique.query <- unique(X$`Sub class`)
  final.table <- data.frame(matrix(NA, nrow=0, ncol=4))
  colnames(final.table)<- c("Classifier","Count.in.list", "p-value", "FDR.q-value")

  for (i in 1:length(unique.query)){
    query.sub <- subsetsubclass (X, subclass = unique.query[i])
    sub.name<- unique.query[i]
    rankingTable.sub<-rankingTable[rankingTable[,1] %in% query.sub$Lipid,]
    if(enrich==FALSE){
      result.sub<-intact.KS(query.sub, colname="Double.Bonds", rankingTable.sub, order=order)
    }else{
      result.sub<-intact.KS.enrich(query.sub, colname="Double.Bonds", rankingTable, order=order, p=p, q=q)}
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of unsaturation of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  return(final.table)
}

