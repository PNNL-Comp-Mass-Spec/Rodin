#' subclass.mainclass()
#'
#' Function to evaluate the enrichment of subclasses within each main class
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param pval regular pvalue cutoff (0.05 by default)
#' @param adjpval BH adjusted pvalue cutoff (1.0 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples subclass.mainclass(queryExample.intact,universeExample.intact,test= "Fisher", enrich=FALSE)
#'
#' @author Geremy Clair
#' @export

subclass.mainclass<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Main class`)

  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

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
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(query.sub$`Sub class`,universe.sub$`Sub class`, pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(query.sub$`Sub class`, universe.sub$`Sub class`, pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(query.sub$`Sub class`, universe.sub$`Sub class`, pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(query.sub$`Sub class`, universe.sub$`Sub class`, pval=pval, adjpval=adjpval)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste("subclass",result.sub$Classifier,"within the main class", sub.name)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}
