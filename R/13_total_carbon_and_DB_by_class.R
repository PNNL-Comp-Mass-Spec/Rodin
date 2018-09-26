#' Total.Carbon.cat()
#'
#' Total number of chain carbon enriched within each category
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
#' @examples total.carbon.cat(queryExample.intact,universeExample.intact,test="Fisher",enrich=FALSE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.carbon.cat<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Category`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

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
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`), pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`), pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`), pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$`Total Number of Carbon`), as.character(universe.sub$`Total Number of Carbon`), pval=pval, adjpval=adjpval)}
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
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param pval regular pvalue cutoff (0.05 by default)
#' @param adjpval BH adjusted pvalue cutoff (1.0 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.DB.cat(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.DB.cat<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Category`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

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
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`), pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`), pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`), pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$`Double Bonds`), as.character(universe.sub$`Double Bonds`), pval=pval, adjpval=adjpval)}
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
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param pval regular pvalue cutoff (0.05 by default)
#' @param adjpval BH adjusted pvalue cutoff (1.0 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.carbon.main(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.carbon.main<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}

unique.query <- unique(X$`Main class`)
final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

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
    if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), pval=pval, adjpval=adjpval)}
    if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), pval=pval, adjpval=adjpval)}
    if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), pval=pval, adjpval=adjpval)}
    if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), pval=pval, adjpval=adjpval)}
  }
  if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain carbon of",result.sub$Classifier)}
  final.table<- rbind(final.table,result.sub)
  }
if (nrow(final.table)>1){final.table <- final.table[-1,]}
return(final.table)
}

#' Total.Carbon.main()
#'
#' Total number of chain unsaturation enriched within each main lipid class
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
#' @examples total.DB.main(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.DB.main<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Main class`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

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
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), pval=pval, adjpval=adjpval)}
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
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param pval regular pvalue cutoff (0.05 by default)
#' @param adjpval BH adjusted pvalue cutoff (1.0 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.carbon.sub(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.carbon.sub<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}

  unique.query <- unique(X$`Sub class`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

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
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$Total.Number.of.Carbon), as.character(universe.sub$Total.Number.of.Carbon), pval=pval, adjpval=adjpval)}
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
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param pval regular pvalue cutoff (0.05 by default)
#' @param adjpval BH adjusted pvalue cutoff (1.0 by default)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples total.DB.sub(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

total.DB.sub<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}
  unique.query <- unique(X$`Sub class`)
  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

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
      if (test=="Fisher"){result.sub<-intact.fisher.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-intact.binom.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-intact.hyper.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-intact.EASE.enrich(as.character(query.sub$Double.Bonds), as.character(universe.sub$Double.Bonds), pval=pval, adjpval=adjpval)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with a total number of chain unsaturation of",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}
