#' allchains.cat()
#'
#' Function to evaluate the specific chains enriched within each category
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param pval regular pvalue cutoff (0.05 by default)
#' @param adjpval BH adjusted pvalue cutoff (1.0 by default)
#' @param TGcollapse.rm Boolean to indicate if the collapsed TG need to be removed (default TRUE)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples allchains.cat(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

allchains.cat<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0, TGcollapse.rm=TRUE){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}
  if(!exists("TGcollapse.rm")){TGcollapse.rm<-TRUE}

  universe.Allchains.cat<- data.frame(matrix(nrow = 4*nrow(Y), ncol=3))
  universe.Allchains.cat[,1]<- Y$Lipid
  universe.Allchains.cat[,2]<- Y$Category
  universe.Allchains.cat[,3]<-c(as.character(Y$`Chain 1`),as.character(Y$`Chain 2`),as.character(Y$`Chain 3`),as.character(Y$`Chain 4`))
  colnames(universe.Allchains.cat)<- c("Lipid","Category", "Chain")
  universe.Allchains.cat<- universe.Allchains.cat[!is.na(universe.Allchains.cat$Chain),]

  query.Allchains.cat<- data.frame(matrix(nrow = 4*nrow(X), ncol=3))
  query.Allchains.cat[,1]<- X$Lipid
  query.Allchains.cat[,2]<- X$Category
  query.Allchains.cat[,3]<-c(as.character(X$`Chain 1`),as.character(X$`Chain 2`),as.character(X$`Chain 3`),as.character(X$`Chain 4`))
  colnames(query.Allchains.cat)<- c("Lipid","Category", "Chain")
  query.Allchains.cat<- query.Allchains.cat[!is.na(query.Allchains.cat$Chain),]

  if(TGcollapse.rm==TRUE){
    testTGrm<-    !(grepl("TG", query.Allchains.cat$Lipid)&nchar(gsub("\\;.*","",query.Allchains.cat$Lipid))-nchar(gsub("\\/","",gsub("\\;.*","",query.Allchains.cat$Lipid)))!=2)
    query.Allchains.cat<-query.Allchains.cat[testTGrm,]
    testTGrm<-    !(grepl("TG", universe.Allchains.cat$Lipid)&nchar(gsub("\\;.*","",universe.Allchains.cat$Lipid))-nchar(gsub("\\/","",gsub("\\;.*","",universe.Allchains.cat$Lipid)))!=2)
    universe.Allchains.cat<-universe.Allchains.cat[testTGrm,]
  }

  unique.query <- unique(query.Allchains.cat$`Category`)


  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

  for (i in 1:length(unique.query)){
  query.sub <- query.Allchains.cat[query.Allchains.cat$Category==unique.query[i],c(1,3)]
  universe.sub <- universe.Allchains.cat[universe.Allchains.cat$Category==unique.query[i],c(1,3)]
  sub.name<- unique.query[i]
  if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-allchains.fisher(query.sub, universe.sub)}
      if (test=="Binom"){result.sub<-allchains.binom(query.sub, universe.sub)}
      if (test=="Hyper"){result.sub<-allchains.hyper(query.sub, universe.sub)}
      if (test=="EASE"){result.sub<-allchains.EASE(query.sub, universe.sub)}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-allchains.fisher.enrich(query.sub,universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-allchains.binom.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-allchains.hyper.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-allchains.EASE.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with the chain",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}

#' allchains.main()
#'
#' Function to evaluate the specific chains enriched within each main class
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param pval regular pvalue cutoff (0.05 by default)
#' @param adjpval BH adjusted pvalue cutoff (1.0 by default)
#' @param TGcollapse.rm Boolean to indicate if the collapsed TG need to be removed (default TRUE)
#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples allchains.main(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export

allchains.main<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0,TGcollapse.rm=TRUE){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}
  if(!exists("TGcollapse.rm")){TGcollapse.rm<-TRUE}

  universe.Allchains.main<- data.frame(matrix(nrow = 4*nrow(Y), ncol=3))
  universe.Allchains.main[,1]<- Y$Lipid
  universe.Allchains.main[,2]<- Y$`Main class`
  universe.Allchains.main[,3]<-c(as.character(Y$`Chain 1`),as.character(Y$`Chain 2`),as.character(Y$`Chain 3`),as.character(Y$`Chain 4`))
  colnames(universe.Allchains.main)<- c("Lipid","Main class", "Chain")
  universe.Allchains.main<- universe.Allchains.main[!is.na(universe.Allchains.main$Chain),]

  query.Allchains.main<- data.frame(matrix(nrow = 4*nrow(X), ncol=3))
  query.Allchains.main[,1]<- X$Lipid
  query.Allchains.main[,2]<- X$`Main class`
  query.Allchains.main[,3]<-c(as.character(X$`Chain 1`),as.character(X$`Chain 2`),as.character(X$`Chain 3`),as.character(X$`Chain 4`))
  colnames(query.Allchains.main)<- c("Lipid","Main class", "Chain")
  query.Allchains.main<- query.Allchains.main[!is.na(query.Allchains.main$Chain),]

  if(TGcollapse.rm==TRUE){
    testTGrm<-    !(grepl("TG", query.Allchains.main$Lipid)&nchar(gsub("\\;.*","",query.Allchains.main$Lipid))-nchar(gsub("\\/","",gsub("\\;.*","",query.Allchains.main$Lipid)))!=2)
    query.Allchains.main<-query.Allchains.main[testTGrm,]
    testTGrm<-    !(grepl("TG", universe.Allchains.main$Lipid)&nchar(gsub("\\;.*","",universe.Allchains.main$Lipid))-nchar(gsub("\\/","",gsub("\\;.*","",universe.Allchains.main$Lipid)))!=2)
    universe.Allchains.main<-universe.Allchains.main[testTGrm,]
    }

  unique.query <- unique(query.Allchains.main$`Main class`)

  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

  for (i in 1:length(unique.query)){
    query.sub <- query.Allchains.main[query.Allchains.main$`Main class`==unique.query[i],c(1,3)]
    universe.sub <- universe.Allchains.main[universe.Allchains.main$`Main class`==unique.query[i],c(1,3)]
    sub.name<- unique.query[i]
    if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-allchains.fisher(query.sub, universe.sub)}
      if (test=="Binom"){result.sub<-allchains.binom(query.sub, universe.sub)}
      if (test=="Hyper"){result.sub<-allchains.hyper(query.sub, universe.sub)}
      if (test=="EASE"){result.sub<-allchains.EASE(query.sub, universe.sub)}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-allchains.fisher.enrich(query.sub,universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-allchains.binom.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-allchains.hyper.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-allchains.EASE.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with the chain",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}

#' allchains.sub()
#'
#' Function to evaluate the specific chains enriched within each sub class
#'
#' @param X has to be a <Query>.intact object that was created using the lipid.miner function
#' @param Y has to be your <Universe>.intact object that was created using the lipid.miner function
#' @param test has to be a type of test in the list "Fisher", "Binom", "Hyper", or "EASE"
#' @param enrich TRUE or FALSE (FALSE by default) to choose between showing all the results of the test or only the enriched ones based on a pvalue
#' @param pval regular pvalue cutoff (0.05 by default)
#' @param adjpval BH adjusted pvalue cutoff (1.0 by default)
#' @param TGcollapse.rm Boolean to indicate if the collapsed TG need to be removed (default TRUE)

#'
#' @return This function will return a data.frame containing the tested classifiers (rows), the counts and percentage of identifiers/chains falling in this classifier, a binom.test() pvalue, a BH corrected pvalue and a fold change
#'
#' @examples allchains.sub(queryExample.intact,universeExample.intact,test="Hyper",enrich=TRUE, p=0.05)
#'
#' @author Geremy Clair
#' @export
#'

allchains.sub<- function(X,Y, test= "Fisher", enrich=FALSE, pval = 0.05, adjpval = 1.0,TGcollapse.rm=TRUE){
  if (!exists("test")){test="Fisher"}
  if (!exists("pval")){pval=0.05}
  if (!exists("adjpval")){adjpval=1}
  if (!exists("enrich")){enrich=F}
  if(!exists("TGcollapse.rm")){TGcollapse.rm<-TRUE}

  universe.Allchains.sub<- data.frame(matrix(nrow = 4*nrow(Y), ncol=3))
  universe.Allchains.sub[,1]<- Y$Lipid
  universe.Allchains.sub[,2]<- Y$`Sub class`
  universe.Allchains.sub[,3]<-c(as.character(Y$`Chain 1`),as.character(Y$`Chain 2`),as.character(Y$`Chain 3`),as.character(Y$`Chain 4`))
  colnames(universe.Allchains.sub)<- c("Lipid","Sub class", "Chain")
  universe.Allchains.sub<- universe.Allchains.sub[!is.na(universe.Allchains.sub$Chain),]

  query.Allchains.sub<- data.frame(matrix(nrow = 4*nrow(X), ncol=3))
  query.Allchains.sub[,1]<- X$Lipid
  query.Allchains.sub[,2]<- X$`Sub class`
  query.Allchains.sub[,3]<-c(as.character(X$`Chain 1`),as.character(X$`Chain 2`),as.character(X$`Chain 3`),as.character(X$`Chain 4`))
  colnames(query.Allchains.sub)<- c("Lipid","Sub class", "Chain")
  query.Allchains.sub<- query.Allchains.sub[!is.na(query.Allchains.sub$Chain),]

  if(TGcollapse.rm==TRUE){
    testTGrm<-    !(grepl("TG", query.Allchains.sub$Lipid)&nchar(gsub("\\;.*","",query.Allchains.sub$Lipid))-nchar(gsub("\\/","",gsub("\\;.*","",query.Allchains.sub$Lipid)))!=2)
    query.Allchains.sub<-query.Allchains.sub[testTGrm,]
    testTGrm<-    !(grepl("TG", universe.Allchains.sub$Lipid)&nchar(gsub("\\;.*","",universe.Allchains.sub$Lipid))-nchar(gsub("\\/","",gsub("\\;.*","",universe.Allchains.sub$Lipid)))!=2)
    universe.Allchains.sub<-universe.Allchains.sub[testTGrm,]
  }

  unique.query <- unique(query.Allchains.sub$`Sub class`)

  final.table <- data.frame(matrix(NA, nrow=length(0), ncol=8))
  colnames(final.table)<- c("Classifier","Count.query","Count.universe","%.query","%.universe","Pvalue","BHadjustPvalue","fold.change")

  for (i in 1:length(unique.query)){
    query.sub <- query.Allchains.sub[query.Allchains.sub$`Sub class`==unique.query[i],c(1,3)]
    universe.sub <- universe.Allchains.sub[universe.Allchains.sub$`Sub class`==unique.query[i],c(1,3)]
    sub.name<- unique.query[i]
    if(enrich==FALSE){
      if (test=="Fisher"){result.sub<-allchains.fisher(query.sub, universe.sub)}
      if (test=="Binom"){result.sub<-allchains.binom(query.sub, universe.sub)}
      if (test=="Hyper"){result.sub<-allchains.hyper(query.sub, universe.sub)}
      if (test=="EASE"){result.sub<-allchains.EASE(query.sub, universe.sub)}
    }
    if(enrich==TRUE){
      if (test=="Fisher"){result.sub<-allchains.fisher.enrich(query.sub,universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="Binom"){result.sub<-allchains.binom.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="Hyper"){result.sub<-allchains.hyper.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
      if (test=="EASE"){result.sub<-allchains.EASE.enrich(query.sub, universe.sub, pval=pval, adjpval=adjpval)}
    }
    if(nrow(result.sub)>0){result.sub$Classifier<- paste(sub.name,"with the chain",result.sub$Classifier)}
    final.table<- rbind(final.table,result.sub)
  }
  if (nrow(final.table)>1){final.table <- final.table[-1,]}
  return(final.table)
}
